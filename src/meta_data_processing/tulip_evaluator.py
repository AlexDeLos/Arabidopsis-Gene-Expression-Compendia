"""
tulip_evaluator.py
------------------
Standalone TULIP-based evaluator for pipeline-assigned GEO sample labels.

What it does
------------
For each labelled sample it calls the TULIP LLM with:
  - the raw GEO metadata fields (title, characteristics_ch1, source_name, etc.)
  - the label axes and canonical options
  - the pipeline's assigned labels

TULIP returns a structured JSON verdict per label axis:
  {
    "treatment": {"verdict": "correct"|"incorrect"|"uncertain",
                  "suggestion": "<canonical value or null>",
                  "reason": "<one sentence>"},
    "tissue":    { ... },
    ...
  }

Verdicts are aggregated into:
  - per-sample results  (saved alongside the label files)
  - per-study summary   (accuracy per axis, counts)
  - global summary      (across all evaluated studies)

Usage
-----
From the command line:

    python tulip_evaluator.py \\
        --labels-path  /path/to/labels \\
        --raw-data-path /path/to/processed_microarray_data \\
        --output-path  /path/to/eval_results \\
        [--studies GSE5622 GSE9415] \\
        [--max-samples 20]           # per study, for cost control

Or import and call programmatically:

    from tulip_evaluator import TulipEvaluator
    evaluator = TulipEvaluator(labels_path=..., raw_data_path=..., output_path=...)
    evaluator.run(studies=["GSE5622"])

Design principles
-----------------
* Zero coupling to the rest of the pipeline — only standard-library + openai.
  No imports from constants, groundingOptimizer, universal_extractor, etc.
* The canonical options list is derived at runtime from the label files themselves
  (unique values seen across all loaded labels), so it stays in sync with whatever
  ontology version was used when the labels were produced.
* One TULIP call per sample, returning verdicts for ALL label axes together.
  This is cheaper and gives the model cross-axis context (e.g. it can notice that
  the tissue says "root" but the source_name says "leaf" and flag the inconsistency).
* Results are written incrementally so a crashed run can be resumed.
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import re
from collections import defaultdict

from openai import OpenAI

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)

# ── TULIP constants ────────────────────────────────────────────────────────────
TULIP_BASE_URL = "https://api.tulip.tudelft.nl/chat/v1"
TULIP_CHAT_MODEL = "chat"
TULIP_API_KEY = os.environ.get("TULIP_API_KEY", "DUMMY_API_KEY")
TULIP_MAX_TOKENS = 2048  # Must be large enough for reasoning trace + JSON output.
# The TULIP chat model is a reasoning model — it writes a
# chain-of-thought into reasoning_content before producing
# the final answer in content. 512 was too small; the model
# exhausted the budget on thinking and left content empty.

# ── Prompts ────────────────────────────────────────────────────────────────────
_SYSTEM_PROMPT = """\
You are an expert curator of Arabidopsis thaliana gene-expression metadata from NCBI GEO.

You will be given:
1. Raw GEO metadata fields for one sample (title, characteristics, source name, etc.)
2. The canonical label options available for each label axis
3. The labels that an automated pipeline has already assigned to this sample

Your task is to evaluate whether each assigned label is correct given the raw metadata.

Key rules:
- "unspecified" is a valid pipeline output meaning the metadata did not contain enough
  information to assign a label. Only mark it "incorrect" if the metadata CLEARLY and
  EXPLICITLY states a value for that axis. If it is merely implied or ambiguous, verdict
  should be "uncertain", not "incorrect".
- "suggestion" must be EXACTLY one of the canonical values listed for that axis.
  Never suggest a value that is not in the canonical options list. If the correct value
  is not in the list, use verdict "uncertain" instead.
- Only use "incorrect" when you are confident the metadata contradicts the assigned label.
- Use "uncertain" when the metadata is ambiguous, missing, or when you are not sure.

IMPORTANT: Write the JSON object FIRST, before any reasoning. This ensures a complete
response even if output is truncated. No markdown fences around the JSON.

Format:
{
  "<label_axis>": {"verdict": "correct"|"incorrect"|"uncertain",
                   "suggestion": "<exact canonical value or null>",
                   "reason": "<one sentence>"},
  ...
}
"""

_USER_PROMPT_TEMPLATE = """\
=== RAW METADATA ===
{metadata_block}

=== CANONICAL OPTIONS ===
{options_block}

=== PIPELINE-ASSIGNED LABELS ===
{assigned_block}

Evaluate each assigned label and return your verdict as a JSON object.
"""

# Metadata fields to pull from the raw sample JSON, in priority order.
# We keep this list lean — enough context for the LLM without hitting token limits.
_METADATA_FIELDS = [
    "title",
    "source_name_ch1",
    "characteristics_ch1",
    "description",
    "molecule_ch1",
    "data_processing",
    "extract_protocol_ch1",
]
_STUDY_FIELDS = ["summary", "overall_design"]

# Verdict constants
VERDICT_CORRECT = "correct"
VERDICT_INCORRECT = "incorrect"
VERDICT_UNCERTAIN = "uncertain"
VALID_VERDICTS = {VERDICT_CORRECT, VERDICT_INCORRECT, VERDICT_UNCERTAIN}


# ── TulipEvaluator ─────────────────────────────────────────────────────────────


class TulipEvaluator:
    """
    Evaluates pipeline-assigned GEO sample labels using the TULIP LLM.

    Parameters
    ----------
    labels_path : str
        Directory containing per-study label JSON files produced by the pipeline.
        Each file is named ``{GSE_ID}.json`` and contains
        ``{GSM_ID: {label_axis: [canonical_value, ...]}}`` mappings.
    raw_data_path : str
        Root directory of the processed microarray data, containing one
        subdirectory per study (``{GSE_ID}/{GSM_ID}.json``).
    output_path : str
        Directory where evaluation results will be written.
    model : str
        TULIP model to use.
    timeout : int
        HTTP request timeout in seconds.
    """

    def __init__(
        self,
        labels_path: str,
        raw_data_path: str,
        output_path: str,
        model: str = TULIP_CHAT_MODEL,
        timeout: int = 60,
    ) -> None:

        self.labels_path = labels_path
        self.raw_data_path = raw_data_path
        self.output_path = output_path
        self.model = model
        self._client = OpenAI(
            base_url=TULIP_BASE_URL,
            api_key=TULIP_API_KEY,
            timeout=timeout,
        )
        os.makedirs(output_path, exist_ok=True)
        logger.info("TulipEvaluator ready (model=%s)", self.model)

    # ── Public entry point ─────────────────────────────────────────────────────

    def run(
        self,
        studies: list[str] | None = None,
        max_samples: int | None = None,
    ) -> dict:
        """
        Evaluate all studies (or a subset) and write results to output_path.

        Parameters
        ----------
        studies : list[str], optional
            GSE IDs to evaluate.  If None, all studies found in labels_path
            are evaluated.
        max_samples : int, optional
            Maximum number of samples to evaluate per study.  Useful for a
            quick spot-check without burning through the full dataset.
            Samples are taken from the start of the label file (not randomly,
            so results are reproducible).

        Returns
        -------
        dict
            Global summary statistics.
        """
        study_files = self._discover_studies(studies)
        if not study_files:
            logger.warning("No label files found in %s", self.labels_path)
            return {}

        # Derive canonical options from the full label set so they match the
        # ontology version that was actually used to produce the labels.
        canonical_options = self._derive_canonical_options(study_files)
        logger.info(
            "Canonical options derived: %s",
            {k: len(v) for k, v in canonical_options.items()},
        )

        global_counts: dict[str, dict[str, int]] = {}  # label → verdict → count

        for study_id, label_file in study_files.items():
            study_result_path = os.path.join(self.output_path, f"{study_id}_eval.json")

            # Resume: skip if already evaluated
            if os.path.exists(study_result_path):
                logger.info("Skipping %s (already evaluated)", study_id)
                existing = self._load_json(study_result_path)
                self._accumulate_global(existing.get("per_sample", {}), global_counts)
                continue

            logger.info("Evaluating %s ...", study_id)
            study_result = self._evaluate_study(
                study_id=study_id,
                label_file=label_file,
                canonical_options=canonical_options,
                max_samples=max_samples,
            )

            # Write per-study result immediately (crash-safe)
            self._save_json(study_result_path, study_result)
            self._accumulate_global(study_result.get("per_sample", {}), global_counts)
            self._log_study_summary(study_id, study_result.get("summary", {}))

        global_summary = self._build_global_summary(global_counts)
        self._save_json(os.path.join(self.output_path, "global_summary.json"), global_summary)
        self._print_global_summary(global_summary)
        return global_summary

    # ── Study-level evaluation ─────────────────────────────────────────────────

    def _evaluate_study(
        self,
        study_id: str,
        label_file: str,
        canonical_options: dict[str, list[str]],
        max_samples: int | None,
    ) -> dict:
        """Evaluate all samples in one study and return the full result dict."""
        labels = self._load_json(label_file)
        raw_study_dir = os.path.join(self.raw_data_path, study_id)

        # Load the study-level metadata once (used as context for all samples)
        study_meta = self._load_study_metadata(raw_study_dir)

        per_sample: dict[str, dict] = {}
        sample_ids = list(labels.keys())
        if max_samples:
            sample_ids = sample_ids[:max_samples]

        for sample_id in sample_ids:
            assigned = labels[sample_id]
            raw_sample_meta = self._load_sample_metadata(raw_study_dir, sample_id)

            verdict = self._evaluate_sample(
                sample_id=sample_id,
                assigned=assigned,
                sample_metadata=raw_sample_meta,
                study_metadata=study_meta,
                canonical_options=canonical_options,
            )
            per_sample[sample_id] = verdict

        summary = self._build_study_summary(per_sample)
        return {"study_id": study_id, "per_sample": per_sample, "summary": summary}

    # ── Sample-level evaluation ────────────────────────────────────────────────

    def _evaluate_sample(
        self,
        sample_id: str,
        assigned: dict[str, list[str]],
        sample_metadata: dict,
        study_metadata: dict,
        canonical_options: dict[str, list[str]],
    ) -> dict:
        """
        Ask TULIP to evaluate the assigned labels for one sample.
        Returns the raw verdicts dict plus any parse error info.
        """
        prompt = self._build_eval_prompt(
            assigned=assigned,
            sample_metadata=sample_metadata,
            study_metadata=study_metadata,
            canonical_options=canonical_options,
        )

        try:
            response = self._client.chat.completions.create(
                model=self.model,
                messages=[
                    {"role": "system", "content": _SYSTEM_PROMPT},
                    {"role": "user", "content": prompt},
                ],
                max_tokens=TULIP_MAX_TOKENS,
                temperature=0.0,
            )
            msg = response.choices[0].message

            # The TULIP model is a reasoning model: it writes a chain-of-thought
            # into reasoning_content and the final answer into content.
            # If content is None (budget exhausted mid-reasoning), fall back to
            # extracting a JSON block from reasoning_content so we still get a
            # partial result rather than an empty verdict.
            raw_text = msg.content or ""
            if not raw_text:
                reasoning = getattr(msg, "reasoning_content", None) or ""
                if reasoning:
                    logger.warning(
                        "content was empty (finish_reason=%s); attempting to extract JSON from reasoning_content",
                        response.choices[0].finish_reason,
                    )
                    raw_text = reasoning

            verdicts = self._parse_verdicts(raw_text, list(assigned.keys()), canonical_options=canonical_options)
            return {
                "assigned": assigned,
                "verdicts": verdicts,
                "raw_response": raw_text,
                "error": None,
            }

        except Exception as exc:
            logger.warning("TULIP call failed for %s/%s: %s", sample_id, sample_id, exc)
            return {
                "assigned": assigned,
                "verdicts": {},
                "raw_response": None,
                "error": str(exc),
            }

    # ── Prompt construction ────────────────────────────────────────────────────

    def _build_eval_prompt(
        self,
        assigned: dict[str, list[str]],
        sample_metadata: dict,
        study_metadata: dict,
        canonical_options: dict[str, list[str]],
    ) -> str:
        # --- Metadata block ---
        meta_lines = []
        for field in _METADATA_FIELDS:
            val = sample_metadata.get(field)
            if val:
                if isinstance(val, list):
                    for item in val:
                        meta_lines.append(f"  {field}: {item}")
                else:
                    meta_lines.append(f"  {field}: {val}")
        # Append abbreviated study context
        for field in _STUDY_FIELDS:
            val = study_metadata.get(field, "")
            if val:
                # Truncate long study summaries to keep prompt size reasonable
                truncated = str(val)[:400] + ("..." if len(str(val)) > 400 else "")
                meta_lines.append(f"  study_{field}: {truncated}")
        metadata_block = "\n".join(meta_lines) if meta_lines else "  (no metadata available)"

        # --- Options block ---
        options_lines = []
        for label in assigned:
            opts = canonical_options.get(label, [])
            options_lines.append(f"  {label}: {', '.join(opts)}")
        options_block = "\n".join(options_lines)

        # --- Assigned block ---
        assigned_lines = []
        for label, values in assigned.items():
            assigned_lines.append(f"  {label}: {', '.join(values)}")
        assigned_block = "\n".join(assigned_lines)

        return _USER_PROMPT_TEMPLATE.format(
            metadata_block=metadata_block,
            options_block=options_block,
            assigned_block=assigned_block,
        )

    # ── Response parsing ───────────────────────────────────────────────────────

    def _parse_verdicts(self, raw_text: str, expected_labels: list[str], canonical_options: dict[str, list[str]] | None = None) -> dict:
        """
        Parse the JSON verdict block from the model response.
        Falls back gracefully if the model returns malformed output.
        """
        # Strip markdown fences if the model added them despite instructions
        cleaned = re.sub(r"```(?:json)?\s*|\s*```", "", raw_text).strip()

        try:
            parsed = json.loads(cleaned)
        except json.JSONDecodeError:
            # Try to extract a JSON object with a regex as a last resort
            match = re.search(r"\{.*\}", cleaned, re.DOTALL)
            if match:
                try:
                    parsed = json.loads(match.group())
                except json.JSONDecodeError:
                    logger.debug("Could not parse TULIP verdict JSON: %s", cleaned[:200])
                    return self._unknown_verdicts(expected_labels, reason="JSON parse failed")
            else:
                return self._unknown_verdicts(expected_labels, reason="No JSON found in response")

        # Validate and normalise each label's verdict
        result = {}
        for label in expected_labels:
            raw_verdict = parsed.get(label, {})
            allowed = (canonical_options or {}).get(label, [])
            result[label] = self._normalise_verdict(raw_verdict, label, allowed)

        return result

    @staticmethod
    def _normalise_verdict(raw: dict, label: str, allowed_values: list[str] | None = None) -> dict:
        """
        Ensure a single verdict dict has the expected structure.

        If allowed_values is provided, suggestions are validated against it
        using a case-insensitive match — so "whole plant" is corrected to
        "Whole Plant" (the exact canonical form) rather than rejected.
        If no case-insensitive match exists the verdict is downgraded to
        "uncertain" to avoid storing hallucinated suggestions.
        """
        if not isinstance(raw, dict):
            return {
                "verdict": VERDICT_UNCERTAIN,
                "suggestion": None,
                "reason": "Malformed verdict from model",
            }
        verdict = raw.get("verdict", VERDICT_UNCERTAIN)
        if verdict not in VALID_VERDICTS:
            verdict = VERDICT_UNCERTAIN

        suggestion = raw.get("suggestion")
        if suggestion and not isinstance(suggestion, str):
            suggestion = str(suggestion)

        # Validate and canonicalise the suggestion when we have an allowed list
        if verdict == VERDICT_INCORRECT and suggestion and allowed_values:
            lower_map = {v.lower(): v for v in allowed_values}
            canonical = lower_map.get(suggestion.lower())
            if canonical:
                suggestion = canonical  # correct capitalisation
            else:
                # Model suggested a value outside the ontology — downgrade
                logger.debug("Suggestion '%s' for label '%s' not in canonical options; downgrading verdict to uncertain", suggestion, label)
                verdict = VERDICT_UNCERTAIN
                suggestion = None

        return {
            "verdict": verdict,
            "suggestion": suggestion if verdict == VERDICT_INCORRECT else None,
            "reason": str(raw.get("reason", ""))[:300],
        }

    @staticmethod
    def _unknown_verdicts(labels: list[str], reason: str = "") -> dict:
        return {label: {"verdict": VERDICT_UNCERTAIN, "suggestion": None, "reason": reason} for label in labels}

    # ── Summary helpers ────────────────────────────────────────────────────────

    @staticmethod
    def _build_study_summary(per_sample: dict) -> dict:
        """
        Aggregate verdicts across all samples in a study.
        Returns per-label counts and accuracy (correct / (correct + incorrect)).
        """
        counts: dict[str, dict[str, int]] = defaultdict(lambda: {VERDICT_CORRECT: 0, VERDICT_INCORRECT: 0, VERDICT_UNCERTAIN: 0})
        for sample_data in per_sample.values():
            for label, verdict_obj in sample_data.get("verdicts", {}).items():
                v = verdict_obj.get("verdict", VERDICT_UNCERTAIN)
                counts[label][v] += 1

        summary = {}
        for label, vc in counts.items():
            total_decided = vc[VERDICT_CORRECT] + vc[VERDICT_INCORRECT]
            accuracy = round(vc[VERDICT_CORRECT] / total_decided, 3) if total_decided > 0 else None
            summary[label] = {**vc, "accuracy": accuracy}
        return summary

    @staticmethod
    def _accumulate_global(per_sample: dict, global_counts: dict) -> None:
        """Merge a study's per-sample verdicts into the global counter."""
        for sample_data in per_sample.values():
            for label, verdict_obj in sample_data.get("verdicts", {}).items():
                if label not in global_counts:
                    global_counts[label] = {
                        VERDICT_CORRECT: 0,
                        VERDICT_INCORRECT: 0,
                        VERDICT_UNCERTAIN: 0,
                    }
                v = verdict_obj.get("verdict", VERDICT_UNCERTAIN)
                global_counts[label][v] += 1

    @staticmethod
    def _build_global_summary(global_counts: dict) -> dict:
        summary = {}
        for label, vc in global_counts.items():
            total_decided = vc[VERDICT_CORRECT] + vc[VERDICT_INCORRECT]
            accuracy = round(vc[VERDICT_CORRECT] / total_decided, 3) if total_decided > 0 else None
            summary[label] = {**vc, "accuracy": accuracy}
        return summary

    # ── Canonical options ──────────────────────────────────────────────────────

    @staticmethod
    def _derive_canonical_options(study_files: dict[str, str]) -> dict[str, list[str]]:
        """
        Scan all label files and collect unique values per label axis.
        This keeps the evaluator in sync with whatever ontology the pipeline used.
        """
        options: dict[str, set] = defaultdict(set)
        for path in study_files.values():
            try:
                data = json.loads(open(path).read())  # noqa: SIM115
                for sample_labels in data.values():
                    if not isinstance(sample_labels, dict):
                        continue
                    for label, values in sample_labels.items():
                        if isinstance(values, list):
                            for v in values:
                                if v not in ("unspecified", "unknown"):
                                    options[label].add(v)
                        elif isinstance(values, str) and values not in ("unspecified", "unknown"):
                            options[label].add(values)
            except Exception:
                pass
        return {label: sorted(vals) for label, vals in options.items()}

    # ── File / path helpers ────────────────────────────────────────────────────

    def _discover_studies(self, studies: list[str] | None) -> dict[str, str]:
        """Return {study_id: label_file_path} for studies to evaluate."""
        result = {}
        for fname in os.listdir(self.labels_path):
            if not fname.endswith(".json"):
                continue
            study_id = fname.replace(".json", "")
            if studies and study_id not in studies:
                continue
            result[study_id] = os.path.join(self.labels_path, fname)
        return result

    def _load_sample_metadata(self, study_dir: str, sample_id: str) -> dict:
        """Load the raw sample JSON and return only the sample_metadata dict."""
        path = os.path.join(study_dir, f"{study_dir.rsplit('/', maxsplit=1)[-1]}_{sample_id}.json")
        try:
            raw = self._load_json(path)
            return raw.get("sample_metadata", {})
        except Exception:
            return {}

    def _load_study_metadata(self, study_dir: str) -> dict:
        """
        Load study-level metadata. The pipeline stores it in each sample file
        under the 'study_metadata' key — we only need to read one sample file.
        """
        try:
            for fname in os.listdir(study_dir):
                if fname.endswith(".json"):
                    raw = self._load_json(os.path.join(study_dir, fname))
                    return raw.get("study_metadata", {})
        except Exception:
            pass
        return {}

    @staticmethod
    def _load_json(path: str) -> dict:
        with open(path) as f:
            return json.load(f)

    @staticmethod
    def _save_json(path: str, data: dict) -> None:
        with open(path, "w") as f:
            json.dump(data, f, indent=2)

    # ── Logging helpers ────────────────────────────────────────────────────────

    @staticmethod
    def _log_study_summary(study_id: str, summary: dict) -> None:
        lines = [f"\n  Study {study_id}:"]
        for label, stats in summary.items():
            acc = f"{stats['accuracy']:.1%}" if stats["accuracy"] is not None else "n/a"
            lines.append(f"    {label:<12} correct={stats[VERDICT_CORRECT]:3d}  incorrect={stats[VERDICT_INCORRECT]:3d}  uncertain={stats[VERDICT_UNCERTAIN]:3d}  accuracy={acc}")
        logger.info("\n".join(lines))

    @staticmethod
    def _print_global_summary(summary: dict) -> None:
        print("\n" + "=" * 60)
        print("GLOBAL EVALUATION SUMMARY")
        print("=" * 60)
        for label, stats in summary.items():
            acc = f"{stats['accuracy']:.1%}" if stats["accuracy"] is not None else "n/a"
            print(f"  {label:<12} correct={stats[VERDICT_CORRECT]:4d}  incorrect={stats[VERDICT_INCORRECT]:4d}  uncertain={stats[VERDICT_UNCERTAIN]:4d}  accuracy={acc}")
        print("=" * 60 + "\n")


# ── CLI entry point ────────────────────────────────────────────────────────────


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Evaluate pipeline-assigned GEO sample labels using TULIP LLM.")
    parser.add_argument("--labels-path", default="new_storage/labels/extractor_semantic_search/4.5", help="Directory containing per-study label JSON files (pipeline output).")
    parser.add_argument("--raw-data-path", default="new_storage/processed_microarray_data/", help="Root directory of processed_microarray_data/ (raw sample JSONs).")
    parser.add_argument("--output-path", default="outputs/eval_results", help="Directory where evaluation results will be written.")
    parser.add_argument("--studies", nargs="+", default=None, help="Optional list of GSE IDs to evaluate (default: all).")
    parser.add_argument("--max-samples", type=int, default=10, help="Max samples per study to evaluate (default: all).")
    parser.add_argument("--model", default=TULIP_CHAT_MODEL, help=f"TULIP model to use (default: {TULIP_CHAT_MODEL}).")
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    evaluator = TulipEvaluator(
        labels_path=args.labels_path,
        raw_data_path=args.raw_data_path,
        output_path=args.output_path,
        model=args.model,
    )
    evaluator.run(studies=args.studies, max_samples=args.max_samples)
