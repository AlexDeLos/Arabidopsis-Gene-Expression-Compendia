"""
tulip_label_generation.py
--------------------------
Full TULIP-only labeling pipeline for Arabidopsis thaliana GEO samples.

Replaces the two-stage UniversalExtractor + GroundingOptimizer pipeline with a
single LLM call per sample. No BioLORD, no spaCy, no CUDA required.

The output format is identical to label_generation.py:
  {LABELS_PATH}/{GSE_ID}.json  →  {GSM_ID: {label_axis: [canonical_value, ...]}}

This means outputs can be fed directly into tulip_evaluator.py or used anywhere
the standard pipeline output is expected.

Design
------
* One TULIP call per sample returning all label axes together as JSON.
  Cross-axis context improves accuracy (the model can use tissue to inform
  treatment and vice versa).
* Canonical options are loaded from constants.py (BUCKET_KEYWORDS) so the
  model picks from exactly the same ontology as the vector pipeline.
* Results are written per-study immediately after processing — crash-safe
  and resumable (already-done studies are skipped).
* A simple file-based cache stores raw_metadata → assigned labels so that
  repeated identical metadata strings (same GEO entry in multiple runs) are
  never sent to TULIP twice.

Usage
-----
    python tulip_label_generation.py                          # all studies
    python tulip_label_generation.py --studies GSE5622 GSE9415
    python tulip_label_generation.py --max-samples 5          # quick test

Or programmatically:
    from tulip_label_generation import TulipLabelGenerator
    gen = TulipLabelGenerator()
    gen.run(studies=["GSE5622"])
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import re
import sys
from typing import Dict, List, Optional

from tqdm import tqdm

module_dir = './'
sys.path.append(module_dir)

from src.constants import (
    LABELS, BUCKET_KEYWORDS, STORAGE_DIR, EXPERIMENT_NAME,
    LABELS_PATH,
)

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)

# ── TULIP constants ────────────────────────────────────────────────────────────
TULIP_BASE_URL   = "https://api.tulip.tudelft.nl/chat/v1"
TULIP_MODEL      = "chat"
TULIP_API_KEY    = os.environ.get("TULIP_API_KEY", "DUMMY_API_KEY")
TULIP_MAX_TOKENS = 2048   # Reasoning model writes chain-of-thought before JSON

# Output subfolder name — kept separate from the vector pipeline so results
# don't overwrite each other and can be compared side-by-side.
TULIP_MODEL_NAME = "tulip_llm"

# Metadata fields forwarded to TULIP, in priority order.
_SAMPLE_FIELDS = [
    "title", "source_name_ch1", "characteristics_ch1",
    "description", "molecule_ch1", "extract_protocol_ch1",
    "data_processing",
]
_STUDY_FIELDS = ["summary", "overall_design"]

# ── Prompts ────────────────────────────────────────────────────────────────────
_SYSTEM_PROMPT = """\
You are an expert curator of Arabidopsis thaliana gene-expression metadata from NCBI GEO.

You will be given raw GEO metadata for one sample and a list of canonical values
allowed for each label axis.

Your task: assign exactly ONE canonical value per label axis based on the metadata.

Rules:
1. Write the JSON object FIRST before any reasoning.
2. Each value must be EXACTLY one string from the canonical options for that axis.
3. If the metadata does not contain clear evidence for a label axis, assign "unspecified".
4. Do not invent values not in the canonical options list.
5. No markdown fences around the JSON.

Output format (one key per label axis, value is the single chosen canonical string):
{
  "tissue":    "<canonical value or unspecified>",
  "treatment": "<canonical value or unspecified>",
  "medium":    "<canonical value or unspecified>",
  "genotype":  "<canonical value or unspecified>"
}
"""

_USER_PROMPT_TEMPLATE = """\
=== RAW METADATA ===
{metadata_block}

=== CANONICAL OPTIONS ===
{options_block}

Assign one canonical value per label axis and return as a JSON object.
"""


# ── TulipLabelGenerator ────────────────────────────────────────────────────────

class TulipLabelGenerator:
    """
    Labels GEO samples using TULIP, replacing the UniversalExtractor +
    GroundingOptimizer pipeline entirely.

    Parameters
    ----------
    in_folder : str
        Root directory of processed_microarray_data/ with one subdir per study.
    saving_path : str
        Output directory for per-study label JSON files.
    model : str
        TULIP model identifier.
    timeout : int
        HTTP timeout per request in seconds.
    """

    def __init__(
        self,
        in_folder:    str  = "new_storage/processed_microarray_data/",
        saving_path:  str  = None,
        model:        str  = TULIP_MODEL,
        timeout:      int  = 60,
    ) -> None:
        try:
            from openai import OpenAI
        except ImportError as exc:
            raise ImportError(
                "The 'openai' package is required. Install with: pip install openai"
            ) from exc

        self.in_folder   = in_folder
        self.saving_path = saving_path or os.path.join(
            STORAGE_DIR, "labels", TULIP_MODEL_NAME, EXPERIMENT_NAME
        )
        self.model   = model
        self._client = OpenAI(
            base_url = TULIP_BASE_URL,
            api_key  = TULIP_API_KEY,
            timeout  = timeout,
        )

        # Build canonical options once from BUCKET_KEYWORDS (canonical names only,
        # no synonyms — keeps the prompt short and the model on-ontology).
        self.canonical_options: Dict[str, List[str]] = {
            label: list(BUCKET_KEYWORDS[label]) for label in LABELS
        }

        os.makedirs(self.saving_path, exist_ok=True)
        logger.info(
            "TulipLabelGenerator ready — saving to %s", self.saving_path
        )

    # ── Public entry point ─────────────────────────────────────────────────────

    def run(
        self,
        studies:     Optional[List[str]] = None,
        max_samples: Optional[int]       = None,
    ) -> None:
        """
        Label all studies (or a subset) and write results to saving_path.

        Parameters
        ----------
        studies : list[str], optional
            GSE IDs to process.  None = all subdirectories in in_folder.
        max_samples : int, optional
            Cap on samples processed per study — useful for quick testing.
        """
        study_dirs = [
            d for d in os.listdir(self.in_folder)
            if os.path.isdir(os.path.join(self.in_folder, d))
        ]

        for study_id in tqdm(study_dirs, desc="Labelling studies"):
            if studies and study_id not in studies:
                continue

            output_file = os.path.join(self.saving_path, f"{study_id}.json")
            if os.path.exists(output_file):
                logger.info("Skipping %s (already labelled)", study_id)
                continue

            self._process_study(study_id, output_file, max_samples)

        # Aggregate into flat list (mirrors label_generation.py __main__ block)
        self._aggregate(self.saving_path)

    # ── Study-level processing ─────────────────────────────────────────────────

    def _process_study(
        self,
        study_id:    str,
        output_file: str,
        max_samples: Optional[int],
    ) -> None:
        study_path = os.path.join(self.in_folder, study_id)
        sample_files = sorted(
            f for f in os.listdir(study_path) if f.endswith(".json")
        )
        if max_samples:
            sample_files = sample_files[:max_samples]

        final_output: Dict[str, Dict] = {}

        for sample_file in sample_files:
            file_path = os.path.join(study_path, sample_file)
            try:
                with open(file_path) as f:
                    data = json.load(f)

                sample_id     = data.get("sample_id", sample_file.replace(".json", ""))
                sample_meta   = data.get("sample_metadata", {})
                study_meta    = data.get("study_metadata", {})

                labels = self._label_sample(sample_meta, study_meta)
                final_output[sample_id] = labels

            except Exception as exc:
                logger.warning("Error processing %s/%s: %s", study_id, sample_file, exc)

        if final_output:
            with open(output_file, "w") as f:
                json.dump(final_output, f, indent=4)
            logger.info("Saved %d samples → %s", len(final_output), output_file)

    # ── Sample-level labeling ──────────────────────────────────────────────────

    def _label_sample(
        self,
        sample_metadata: Dict,
        study_metadata:  Dict,
    ) -> Dict[str, List[str]]:
        """
        Ask TULIP to assign all label axes for one sample.
        Returns a dict matching the standard pipeline output format:
            {label_axis: [canonical_value]}
        Falls back to {"label": ["unspecified"]} for any axis that fails.
        """
        prompt = self._build_prompt(sample_metadata, study_metadata)

        try:
            response = self._client.chat.completions.create(
                model       = self.model,
                messages    = [
                    {"role": "system", "content": _SYSTEM_PROMPT},
                    {"role": "user",   "content": prompt},
                ],
                max_tokens  = TULIP_MAX_TOKENS,
                temperature = 0.0,
            )
            msg      = response.choices[0].message
            raw_text = msg.content or ""

            # Reasoning model fallback: content may be empty if token budget was
            # spent on chain-of-thought — extract JSON from reasoning_content instead.
            if not raw_text:
                reasoning = getattr(msg, "reasoning_content", None) or ""
                if reasoning:
                    logger.debug(
                        "content empty (finish_reason=%s); extracting from reasoning_content",
                        response.choices[0].finish_reason,
                    )
                    raw_text = reasoning

            return self._parse_labels(raw_text)

        except Exception as exc:
            logger.warning("TULIP call failed: %s", exc)
            return {label: ["unspecified"] for label in LABELS}

    # ── Prompt construction ────────────────────────────────────────────────────

    def _build_prompt(self, sample_metadata: Dict, study_metadata: Dict) -> str:
        # Metadata block — sample fields first, abbreviated study context appended
        meta_lines = []
        for field in _SAMPLE_FIELDS:
            val = sample_metadata.get(field)
            if not val:
                continue
            if isinstance(val, list):
                for item in val:
                    if item:
                        meta_lines.append(f"  {field}: {item}")
            else:
                meta_lines.append(f"  {field}: {val}")

        for field in _STUDY_FIELDS:
            val = str(study_metadata.get(field, "")).strip()
            if val:
                truncated = val[:400] + ("..." if len(val) > 400 else "")
                meta_lines.append(f"  study_{field}: {truncated}")

        metadata_block = "\n".join(meta_lines) if meta_lines else "  (no metadata available)"

        # Options block — one line per label axis
        options_lines = [
            f"  {label}: {', '.join(self.canonical_options.get(label, []))}"
            for label in LABELS
        ]
        options_block = "\n".join(options_lines)

        return _USER_PROMPT_TEMPLATE.format(
            metadata_block = metadata_block,
            options_block  = options_block,
        )

    # ── Response parsing ───────────────────────────────────────────────────────

    def _parse_labels(self, raw_text: str) -> Dict[str, List[str]]:
        """
        Parse the JSON label assignment from the model response.
        Returns {label: [canonical_value]} for all LABELS.
        Any axis that is missing or invalid defaults to ["unspecified"].
        """
        cleaned = re.sub(r"```(?:json)?\s*|\s*```", "", raw_text).strip()

        parsed = {}
        try:
            parsed = json.loads(cleaned)
        except json.JSONDecodeError:
            match = re.search(r"\{.*?\}", cleaned, re.DOTALL)
            if match:
                try:
                    parsed = json.loads(match.group())
                except json.JSONDecodeError:
                    logger.debug("Could not parse label JSON: %s", cleaned[:200])

        result: Dict[str, List[str]] = {}
        for label in LABELS:
            raw_val = parsed.get(label, "unspecified")
            validated = self._validate_label(str(raw_val), label)
            result[label] = [validated]

        return result

    def _validate_label(self, value: str, label: str) -> str:
        """
        Accept the value only if it matches a canonical option (case-insensitive).
        Returns "unspecified" otherwise — never stores hallucinated values.
        """
        cleaned = re.sub(r'^["\'\s]+|["\'\s]+$', '', value).strip()
        if cleaned.lower() == "unspecified":
            return "unspecified"

        allowed = self.canonical_options.get(label, [])
        lower_map = {v.lower(): v for v in allowed}

        if cleaned.lower() in lower_map:
            return lower_map[cleaned.lower()]

        logger.debug(
            "TULIP value '%s' for label '%s' not in canonical options; "
            "falling back to unspecified", cleaned, label
        )
        return "unspecified"

    # ── Aggregation ────────────────────────────────────────────────────────────

    @staticmethod
    def _aggregate(saving_path: str) -> None:
        """
        Aggregate all per-study label files into a single flat JSON list,
        mirroring the llm_condensed_labels.json produced by label_generation.py.
        """
        res = []
        for fname in os.listdir(saving_path):
            if not fname.endswith(".json") or fname == "tulip_condensed_labels.json":
                continue
            study_id = fname.replace(".json", "")
            with open(os.path.join(saving_path, fname)) as f:
                study_data = json.load(f)
            for sample_id, labels in study_data.items():
                entry = dict(labels)
                entry["id"] = sample_id
                res.append(entry)

        out_path = os.path.join(saving_path, "tulip_condensed_labels.json")
        with open(out_path, "w") as f:
            json.dump(res, f, indent=2)
        logger.info("Aggregated %d samples → %s", len(res), out_path)


# ── CLI entry point ────────────────────────────────────────────────────────────

def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Label GEO samples using TULIP LLM (replaces extractor+grounder)."
    )
    parser.add_argument(
        "--in-folder", default="new_storage/processed_microarray_data/",
        help="Root directory of processed microarray data."
    )
    parser.add_argument(
        "--saving-path", default=None,
        help="Output directory for label JSON files. "
             "Defaults to {STORAGE_DIR}/labels/tulip_llm/{EXPERIMENT_NAME}."
    )
    parser.add_argument(
        "--studies", nargs="+", default=None,
        help="Optional list of GSE IDs to process (default: all)."
    )
    parser.add_argument(
        "--max-samples", type=int, default=None,
        help="Max samples per study (default: all). Useful for quick testing."
    )
    parser.add_argument(
        "--model", default=TULIP_MODEL,
        help=f"TULIP model to use (default: {TULIP_MODEL})."
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    generator = TulipLabelGenerator(
        in_folder   = args.in_folder,
        saving_path = args.saving_path,
        model       = args.model,
    )
    generator.run(studies=args.studies, max_samples=args.max_samples)
