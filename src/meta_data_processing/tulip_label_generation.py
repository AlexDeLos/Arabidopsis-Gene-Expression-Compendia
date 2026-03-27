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

from src.constants import STORAGE_DIR, EXPERIMENT_NAME, LABELS_PATH
from src.constants_labeling import LABEL_CONFIG, UNIQUE_LABELS, LABELS, BUCKET_KEYWORDS, EXPLICIT_KEYWORDS

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

_SAMPLE_FIELDS = [
    "title", "source_name_ch1", "characteristics_ch1",
    "description", "molecule_ch1", "extract_protocol_ch1",
    "data_processing",
]
_STUDY_FIELDS = ["summary", "overall_design"]

# We use XML tags to clearly separate the data payload from the system instructions.
_USER_PROMPT_TEMPLATE = """\
Please analyze the following GEO metadata and extract the requested labels according to the system instructions.

<input>
<sample_metadata>
{sample_block}
</sample_metadata>

<study_context>
{study_block}
</study_context>
</input>
"""
import json

def generate_system_prompt():
    prompt = (
        "<context>You are an expert plant biology experiment annotator. Extract standardized metadata labels "
        "for Arabidopsis thaliana samples based on the provided GEO description. It is crucial you keep your decision concise and to the point, DO NOT ADD EXPLANATIONS\n</context>\n"
        "<instructions>\n"
        "1. Map the text to the exact canonical values provided in the <ontology>. Use the provided definitions to guide your categorization.\n"
        "2. Return ONLY a strictly valid JSON object. Do NOT wrap the JSON in markdown blocks (e.g., ```json) or include text.\n"
        "3. The response will be automatically parsed by a JSON parser  .\n"
        "</instructions>\n\n"
        "<rules>\n"
    )
    
    # 1. Enforce Unique Labels & Globalize Fallbacks
    for label in LABELS:
        if label in UNIQUE_LABELS:
            prompt += f"- '{label}': EXACTLY ONE canonical value.\n"
        else:
            prompt += f"- '{label}': ONE OR MULTIPLE canonical values.\n"
            
    prompt += "- FALLBACK: If a category is missing/unclear, output 'unspecified'. If explicitly stated as unknown, output 'unknown'.\n"
    prompt += "</rules>\n\n"
            
    # 2. Dynamically build ontology using Definitions
    prompt += "<ontology>\n"
    for label, config in LABEL_CONFIG.items():
        prompt += f"  <{label}>\n"
        
        for enum_item in config['enum']:
            canonical = enum_item.value
            
            # Skip printing unknown/unspecified since we made them global rules
            if canonical.lower() in ['unknown', 'unspecified']:
                continue
                
            # Fetch the definition instead of synonyms
            descriptions_dict = config.get('descriptions', {})
            description = descriptions_dict.get(enum_item, "")
            
            # Format cleanly: "- Canonical - Definition"
            if description:
                prompt += f"    - {canonical} - {description}\n"
            else:
                prompt += f"    - {canonical}\n"

        # Handle Sub-Attributes (like Treatment Intensity)
        if 'sub_attributes' in config:
            for sub_key, sub_config in config['sub_attributes'].items():
                prompt += f"\n    <sub_attribute name='{sub_key}'>\n"
                prompt += f"      {sub_config['instruction']}\n"
                for val_enum, desc in sub_config['descriptions'].items():
                     val_str = getattr(val_enum, 'value', val_enum)
                     prompt += f"        - {val_str} = {desc}\n"
                prompt += f"    </sub_attribute>\n"
                
        prompt += f"  </{label}>\n"
    prompt += "</ontology>\n\n"

    # 3. Output Format
    prompt += "<output_format>\n"
    
    schema_template = {}
    for label, config in LABEL_CONFIG.items():
        if 'sub_attributes' in config:
            obj_schema = {"val": f"CANONICAL_{label.upper()}"}
            for sub_key, sub_config in config['sub_attributes'].items():
                obj_schema[sub_key] = sub_config.get('type_hint', "integer") 
            schema_template[label] = [obj_schema]
        else:
            schema_template[label] = [f"CANONICAL_{label.upper()}"]

    schema_json_str = json.dumps(schema_template, indent=2)
    prompt += f"{schema_json_str}\n"
    prompt += "</output_format>\n\n"
    prompt += "FINAL REMINDER: Output ONLY valid JSON."
    
    return prompt


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
        self.saving_path = LABELS_PATH
        # saving_path or os.path.join(
        #     STORAGE_DIR, "labels", TULIP_MODEL_NAME, EXPERIMENT_NAME
        # )
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
        self.system_prompt = generate_system_prompt()
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
        study_metadata:  Dict):
        """
        Ask TULIP to assign all label axes for one sample.
        Returns a dict matching the standard pipeline output format:
            {label_axis: [canonical_value]}
        Falls back to {"label": ["unspecified"]} for any axis that fails.
        """
        prompt = self._build_prompt(sample_metadata, study_metadata)

        max_retries = 3
        for attempt in range(max_retries):
            try:
                response = self._client.chat.completions.create(
                    model       = self.model,
                    messages    = [
                        {"role": "system", "content": self.system_prompt},
                        {"role": "user",   "content": prompt},
                    ],
                    max_tokens  = TULIP_MAX_TOKENS,
                    temperature = 0.1*attempt,
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
                logger.warning(
                    "TULIP call failed on attempt %d/%d: %s", 
                    attempt + 1, max_retries, exc
                )
                
                if attempt == max_retries - 1:
                    logger.error("All %d attempts to reach TULIP failed. Falling back to default labels.", max_retries)
                    return self._parse_labels("")  # Generates the default unclassified schema
                
                # Optional: Add a standard `import time` at top of file and do:
                # time.sleep(1 * (attempt + 1)) # to give transient network blips time to resolve
            # return {label: ["unspecified"] for label in LABELS}

    # ── Prompt construction ────────────────────────────────────────────────────

    def _build_prompt(self, sample_metadata: Dict, study_metadata: Dict) -> str:
        # --- Format Sample Metadata ---
        sample_lines = []
        for field in _SAMPLE_FIELDS:
            val = sample_metadata.get(field)
            if not val:
                continue
            
            # Format lists cleanly rather than duplicating the field name
            if isinstance(val, list):
                clean_val = "; ".join([str(item).strip() for item in val if item])
                if clean_val:
                    sample_lines.append(f"  [{field.upper()}]: {clean_val}")
            else:
                sample_lines.append(f"  [{field.upper()}]: {str(val).strip()}")

        sample_block = "\n".join(sample_lines) if sample_lines else "  (No sample metadata available)"

        # --- Format Study Metadata (Context) ---
        study_lines = []
        for field in _STUDY_FIELDS:
            # Some GEO fields come back as lists, ensure we handle them safely
            val = study_metadata.get(field, "")
            if isinstance(val, list):
                val = " ".join([str(v) for v in val])
            else:
                val = str(val).strip()

            if val:
                # 400 characters is too short for GEO abstracts where the treatment protocol 
                # is often described at the very end. 1500 is much safer for modern LLMs.
                truncated = val[:1500] + ("..." if len(val) > 1500 else "")
                study_lines.append(f"  [{field.upper()}]: {truncated}")

        study_block = "\n".join(study_lines) if study_lines else "  (No study context available)"

        # Options block is intentionally removed here because the System Prompt 
        # now handles the entire ontology, rules, and JSON schema.
        return _USER_PROMPT_TEMPLATE.format(
            sample_block=sample_block,
            study_block=study_block,
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
                    raise ValueError('invalid response')
                    logger.debug("Could not parse label JSON: %s", cleaned[:200])
            else:

                    raise ValueError('invalid response')

        result: Dict[str, List] = {}  # Changed List[str] to List
        for label in LABELS:
            raw_val = parsed.get(label, "unspecified")
            
            # Coerce into a list if it's a single item or single string
            if not isinstance(raw_val, list):
                raw_val = [raw_val] if raw_val != "unspecified" else []

            validated = []
            sub_attr_cfg = LABEL_CONFIG[label].get('sub_attributes', {})

            for r in raw_val:
                # Handle cases where the LLM might return a flat string instead of a dict
                item_dict = r if isinstance(r, dict) else {'val': str(r)}
                
                # 1. Validate the primary value ('val')
                main_val = self._validate_label(item_dict.get('val', 'unspecified'), label)
                
                if sub_attr_cfg:
                    valid_item = {'val': main_val}
                    
                    # 2. Dynamically validate all sub-attributes
                    for sub_key, sub_cfg in sub_attr_cfg.items():
                        raw_sub_val = item_dict.get(sub_key)
                        enum_cls = sub_cfg['enum']
                        
                        try:
                            # Infer the expected type from the first enum value (e.g. int vs str)
                            expected_type = type(list(enum_cls)[0].value)
                            typed_val = expected_type(raw_sub_val)
                            
                            if any(typed_val == e.value for e in enum_cls):
                                valid_item[sub_key] = typed_val
                            else:
                                valid_item[sub_key] = list(enum_cls)[0].value
                        except (TypeError, ValueError, IndexError):
                            # Fallback if the parsing fails or value is missing
                            valid_item[sub_key] = list(enum_cls)[0].value if list(enum_cls) else "unspecified"
                    
                    validated.append(valid_item)
                else:
                    # Simple label validation without sub-attributes
                    validated.append(main_val)

            # Ensure we always return at least one entry per label
            if not validated:
                validated = ["unspecified"] if not sub_attr_cfg else [{"val": "unspecified", **{k: list(v['enum'])[0].value for k,v in sub_attr_cfg.items()}}]

            result[label] = validated

        return result
    

    def _validate_sub_attribute(self, raw_value, label: str, sub_key: str):
        """
        Dynamically validate a sub-attribute value against its enumeration in LABEL_CONFIG.
        Discovers type constraints and falls back to a standardized default if invalid.
        """
        try:
            sub_config = LABEL_CONFIG[label]['sub_attributes'][sub_key]
            enum_cls = sub_config['enum']
        except KeyError:
            return "unspecified"

        allowed_values = [e.value for e in enum_cls]
        if not allowed_values:
            return "unspecified"

        # Determine target type (e.g., int vs str) based on the first item in the Enum
        expected_type = type(allowed_values[0])

        try:
            cleaned = raw_value
            if isinstance(cleaned, str):
                cleaned = re.sub(r'^["\'\s]+|["\'\s]+$', '', cleaned).strip()

            typed_val = expected_type(cleaned)
            if typed_val in allowed_values:
                return typed_val
        except (ValueError, TypeError):
            pass

        # Smart fallback loop: Try to find a standardized default, or use the first enum
        for item in enum_cls:
            if item.name.lower() in ["unspecified", "control", "unknown"]:
                return item.value

        return allowed_values[0]
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
