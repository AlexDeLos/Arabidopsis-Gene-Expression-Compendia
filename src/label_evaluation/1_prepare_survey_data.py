"""
1_prepare_survey_data.py
------------------------
Run this script ONCE before distributing the survey.
It reads your label files and metadata, selects a random sample of entries per study,
and packages everything into a single self-contained JSON file.

That JSON file is the only thing you need to share with reviewers —
they do not need access to your file system or any other data.

Usage
-----
    python 1_prepare_survey_data.py

Output
------
    src/label_evaluation/data/survey_data.json

Adjust the paths in the CONFIGURATION block below before running.
"""

import glob
import json
import os
import random

# ══════════════════════════════════════════════════════════════════
#  CONFIGURATION  — edit these paths before running
# ══════════════════════════════════════════════════════════════════

# Folder containing your per-study label JSON files (e.g. GSE30720.json)
LABELS_DIR = "./new_storage/labels/TULIP_1.2_RNA/5.0"

# Base folder for metadata. The script looks for:
#   {METADATA_BASE_DIR}/{study_id}/{study_id}_{sample_id}.json   (microarray)
#   {METADATA_BASE_DIR}/{study_id}/{study_id}_sample_metadata.csv  (RNA-seq, fallback)
METADATA_BASE_DIR = "./new_storage/rnaseq_data/metadata"

# Where to write the packaged survey file
OUTPUT_PAYLOAD = "src/label_evaluation/data/survey_data.json"

# How many samples to randomly select per study (set to None for all)
SAMPLES_PER_STUDY = 5

# Random seed for reproducibility (None = different each run)
RANDOM_SEED = 42

# ══════════════════════════════════════════════════════════════════


def _extract_metadata(meta_content: dict) -> tuple[str, str]:
    """
    Pull the two pieces of context the reviewer needs out of a metadata JSON:
      1. Characteristics (sample-level, most informative for labelling)
      2. Study-level context (title, summary, overall_design)

    Returns (characteristics_str, study_context_str).
    """
    # ── Sample characteristics ─────────────────────────────────────────────────
    sample_meta = meta_content.get("sample_metadata", meta_content)

    # characteristics_ch1 may be a list of "key: value" strings, a dict, or a plain string
    char_ch1 = sample_meta.get("characteristics_ch1", "")
    if isinstance(char_ch1, list):
        characteristics = "\n".join(str(c) for c in char_ch1 if c)
    elif isinstance(char_ch1, dict):
        characteristics = "\n".join(f"{k}: {v}" for k, v in char_ch1.items())
    else:
        characteristics = str(char_ch1).strip()

    # Supplement with title and source_name if present
    extra_fields = ["title", "source_name_ch1", "source_name", "description"]
    extra_lines = []
    for field in extra_fields:
        val = sample_meta.get(field, "")
        if val:
            val_str = val if isinstance(val, str) else "; ".join(str(v) for v in val)
            extra_lines.append(f"{field}: {val_str}")

    if extra_lines:
        characteristics = "\n".join(extra_lines) + "\n" + characteristics

    characteristics = characteristics.strip() or "N/A"

    # ── Study context ──────────────────────────────────────────────────────────
    study_meta = meta_content.get("study_metadata", {})
    priority_fields = ["title", "summary", "overall_design"]
    study_lines = []
    for field in priority_fields:
        val = study_meta.get(field, "")
        if val:
            val_str = val if isinstance(val, str) else " ".join(str(v) for v in val)
            # Truncate very long summaries so the display stays readable
            if len(val_str) > 1200:
                val_str = val_str[:1200] + "..."
            study_lines.append(f"{field.upper()}:\n{val_str}")

    study_context = "\n\n".join(study_lines).strip() or "N/A"

    return characteristics, study_context


def _format_label_for_display(labels: dict) -> list[dict]:
    """
    Convert the raw labels dict into a flat list of per-label entries
    ready for the reviewer to score.

    Each entry:
        {
            "label_category":  "treatment",
            "display_value":   "Heat (intensity: 2)",
            "raw_value":       [{"val": "Heat", "intensity": 2}]
        }

    Handles both:
      - Simple labels:     {"tissue": ["root"]}
      - Sub-attribute:     {"treatment": [{"val": "Heat", "intensity": 2}]}
    """
    entries = []
    for category, value in labels.items():
        raw = value

        # ── Format for human-readable display ─────────────────────────────────
        if isinstance(value, list):
            display_parts = []
            for item in value:
                if isinstance(item, dict):
                    # Sub-attribute format (e.g. treatment with intensity)
                    val_str = item.get("val", "?")
                    extras = {k: v for k, v in item.items() if k != "val"}
                    if extras:
                        extra_str = ", ".join(f"{k}: {v}" for k, v in extras.items())
                        display_parts.append(f"{val_str}  ({extra_str})")
                    else:
                        display_parts.append(str(val_str))
                else:
                    display_parts.append(str(item))
            display_value = " | ".join(display_parts)
        elif isinstance(value, dict):
            val_str = value.get("val", "?")
            extras = {k: v for k, v in value.items() if k != "val"}
            extra_str = ", ".join(f"{k}: {v}" for k, v in extras.items())
            display_value = f"{val_str}  ({extra_str})" if extras else val_str
        else:
            display_value = str(value)

        entries.append(
            {
                "label_category": category,
                "display_value": display_value,
                "raw_value": raw,
            }
        )

    return entries


def main():
    if RANDOM_SEED is not None:
        random.seed(RANDOM_SEED)

    os.makedirs(os.path.dirname(OUTPUT_PAYLOAD), exist_ok=True)

    survey_payload = []
    label_files = sorted(glob.glob(os.path.join(LABELS_DIR, "*.json")))

    # Filter out the aggregated file
    label_files = [f for f in label_files if os.path.basename(f).replace(".json", "") != "tulip_condensed_labels"]

    print(f"Found {len(label_files)} study label files.")
    print(f"Selecting up to {SAMPLES_PER_STUDY} random samples per study...\n")

    for label_path in label_files:
        study_id = os.path.basename(label_path).replace(".json", "")

        try:
            with open(label_path) as f:
                study_labels = json.load(f)
        except Exception as e:
            print(f"  [SKIP] {study_id}: could not read label file ({e})")
            continue

        if not isinstance(study_labels, dict):
            print(f"  [SKIP] {study_id}: unexpected label format")
            continue

        sample_ids = list(study_labels.keys())
        if SAMPLES_PER_STUDY and len(sample_ids) > SAMPLES_PER_STUDY:
            sample_ids = random.sample(sample_ids, SAMPLES_PER_STUDY)

        study_meta_dir = os.path.join(METADATA_BASE_DIR, study_id)
        samples_added = 0

        for sample_id in sample_ids:
            raw_labels = study_labels[sample_id]

            # ── Find metadata ──────────────────────────────────────────────────
            characteristics = "N/A"
            study_context = "N/A"
            sample_meta = {}
            meta_search = glob.glob(os.path.join(study_meta_dir, f"*{sample_id}*.json"))
            if meta_search:
                try:
                    with open(meta_search[0]) as mf:
                        meta_content = json.load(mf)
                    characteristics, study_context = _extract_metadata(meta_content)
                    sample_meta = meta_content.get("sample_metadata", meta_content)
                except Exception as e:
                    print(f"    [WARN] Could not read metadata for {sample_id}: {e}")
            else:
                print(f"    [WARN] No metadata file found for {sample_id} in {study_meta_dir}")

            # ── Format labels for per-label scoring ────────────────────────────
            label_entries = _format_label_for_display(raw_labels)

            survey_payload.append(
                {
                    "study_id": study_id,
                    "sample_id": sample_id,
                    "characteristics": characteristics,
                    "study_context": study_context,
                    "full_sample_metadata": sample_meta,  # <--- ADD THIS LINE
                    "label_entries": label_entries,
                }
            )
            samples_added += 1

        print(f"  {study_id}: {samples_added} samples added")

    with open(OUTPUT_PAYLOAD, "w", encoding="utf-8") as out_f:
        json.dump(survey_payload, out_f, indent=4, ensure_ascii=False)

    print(f"\n✅  Wrote {len(survey_payload)} samples → {OUTPUT_PAYLOAD}")
    print("    Share this file + 2_survey_app.py with your reviewers.")


if __name__ == "__main__":
    main()
