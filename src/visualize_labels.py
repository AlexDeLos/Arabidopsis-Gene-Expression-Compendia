import ast
import json
import os
import sys

import matplotlib as mpl

# Force matplotlib to use a non-interactive backend so it doesn't crash looking for a GUI
mpl.use("Agg")
from collections import Counter

import matplotlib.pyplot as plt

# Ensure the script can find your src modules if run from the root directory
sys.path.append(os.path.abspath("./"))
from src.constants import FIGURES_DIR, LABELS_PATH, RNA_USED
from src.constants_labeling import LABELS

# Keys (in priority order) that hold the primary label value inside a dict item.
# "val" is checked first to match the {"val": "Control", "intensity": 0} format.
_VALUE_KEYS = ("val", "value")

# Keys inside dicts that carry metadata rather than a label value
_IGNORE_KEYS = {"val", "value", "confidence", "reasoning", "evidence", "explanation", "source"}

# Strings that should be treated as missing data
_IGNORE_WORDS = {"none", "unspecified", "unknown", "na", "n/a", "null", "", "false", "undefined"}


def load_all_labels(labels_dir):
    """
    Loads all sample JSON files from the labels directory.
    Supports both dict-keyed files {sample_id: {...}} and flat list files [{...}].
    Files starting with 'map_' are skipped (they are lookup tables, not sample data).
    """
    all_samples = []
    if not os.path.exists(labels_dir):
        print(f"Error: Directory {labels_dir} does not exist.")
        return all_samples

    files = sorted(f for f in os.listdir(labels_dir) if f.endswith(".json") and not f.startswith("map_"))
    for file in files:
        file_path = os.path.join(labels_dir, file)
        with open(file_path) as f:
            try:
                study_data = json.load(f)
                if isinstance(study_data, dict):
                    samples = list(study_data.values())
                elif isinstance(study_data, list):
                    samples = study_data
                else:
                    print(f"  [WARN] Unexpected format in {file}, skipping.")
                    continue

                valid = [s for s in samples if isinstance(s, dict)]
                all_samples.extend(valid)
                print(f"  Loaded {len(valid):>5} samples  ← {file}")

            except Exception as e:
                print(f"  [WARN] Could not read {file_path}: {e}")

    return all_samples


def clean_and_format_str(s):
    """Safely cleans strings and applies Title Case to merge casing fractures."""
    if not s:
        return ""
    s = str(s).strip("[]'\" {}")
    if not s or s.lower() in _IGNORE_WORDS:
        return ""
    return s.title()


def aggregate_label_counts(samples, label_category):
    """Counts the occurrences of each grounded term in a specific category."""
    counter = Counter()

    for sample in samples:
        if not isinstance(sample, dict):
            continue

        raw_items = sample.get(label_category, [])

        # Normalise to a list
        if isinstance(raw_items, (dict, str)):
            items = [raw_items]
        elif raw_items is None:
            items = []
        else:
            items = list(raw_items)

        if not items:
            counter["Unspecified"] += 1
            continue

        for item in items:
            # Catch stringified JSON dicts from the LLM
            if isinstance(item, str) and item.strip().startswith("{"):
                try:  # noqa: SIM105
                    item = ast.literal_eval(item)  # noqa: PLW2901
                except (ValueError, SyntaxError):
                    pass

            if isinstance(item, dict):
                # Resolve primary value — check "val" before "value"
                base_raw = None
                for key in _VALUE_KEYS:
                    if key in item:
                        base_raw = item[key]
                        break

                if isinstance(base_raw, list):
                    base_raw = base_raw[0] if base_raw else None

                base_value = clean_and_format_str(base_raw) or "Unspecified"

                # Collect sub-attributes (e.g. intensity), skipping metadata keys
                sub_vals = []
                for key, val in item.items():
                    if key.lower() in _IGNORE_KEYS:
                        continue
                    if isinstance(val, list):
                        val = ", ".join(str(v) for v in val)  # noqa: PLW2901
                    cleaned = clean_and_format_str(str(val))
                    if cleaned:
                        sub_vals.append(f"{key.capitalize()}: {cleaned}")

                formatted_term = f"{base_value} ({', '.join(sub_vals)})" if sub_vals else base_value
                counter[formatted_term] += 1

            else:
                clean_item = clean_and_format_str(item)
                counter[clean_item if clean_item else "Unspecified"] += 1

    return counter


if __name__ == "__main__":
    LABELS_DIR = LABELS_PATH

    print(f"Loading data from: {LABELS_DIR}")
    samples = load_all_labels(LABELS_DIR)
    print(f"\nTotal samples pooled: {len(samples)}")

    if not samples:
        print("No samples found. Exiting.")
        sys.exit(0)

    prefix = "RNA" if RNA_USED else "Microarray"
    output_dir = os.path.join(FIGURES_DIR, f"label_distributions_{prefix}")
    os.makedirs(output_dir, exist_ok=True)

    print(f"Saving plots to: {output_dir}\n")

    for label in LABELS:
        counts = aggregate_label_counts(samples, label)
        top_counts = counts.most_common(20)

        if not top_counts:
            print(f"  [SKIP] No data for label '{label}'.")
            continue

        terms = [t[0] for t in top_counts][::-1]
        frequencies = [t[1] for t in top_counts][::-1]

        # Scale figure height to the number of bars so labels are never cramped
        fig_height = max(5, len(terms) * 0.45 + 1.5)
        fig, ax = plt.subplots(figsize=(11, fig_height))

        colors = [
            "#cccccc" if any(w in t.lower() for w in ("unspecified", "unknown")) else "#4C72B0"
            for t in terms
        ]

        ax.barh(terms, frequencies, color=colors, edgecolor="none", alpha=0.85)

        x_pad = max(frequencies) * 0.012
        for i, freq in enumerate(frequencies):
            ax.text(freq + x_pad, i, str(freq), va="center", ha="left", fontsize=9, color="#444444")

        ax.set_title(f"{prefix} — Distribution of {label.replace('_', ' ').capitalize()} labels",
                     fontsize=13, pad=12)
        ax.set_xlabel("Number of samples", fontsize=11)
        ax.set_ylabel(f"{label.replace('_', ' ').capitalize()} term", fontsize=11)

        # Extend x-axis so inline count labels are never clipped
        ax.set_xlim(right=max(frequencies) * 1.12)

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(axis="both", labelsize=9)

        plt.tight_layout()

        save_path = os.path.join(output_dir, f"{prefix}_{label}_distribution.png")
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"  Saved → {save_path}")

    print("\nDone.")