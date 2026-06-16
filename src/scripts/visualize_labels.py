import ast
import json
import os
import sys
import pandas as pd
from collections import Counter, defaultdict

import matplotlib as mpl

# Force matplotlib to use a non-interactive backend so it doesn't crash looking for a GUI
mpl.use("Agg")

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

# Ensure the script can find your src modules if run from the root directory
sys.path.append(os.path.abspath("./"))
from src.constants import FIGURES_DIR, LABELS_PATH, RNA_USED, EXPR_PATH
from src.constants_labeling import LABELS, IntensityEnum
from src.data_analisys.utils.cluster_exploration_utils_final import get_gsm_id

# Keys (in priority order) that hold the primary label value inside a dict item.
# "val" is checked first to match the {"val": "Control", "intensity": 0} format.
_VALUE_KEYS = ("val", "value")

# Keys inside dicts that carry metadata rather than a label value
_IGNORE_KEYS = {"val", "value", "confidence", "reasoning", "evidence", "explanation", "source"}

# Strings that should be treated as missing data
_IGNORE_WORDS = {"none", "unspecified", "unknown", "na", "n/a", "null", "", "false", "undefined"}

# Intensity level display config: IntensityEnum int value -> (legend label, hex colour)
_INTENSITY_STYLE = {
    IntensityEnum.CONTROL.value:  ("Control (0)",  "#4CAF50"),  # green
    IntensityEnum.MILD.value:     ("Mild (1)",      "#FF9800"),  # orange
    IntensityEnum.MODERATE.value: ("Moderate (2)",  "#F44336"),  # red
    IntensityEnum.SEVERE.value:   ("Severe (3)",    "#7B0000"),  # dark red
}
_INTENSITY_ORDER = [
    IntensityEnum.CONTROL.value,
    IntensityEnum.MILD.value,
    IntensityEnum.MODERATE.value,
    IntensityEnum.SEVERE.value,
]


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def load_all_labels(labels_dir):
    """
    Loads all sample JSON files and injects the 'Sample_ID' from dictionary keys.
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
                processed_samples = []

                if isinstance(study_data, dict):
                    # Iterate through the dictionary to keep the GSM ID (the key)
                    for sample_id, sample_meta in study_data.items():
                        if isinstance(sample_meta, dict):
                            # Inject the ID into the dictionary so it's accessible later
                            sample_meta['Sample_ID'] = sample_id
                            processed_samples.append(sample_meta)
                
                elif isinstance(study_data, list):
                    # If it's already a list, ensure each item is a dict
                    processed_samples = [s for s in study_data if isinstance(s, dict)]

                all_samples.extend(processed_samples)
                print(f"  Loaded {len(processed_samples):>5} samples  <- {file}")

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


def _resolve_item_base(item):
    """Extract the primary string value from a dict or plain-string item.
    Returns (item_as_dict_or_None, base_value_string).
    """
    if isinstance(item, str) and item.strip().startswith("{"):
        try:
            item = ast.literal_eval(item)
        except (ValueError, SyntaxError):
            pass

    if isinstance(item, dict):
        base_raw = None
        for key in _VALUE_KEYS:
            if key in item:
                base_raw = item[key]
                break
        if isinstance(base_raw, list):
            base_raw = base_raw[0] if base_raw else None
        return item, clean_and_format_str(base_raw) or "Unspecified"

    return None, clean_and_format_str(item) or "Unspecified"


# ---------------------------------------------------------------------------
# Standard label aggregation (all labels except treatment)
# ---------------------------------------------------------------------------

def aggregate_label_counts(samples, label_category):
    """Counts the occurrences of each grounded term in a specific category."""
    counter = Counter()

    for sample in samples:
        if not isinstance(sample, dict):
            continue

        raw_items = sample.get(label_category, [])

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
            item_dict, base_value = _resolve_item_base(item)

            if item_dict is not None:
                sub_vals = []
                for key, val in item_dict.items():
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
                counter[base_value] += 1

    return counter


# ---------------------------------------------------------------------------
# Treatment-specific aggregation (groups by name, splits by intensity)
# ---------------------------------------------------------------------------

def aggregate_treatment_intensity_counts(samples):
    """
    Groups treatment samples by treatment name and counts how many samples
    fall into each intensity level (0-3) per treatment.

    Returns: {treatment_name: {intensity_int: count}}
    """
    counts = defaultdict(lambda: defaultdict(int))

    for sample in samples:
        if not isinstance(sample, dict):
            continue

        raw_items = sample.get("treatment", [])
        if isinstance(raw_items, (dict, str)):
            items = [raw_items]
        elif raw_items is None:
            items = []
        else:
            items = list(raw_items)

        if not items:
            counts["Unspecified"][IntensityEnum.CONTROL.value] += 1
            continue

        for item in items:
            item_dict, name = _resolve_item_base(item)

            intensity = IntensityEnum.CONTROL.value
            if item_dict is not None:
                try:
                    raw_intensity = item_dict.get("intensity", IntensityEnum.CONTROL.value)
                    intensity = int(raw_intensity)
                    if intensity not in _INTENSITY_STYLE:
                        intensity = IntensityEnum.CONTROL.value
                except (TypeError, ValueError):
                    pass

            counts[name][intensity] += 1

    return counts


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_label_distribution(counts, label, top_n, output_path, prefix=""):
    """Standard horizontal bar chart for a single label category."""
    top_counts = counts.most_common(top_n)
    if not top_counts:
        print(f"  [SKIP] No data for label '{label}'.")
        return

    terms = [t[0] for t in top_counts][::-1]
    frequencies = [t[1] for t in top_counts][::-1]

    fig_height = max(5, len(terms) * 0.45 + 1.5)
    fig, ax = plt.subplots(figsize=(14, fig_height))

    colors = [
        "#cccccc" if any(w in t.lower() for w in ("unspecified", "unknown")) else "#4C72B0"
        for t in terms
    ]

    ax.barh(terms, frequencies, color=colors, edgecolor="none", alpha=0.85)

    x_pad = max(frequencies) * 0.012
    for i, freq in enumerate(frequencies):
        ax.text(
            freq + x_pad,
            i,
            str(freq),
            va="center",
            ha="left",
            fontsize=12,
            color="#444444",
        )
    prefix_str = f"{prefix} -- " if prefix else ""
    ax.set_title(f"{prefix_str}Distribution of {label.replace('_', ' ').capitalize()} labels",
                 fontsize=18, pad=15)
    ax.set_xlabel("Number of samples", fontsize=11)
    ax.set_ylabel(f"{label.replace('_', ' ').capitalize()} term", fontsize=11)
    ax.set_xlim(right=max(frequencies) * 1.12)
    
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", labelsize=9)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved -> {output_path}")


def plot_treatment_intensity(treatment_counts, top_n, output_path, prefix=""):
    """
    Horizontal stacked bar chart: each bar = one treatment type, split by intensity.
      green  = Control (0)
      orange = Mild (1)
      red    = Moderate (2)
      dark red = Severe (3)
    Bars are sorted by total sample count, highest at the top.
    """
    totals = {name: sum(ic.values()) for name, ic in treatment_counts.items()}
    sorted_names = sorted(totals, key=totals.get, reverse=True)[:top_n]
    sorted_names = sorted_names[::-1]  # reverse so highest bar sits at top

    if not sorted_names:
        print("  [SKIP] No treatment data found.")
        return

    max_total = max(totals[n] for n in sorted_names)
    fig_height = max(5, len(sorted_names) * 0.45 + 1.5)
    fig, ax = plt.subplots(figsize=(14, fig_height))

    lefts = [0] * len(sorted_names)
    legend_patches = []

    for intensity in _INTENSITY_ORDER:
        label_str, color = _INTENSITY_STYLE[intensity]
        values = [treatment_counts[name].get(intensity, 0) for name in sorted_names]

        bars = ax.barh(sorted_names, values, left=lefts, color=color,
                       edgecolor="none", alpha=0.88)

        # White count label inside segment — only if segment is wide enough to read
        min_width_for_label = max_total * 0.04
        for bar, val in zip(bars, values):
            if val > 0 and bar.get_width() >= min_width_for_label:
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_y() + bar.get_height() / 2,
                    str(val), va="center", ha="center",
                    fontsize=11, color="white", fontweight="bold",
                )

        lefts = [left + val for left, val in zip(lefts, values)]
        legend_patches.append(mpatches.Patch(color=color, label=label_str))

    # Total count label to the right of each full bar
    x_pad = max_total * 0.012
    for i, name in enumerate(sorted_names):
        ax.text(totals[name] + x_pad, i, str(totals[name]),
                va="center", ha="left", fontsize=12, color="#444444")

    prefix_str = f"{prefix} -- " if prefix else ""
    ax.set_title(f"{prefix_str}Treatment distribution by intensity", fontsize=18, pad=15)
    ax.set_xlabel("Number of samples", fontsize=16)
    ax.set_ylabel("Treatment", fontsize=16)
    ax.set_xlim(right=max_total * 1.14)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", labelsize=13)
    ax.legend(handles=legend_patches, loc="lower right", fontsize=12,
              title="Intensity", title_fontsize=12, framealpha=0.7)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved -> {output_path}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print(f"Loading data from: {LABELS_PATH}")
    samples = load_all_labels(LABELS_PATH)
    print(f"\nTotal samples pooled: {len(samples)}")

    # 1. Get valid IDs from the actual expression data
    data_cols = pd.read_csv(EXPR_PATH, index_col=0, nrows=0).columns
    
    if RNA_USED:
        valid_data_ids = set([get_gsm_id(col.split('_')[-1]) for col in data_cols])
    else:
        valid_data_ids = set(data_cols)

    # 2. Extract IDs, Filter for data presence, and Drop Duplicates
    label_id_key = 'Sample_ID' 
    filtered_samples = []
    seen_ids = set()
    duplicate_count = 0

    for s in samples:
        sid = s.get(label_id_key)
        
        # Condition 1: Must be in the expression matrix
        if sid in valid_data_ids:
            # Condition 2: Must not have been added already
            if sid not in seen_ids:
                filtered_samples.append(s)
                seen_ids.add(sid)
            else:
                duplicate_count += 1
    
    # 3. Validation Logic (Reporting Mismatches)
    pooled_label_ids = set([s.get(label_id_key) for s in samples])
    overlap = seen_ids # After loop, seen_ids is the intersection of unique labels and data
    
    print("--- FLORA Label Synchronization Report ---")
    print(f"  [1] Expression IDs found:             {len(valid_data_ids)}")
    print(f"  [1] Unique Labels found in JSONs:     {len(pooled_label_ids)}")
    print(f"  [1] Duplicate Labels dropped:         {duplicate_count}")
    print(f"  [1] Confirmed Unique Overlap:         {len(overlap)}")
    
    # Report mismatches
    in_labels_not_data = pooled_label_ids - valid_data_ids
    in_data_not_labels = valid_data_ids - pooled_label_ids
    
    if in_labels_not_data:
        print(f"  [!] Labels UNUSED (not in data):     {len(in_labels_not_data)}")
    if in_data_not_labels:
        print(f"  [!] Data UNLABELED (missing JSON):   {len(in_data_not_labels)}")

    samples = filtered_samples
    print(f"Final unique samples retained: {len(samples)}\n")
    if not samples:
        print("No samples found. Exiting.")
        sys.exit(0)

    prefix = "RNA" if RNA_USED else "Microarray"
    output_dir = os.path.join(FIGURES_DIR, f"label_distributions_{prefix}")
    os.makedirs(output_dir, exist_ok=True)

    print(f"Saving plots to: {output_dir}\n")

    for label in LABELS:
        save_path = os.path.join(output_dir, f"{prefix}_{label}_distribution.pdf")

        if label == "treatment":
            # Treatment gets a stacked-intensity chart instead of a plain bar chart
            treatment_counts = aggregate_treatment_intensity_counts(samples)
            plot_treatment_intensity(treatment_counts, top_n=20,
                                     output_path=save_path, prefix=prefix)
        else:
            counts = aggregate_label_counts(samples, label)
            plot_label_distribution(counts, label, top_n=20,
                                    output_path=save_path, prefix=prefix)

    print("\nDone.")