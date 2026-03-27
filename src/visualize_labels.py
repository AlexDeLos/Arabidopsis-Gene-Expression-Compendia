import os
import json
import ast
import matplotlib.pyplot as plt
from collections import Counter
import sys

# Ensure the script can find your src modules if run from the root directory
sys.path.append(os.path.abspath('./'))
from src.constants import LABELS_PATH, FIGURES_DIR
from src.constants_labeling import LABELS

def load_all_labels(labels_dir):
    """Loads all sample JSON files from the labels directory."""
    all_samples = []
    if not os.path.exists(labels_dir):
        print(f"Error: Directory {labels_dir} does not exist.")
        return all_samples

    for file in os.listdir(labels_dir):
        if file.endswith('.json') and not file.startswith('map_'):
            file_path = os.path.join(labels_dir, file)
            with open(file_path, 'r') as f:
                try:
                    study_data = json.load(f)
                    for sample_id, sample_data in study_data.items():
                        all_samples.append(sample_data)
                except Exception as e:
                    print(f"Error reading {file_path}: {e}")
                    
    return all_samples

def clean_and_format_str(s):
    """Safely cleans strings and applies Title Case to merge casing fractures."""
    if not s:
        return ""
    s = str(s).strip("[]'\" {}")
    if not s or s.lower() in ["none", "null", "na"]:
        return ""
    return s.title() # Title Case merges "heat stress" and "Heat Stress"

def aggregate_label_counts(samples, label_category):
    """Counts the occurrences of each grounded term in a specific category."""
    counter = Counter()
    
    ignore_words = {"none", "unspecified", "unknown", "na", "n/a", "null", "", "false", "undefined"}
    ignore_keys = {"value", "confidence", "reasoning", "evidence", "explanation", "source"}

    for sample in samples:
        raw_items = sample.get(label_category, [])
        
        # NORMALIZATION: Prevent iterating over dict keys or single strings
        if isinstance(raw_items, dict) or isinstance(raw_items, str):
            items = [raw_items]
        elif raw_items is None:
            items = []
        else:
            items = list(raw_items)
            
        if not items:
            items = ["Unspecified"]
        
        for item in items:
            # Catch stringified JSON dicts from the LLM
            if isinstance(item, str) and item.strip().startswith('{'):
                try:
                    item = ast.literal_eval(item)
                except (ValueError, SyntaxError):
                    pass # Fallback to treating it as a string
            
            if isinstance(item, dict):
                base_val_raw = item.get("value", "Unspecified")
                if isinstance(base_val_raw, list):
                    base_val_raw = base_val_raw[0] if base_val_raw else "Unspecified"
                    
                base_value = clean_and_format_str(base_val_raw)
                
                sub_vals = []
                for key, val in item.items():
                    if key.lower() in ignore_keys:
                        continue
                    
                    if isinstance(val, list):
                        val = ", ".join([str(v) for v in val])
                        
                    clean_val = str(val).strip("[]'\" ")
                    if clean_val.lower() not in ignore_words:
                        sub_vals.append(clean_and_format_str(clean_val))
                
                if sub_vals:
                    formatted_term = f"{base_value} ({', '.join(sub_vals)})"
                else:
                    formatted_term = base_value
                    
                counter.update([formatted_term])
            else:
                # Standard string lists
                clean_item = clean_and_format_str(item)
                counter.update([clean_item if clean_item else "Unspecified"])
                
    return counter

if __name__ == '__main__':
    LABELS_DIR = LABELS_PATH
    
    print(f"Loading data from: {LABELS_DIR}")
    samples = load_all_labels(LABELS_DIR)
    print(f"Loaded {len(samples)} samples.")

    if not samples:
        print("No samples found. Exiting.")
        sys.exit(0)

    output_dir = os.path.join(FIGURES_DIR, "label_distributions")
    os.makedirs(output_dir, exist_ok=True)

    for label in LABELS:
        counts = aggregate_label_counts(samples, label)
        top_counts = counts.most_common(20)
        
        if not top_counts:
            continue
            
        terms = [t[0] for t in top_counts][::-1] 
        frequencies = [t[1] for t in top_counts][::-1]

        fig, ax = plt.subplots(figsize=(10, 8))
        colors = ['#cccccc' if 'unspecified' in t.lower() or 'unknown' in t.lower() else '#4C72B0' for t in terms]
        
        bars = ax.barh(terms, frequencies, color=colors, edgecolor='black', alpha=0.8)
        
        for bar in bars:
            ax.text(bar.get_width() + (max(frequencies)*0.01), 
                    bar.get_y() + bar.get_height()/2, 
                    f'{int(bar.get_width())}', 
                    va='center', ha='left', fontsize=10)

        ax.set_title(f'Distribution of Grounded {label.capitalize()} Labels', fontsize=14, pad=15)
        ax.set_xlabel('Number of Samples', fontsize=12)
        ax.set_ylabel(f'{label.capitalize()} Terms', fontsize=12)
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        plt.tight_layout()
        
        save_path = os.path.join(output_dir, f'{label}_distribution.png')
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved plot: {save_path}")