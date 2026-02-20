import os
import json
import matplotlib.pyplot as plt
from collections import Counter
import sys

# Ensure the script can find your src modules if run from the root directory
sys.path.append(os.path.abspath('./'))
from src.constants import LABELS_PATH, LABELS, FIGURES_DIR

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
                    # The JSON structure is { "sample_id": { "tissue": [...], ... } }
                    for sample_id, sample_data in study_data.items():
                        all_samples.append(sample_data)
                except Exception as e:
                    print(f"Error reading {file_path}: {e}")
                    
    return all_samples

def aggregate_label_counts(samples, label_category):
    """Counts the occurrences of each grounded term in a specific category."""
    counter = Counter()
    for sample in samples:
        # Values are lists (e.g., ["Leaf", "Stem"] or ["Control"])
        items = sample.get(label_category, ["unspecified"])
        
        # Handle cases where the item might be a string instead of a list
        if isinstance(items, str):
            items = [items]
            
        for item in items:
            counter[item] += 1
            
    return counter

def plot_label_distributions(samples, output_dir):
    """Generates and saves a bar chart for each label category."""
    os.makedirs(output_dir, exist_ok=True)
    
    for label in LABELS:
        counts = aggregate_label_counts(samples, label)
        
        if not counts:
            print(f"No data found for {label}. Skipping plot.")
            continue
            
        # Sort counts descending
        sorted_items = sorted(counts.items(), key=lambda x: x[1])
        terms = [x[0] for x in sorted_items]
        frequencies = [x[1] for x in sorted_items]
        
        # Create Plot
        fig, ax = plt.subplots(figsize=(10, max(6, len(terms) * 0.3))) # Dynamic height
        
        # Color "unspecified" and "unknown" differently (e.g., grey)
        colors = ['#cccccc' if t in ['unspecified', 'unknown'] else '#4C72B0' for t in terms]
        
        bars = ax.barh(terms, frequencies, color=colors, edgecolor='black', alpha=0.8)
        
        # Add values on the end of the bars
        for bar in bars:
            ax.text(bar.get_width() + (max(frequencies)*0.01), 
                    bar.get_y() + bar.get_height()/2, 
                    f'{int(bar.get_width())}', 
                    va='center', ha='left', fontsize=10)

        ax.set_title(f'Distribution of Grounded {label.capitalize()} Labels', fontsize=14, pad=15)
        ax.set_xlabel('Number of Samples', fontsize=12)
        ax.set_ylabel(f'{label.capitalize()} Terms', fontsize=12)
        
        # Formatting
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        plt.tight_layout()
        
        # Save figure
        save_path = os.path.join(output_dir, f'{label}_distribution.png')
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved plot: {save_path}")

if __name__ == '__main__':
    print(f"Loading data from: {LABELS_PATH}")
    samples = load_all_labels(LABELS_PATH)
    
    if samples:
        print(f"Loaded {len(samples)} total samples. Generating plots...")
        plot_label_distributions(samples, output_dir=FIGURES_DIR+'vis/')
        print("Done!")
    else:
        print("No samples loaded. Check your LABELS_PATH.")