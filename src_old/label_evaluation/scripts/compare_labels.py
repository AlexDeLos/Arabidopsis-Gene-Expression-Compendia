import os
import json
import argparse
from collections import Counter
from pathlib import Path
import os
import json
from collections import Counter, defaultdict
from pathlib import Path
import pandas as pd

def canonicalize(v):
    if isinstance(v, list):
        if len(v) > 0 and isinstance(v[0], dict):
            return tuple(sorted([tuple(sorted(d.items())) for d in v]))
        return tuple(sorted([str(i).lower().strip() for i in v]))
    return str(v).lower().strip()

def analyze_divergence_patterns(dir1, dir2):
    path1, path2 = Path(dir1), Path(dir2)
    common_files = set(f.name for f in path1.glob("*.json")) & set(f.name for f in path2.glob("*.json"))
    
    study_stats = []
    global_swaps = defaultdict(Counter) # Tracks exactly what changed (e.g., 'root' -> 'roots')

    for filename in common_files:
        with open(path1 / filename) as f1, open(path2 / filename) as f2:
            d1, d2 = json.load(f1), json.load(f2)
        
        gse_id = filename.replace('.json', '')
        common_gsms = set(d1.keys()) & set(d2.keys())
        
        study_diff_count = 0
        study_cat_counts = Counter()

        for gsm in common_gsms:
            sample_diff = False
            for cat in ['tissue', 'treatment']: # Focusing on your priorities
                v1, v2 = canonicalize(d1[gsm].get(cat)), canonicalize(d2[gsm].get(cat))
                if v1 != v2:
                    study_cat_counts[cat] += 1
                    global_swaps[cat][(v1, v2)] += 1
                    sample_diff = True
            if sample_diff: study_diff_count += 1

        if common_gsms:
            study_stats.append({
                'GSE': gse_id,
                'Samples': len(common_gsms),
                'Diff_Samples': study_diff_count,
                'Conflict_Rate': study_diff_count / len(common_gsms),
                'Tissue_Diffs': study_cat_counts['tissue'],
                'Treatment_Diffs': study_cat_counts['treatment']
            })

    df = pd.DataFrame(study_stats)
    
    # 1. Study Concentration Analysis
    print("\n--- TOP 5 MOST CONFLICTED STUDIES ---")
    print(df.sort_values('Conflict_Rate', ascending=False).head(5)[['GSE', 'Samples', 'Conflict_Rate']])

    # 2. Value Swap Analysis (The "Nature" of the change)
    for cat in ['tissue', 'treatment']:
        print(f"\n--- TOP 3 COMMON {cat.upper()} DIVERGENCES ---")
        for (v1, v2), count in global_swaps[cat].most_common(3):
            print(f"  {v1}  ==>  {v2}  ({count} occurrences)")

    # 3. Randomness Test
    high_conflict_studies = len(df[df['Conflict_Rate'] > 0.5])
    zero_conflict_studies = len(df[df['Conflict_Rate'] == 0])
    print(f"\n--- DISTRIBUTION SUMMARY ---")
    print(f"Studies with >50% conflict: {high_conflict_studies}")
    print(f"Studies with 0% conflict  : {zero_conflict_studies}")

def canonicalize(value):
    """
    Converts lists and dictionaries into a stable, sortable format for comparison.
    This handles the [{'val': 'X', 'intensity': 0}] structure vs simple ['X'].
    """
    if isinstance(value, list):
        # If it's a list of dicts (like treatment), convert dicts to sorted tuples
        if len(value) > 0 and isinstance(value[0], dict):
            return tuple(sorted([tuple(sorted(d.items())) for d in value]))
        # Standard list of strings
        return tuple(sorted(value))
    return value

def compare_labels(dir1, dir2):
    path1 = Path(dir1)
    path2 = Path(dir2)

    # Find common filenames
    files1 = set(f.name for f in path1.glob("*.json"))
    files2 = set(f.name for f in path2.glob("*.json"))
    common_files = files1.intersection(files2)

    print(f"Found {len(common_files)} common files between directories.")
    
    total_samples_compared = 0
    total_differences = 0
    category_diffs = Counter()
    
    for filename in common_files:
        with open(path1 / filename, 'r') as f1, open(path2 / filename, 'r') as f2:
            data1 = json.load(f1)
            data2 = json.load(f2)
            
        # Only compare GSMs present in both files
        common_gsms = set(data1.keys()).intersection(data2.keys())
        
        for gsm in common_gsms:
            total_samples_compared += 1
            sample_has_diff = False
            
            # Categories are the keys inside the GSM dict (tissue, treatment, etc.)
            categories = set(data1[gsm].keys()).union(data2[gsm].keys())
            
            for cat in categories:
                val1 = canonicalize(data1[gsm].get(cat))
                val2 = canonicalize(data2[gsm].get(cat))
                
                if val1 != val2:
                    if val1 == 'unspecified':
                        pass
                    category_diffs[cat] += 1
                    sample_has_diff = True
            
            if sample_has_diff:
                total_differences += 1

    # --- Output Statistics ---
    print("\n" + "="*40)
    print("LABEL COMPARISON STATISTICS")
    print("="*40)
    print(f"Total Shared Files Processed: {len(common_files)}")
    print(f"Total Shared Samples (GSMs): {total_samples_compared}")
    print(f"Samples with at least one difference: {total_differences}")
    
    if total_samples_compared > 0:
        error_rate = (total_differences / total_samples_compared) * 100
        print(f"Overall Divergence Rate: {error_rate:.2f}%")
        
        print("\nDifferences by Category (Ranked):")
        for cat, count in category_diffs.most_common():
            cat_perc = (count / total_samples_compared) * 100
            print(f" - {cat:20}: {count:5} differences ({cat_perc:.2f}%)")
    else:
        print("No overlapping samples found to compare.")

if __name__ == "__main__":
    # parser = argparse.ArgumentParser(description="Compare label JSONs in two directories.")
    # parser.add_argument("dir1", help="Path to the first directory")
    # parser.add_argument("dir2", help="Path to the second directory")
    
    # args = parser.parse_args()
    compare_labels('new_storage/labels/TULIP_1.2_old/5.0', 'new_storage/labels/TULIP_1.2/5.0')
    analyze_divergence_patterns('new_storage/labels/TULIP_1.2_old/5.0', 'new_storage/labels/TULIP_1.2/5.0')