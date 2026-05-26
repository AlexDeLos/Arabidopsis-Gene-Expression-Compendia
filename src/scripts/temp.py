import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import sys

sys.path.append(os.path.abspath("./"))
from src.constants import STORAGE_DIR

def load_and_summarize(path, name):
    if not os.path.exists(path):
        print(f"Skipping {name}: File not found at {path}")
        return None
        
    print(f"Loading {name} from {path}...")
    df = pd.read_csv(path, index_col=0)
    
    # --- TRANSPOSE FIX ---
    # PCA and Distributions require Genes as Columns. 
    # Your logs showed [34858 rows x 7582 columns] where rows were genes.
    if str(df.index[0]).upper().startswith('AT'):
        print(f"  -> Detected Genes-as-Rows for {name}. Transposing...")
        df = df.T
    
    # Standardize gene names (columns)
    df.columns = [str(c).upper().strip() for c in df.columns]
    
    summary = {
        "name": name,
        "samples": df.shape[0],
        "genes": df.shape[1],
        "sparsity": (df == 0).sum().sum() / df.size,
        "mean": df.values.mean(),
        "max": df.values.max(),
        "df": df
    }
    print(f"  -> Final Shape: {df.shape[0]} samples, {df.shape[1]} genes")
    return summary

def plot_comparisons(data_dicts, group_name):
    if not data_dicts:
        return
        
    output_dir = f"{STORAGE_DIR}plots_{group_name}"
    os.makedirs(output_dir, exist_ok=True)
    sns.set_theme(style="whitegrid")
    num_files = len(data_dicts)
    
    print(f"\nGenerating plots for {group_name}...")

    # 1. Distribution Comparison
    plt.figure(figsize=(10, 6))
    for d in data_dicts:
        vals = d['df'].values.flatten()
        sampled_vals = np.random.choice(vals, min(100000, len(vals)), replace=False)
        sns.kdeplot(sampled_vals, label=f"{d['name']}", fill=True, alpha=0.3)
    plt.title(f"Distributions - {group_name}")
    plt.legend()
    plt.savefig(f"{output_dir}/01_distributions.svg")
    plt.close()

    # 2. Mean vs Std Dev (Squeeze fix for single-file folders)
    fig, axes = plt.subplots(1, num_files, figsize=(6 * num_files, 5), squeeze=False)
    for i, d in enumerate(data_dicts):
        means = d['df'].mean(axis=0)
        stds = d['df'].std(axis=0)
        axes[0, i].scatter(means, stds, s=1, alpha=0.2, color='teal')
        axes[0, i].set_title(f"Mean-Var: {d['name']}")
        if i == 0: axes[0, i].set_ylabel("Std Dev")
    plt.tight_layout()
    plt.savefig(f"{output_dir}/02_mean_variance.svg")
    plt.close()

    # 3. PCA (Only if there are at least 2 files to compare)
    if num_files > 1:
        common_genes = list(set.intersection(*[set(d['df'].columns) for d in data_dicts]))
        print(f"  -> PCA: Comparing across {len(common_genes)} common genes...")
        
        if len(common_genes) > 10:
            plt.figure(figsize=(10, 7))
            for d in data_dicts:
                subset = d['df'][common_genes]
                scaled = StandardScaler().fit_transform(subset)
                pca = PCA(n_components=2)
                coords = pca.fit_transform(scaled)
                plt.scatter(coords[:, 0], coords[:, 1], label=d['name'], alpha=0.6, s=20)
            plt.title(f"PCA - {group_name}")
            plt.legend()
            plt.savefig(f"{output_dir}/03_pca.svg")
            plt.close()

    # 4. Sparsity
    plt.figure(figsize=(8, 6))
    names = [d['name'] for d in data_dicts]
    sparsities = [d['sparsity'] * 100 for d in data_dicts]
    sns.barplot(x=names, y=sparsities, hue=names, palette="viridis", legend=False)
    plt.title(f"Sparsity % - {group_name}")
    plt.savefig(f"{output_dir}/04_sparsity.svg")
    plt.close()
    
    print(f"Done. Plots saved in {output_dir}/")

if __name__ == "__main__":
    FOLDER_A = f"{STORAGE_DIR}final_data/rnaseq_processed"
    FOLDER_B = f"{STORAGE_DIR}final_data/"

    # Group A Comparison
    group_a_files = [
        (os.path.join(FOLDER_A, "rankin.csv"), "Rankin Norm"),
        (os.path.join(FOLDER_A, "filter_norm.csv"), "Filter Norm"),
        (os.path.join(FOLDER_B, "combat_norm_cov.csv"), "ComBat Cov"),
        (os.path.join(FOLDER_A, "combat_norm.csv"), "ComBat Norm")
    ]
    
    # Group B Comparison
    group_b_files = [
        (os.path.join(FOLDER_B, "rankin.csv"), "Rankin Norm"),
        (os.path.join(FOLDER_B, "filter.csv"), "Filter Norm"),
        (os.path.join(FOLDER_B, "combat_norm.csv"), "ComBat Norm"),
        (os.path.join(FOLDER_B, "combat_norm_cov.csv"), "ComBat Cov")

        # (os.path.join(FOLDER_B, "rankin_old.csv"), "Rankin OLD"),
        # (os.path.join(FOLDER_B, "filter_old.csv"), "Filter OLD"),
    ]

    # Process Group A
    summaries_a = [load_and_summarize(p, n) for p, n in group_a_files]
    summaries_a = [s for s in summaries_a if s is not None]
    plot_comparisons(summaries_a, "Folder_RNA_seq")

    # Process Group B
    summaries_b = [load_and_summarize(p, n) for p, n in group_b_files]
    summaries_b = [s for s in summaries_b if s is not None]
    plot_comparisons(summaries_b, "Folder_Microarray")