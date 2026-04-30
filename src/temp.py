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
    print(f"Loading {name} from {path}...")
    df = pd.read_csv(path, index_col=0)
    print(f"path {path}:, {df.head}")
    # --- FIX 1: Standardize Column Names ---
    # Convert all gene names to uppercase and strip spaces to ensure matching
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
    return summary

def plot_comparisons(data_dicts, output_dir="plots_comparison"):
    os.makedirs(output_dir, exist_ok=True)
    sns.set_theme(style="whitegrid")
    num_files = len(data_dicts)
    
    # 1. Distribution Comparison
    plt.figure(figsize=(10, 6))
    for d in data_dicts:
        vals = d['df'].values.flatten()
        sampled_vals = np.random.choice(vals, min(100000, len(vals)), replace=False)
        sns.kdeplot(sampled_vals, label=f"{d['name']} (Max: {d['max']:.2f})", fill=True, alpha=0.3)
    plt.title("Expression Value Distributions")
    plt.legend()
    plt.savefig(f"{output_dir}/01_distributions.png")
    plt.close()

    # 2. Mean vs Standard Deviation (with squeeze fix)
    fig, axes = plt.subplots(1, num_files, figsize=(6 * num_files, 5), squeeze=False)
    for i, d in enumerate(data_dicts):
        means = d['df'].mean(axis=0)
        stds = d['df'].std(axis=0)
        axes[0, i].scatter(means, stds, s=1, alpha=0.2, color='teal')
        axes[0, i].set_title(f"Mean-Var: {d['name']}")
        axes[0, i].set_xlabel("Mean Expression")
        if i == 0: axes[0, i].set_ylabel("Std Dev")
    plt.tight_layout()
    plt.savefig(f"{output_dir}/02_mean_variance.png")
    plt.close()

    # 3. PCA (Global Structure)
    # --- FIX 2: Better Intersection Logic ---
    common_genes = list(set.intersection(*[set(d['df'].columns) for d in data_dicts]))
    print(f"Comparing across {len(common_genes)} common genes...")
    
    if len(common_genes) > 10: # Need at least a few genes for PCA
        plt.figure(figsize=(10, 7))
        for d in data_dicts:
            subset = d['df'][common_genes]
            scaled = StandardScaler().fit_transform(subset)
            pca = PCA(n_components=2)
            coords = pca.fit_transform(scaled)
            plt.scatter(coords[:, 0], coords[:, 1], label=d['name'], alpha=0.6, s=20)
        plt.title("PCA Projection")
        plt.legend()
        plt.savefig(f"{output_dir}/03_pca_comparison.png")
        plt.close()
    else:
        print("!!! WARNING: Skipping PCA. Common genes found <= 10. Check gene ID formats.")

    # 4. Sparsity
    plt.figure(figsize=(8, 6))
    names = [d['name'] for d in data_dicts]
    sparsities = [d['sparsity'] * 100 for d in data_dicts]
    sns.barplot(x=names, y=sparsities, palette="viridis")
    plt.title("Percentage of Zero Values")
    plt.ylabel("% Zeros")
    plt.savefig(f"{output_dir}/04_sparsity.png")
    plt.close()

if __name__ == "__main__":
    FOLDER_A = f"{STORAGE_DIR}final_data/rnaseq_processed"
    FOLDER_B = f"{STORAGE_DIR}final_data/"

    files = [
        (os.path.join(FOLDER_A, "rankin.csv"), "Rankin"),
        (os.path.join(FOLDER_B, "filter.csv"), "Filter"),
        (os.path.join(FOLDER_B, "combat.csv"), "ComBat")
    ]

    data_summaries = []
    for path, name in files:
        if os.path.exists(path):
            data_summaries.append(load_and_summarize(path, name))
        else:
            print(f"Skipping {name}: File not found at {path}")

    if data_summaries:
        plot_comparisons(data_summaries)