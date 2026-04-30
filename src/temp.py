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
    # Basic Stats
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
    
    # 1. Distribution Comparison (Density Plot)
    plt.figure(figsize=(10, 6))
    for d in data_dicts:
        # We sample 100k points to keep it fast
        vals = d['df'].values.flatten()
        sampled_vals = np.random.choice(vals, min(100000, len(vals)), replace=False)
        sns.kdeplot(sampled_vals, label=f"{d['name']} (Max: {d['max']:.2f})", fill=True, alpha=0.3)
    
    plt.title("Expression Value Distributions")
    plt.xlabel("Value")
    plt.ylabel("Density")
    plt.legend()
    plt.savefig(f"{output_dir}/01_distributions.png")
    plt.close()

    # 2. Mean vs Standard Deviation (Biological Signal check)
    fig, axes = plt.subplots(1, len(data_dicts), figsize=(18, 5), sharey=True)
    for i, d in enumerate(data_dicts):
        means = d['df'].mean(axis=0)
        stds = d['df'].std(axis=0)
        axes[i].scatter(means, stds, s=1, alpha=0.2, color='teal')
        axes[i].set_title(f"Mean-Var: {d['name']}")
        axes[i].set_xlabel("Mean Expression")
        if i == 0: axes[i].set_ylabel("Std Dev")
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/02_mean_variance.png")
    plt.close()

    # 3. PCA (Global Structure / Batch Effects)
    plt.figure(figsize=(10, 7))
    combined_data = []
    labels = []
    
    # We find common genes to make PCA comparable
    common_genes = set.intersection(*[set(d['df'].columns) for d in data_dicts])
    print(f"Comparing across {len(common_genes)} common genes...")
    
    for d in data_dicts:
        subset = d['df'][list(common_genes)]
        # Scale for PCA
        scaled = StandardScaler().fit_transform(subset)
        pca = PCA(n_components=2)
        coords = pca.fit_transform(scaled)
        
        plt.scatter(coords[:, 0], coords[:, 1], label=d['name'], alpha=0.6, s=20)
    
    plt.title("PCA Projection: Global Structural Differences")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.legend()
    plt.savefig(f"{output_dir}/03_pca_comparison.png")
    plt.close()

    # 4. Sparsity (Zero counts)
    plt.figure(figsize=(8, 6))
    names = [d['name'] for d in data_dicts]
    sparsities = [d['sparsity'] * 100 for d in data_dicts]
    sns.barplot(x=names, y=sparsities, palette="viridis")
    plt.title("Percentage of Zero Values (Sparsity)")
    plt.ylabel("% Zeros")
    plt.savefig(f"{output_dir}/04_sparsity.png")
    plt.close()

    print(f"All plots saved to {output_dir}/")
if __name__ == "__main__":
    # Update these paths to your actual locations
    FOLDER_A = f"{STORAGE_DIR}final_data/rnaseq_processed"
    FOLDER_B = f"{STORAGE_DIR}final_data/"

    # Define your specific files
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