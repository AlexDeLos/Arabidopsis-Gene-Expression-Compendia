import glob
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA

# 1. CRITICAL FOR CLUSTERS: Use non-interactive backend
plt.switch_backend("Agg")


def check_rnaseq_quality(file_path, output_dir, study_id):
    """
    Runs visualization tests and saves output to the centralized output_dir.
    Files will be named: {output_dir}/{study_id}_qc_{test_name}.png
    """
    print(f"  > Processing {study_id}...")

    try:
        # Load Data
        df = pd.read_csv(file_path, sep="\t", index_col=0)

        # Drop 'gene_name' if present (common in nf-core output)
        if "gene_name" in df.columns:
            df = df.drop(columns=["gene_name"])

        if df.empty or df.shape[1] < 2:
            print(f"    [SKIP] {study_id}: Matrix is empty or has too few samples.")
            return

        # Setup output prefix inside the central folder
        # Result: ./outputs/plots/GSE108118_qc
        prefix = os.path.join(output_dir, f"{study_id}_qc")

        # --- PRE-PROCESSING ---
        # Filter noise: genes with > 10 counts in at least 3 samples
        mask = (df > 10).sum(axis=1) >= 3
        df_clean = df[mask]

        if df_clean.shape[0] < 50:
            print(f"    [SKIP] {study_id}: Too few valid genes found ({df_clean.shape[0]}).")
            return

        # Normalize (Log2 + 1)
        log_counts = np.log2(df_clean + 1)

        # --- PLOT 1: Library Sizes ---
        plt.figure(figsize=(10, 6))
        library_sizes = df.sum(axis=0) / 1e6  # Millions
        sns.barplot(x=library_sizes.index, y=library_sizes.values, color="skyblue")
        plt.title(f"{study_id}: Library Sizes")
        plt.ylabel("Reads (Millions)")
        plt.xticks(rotation=90, fontsize=8)
        plt.tight_layout()
        plt.savefig(f"{prefix}_1_lib_size.png")
        plt.close()

        # --- PLOT 2: PCA ---
        if df.shape[1] >= 3:
            pca = PCA(n_components=2)
            pca_result = pca.fit_transform(log_counts.T)

            plt.figure(figsize=(8, 8))
            sns.scatterplot(x=pca_result[:, 0], y=pca_result[:, 1], s=100)

            # Label points
            for i, txt in enumerate(log_counts.columns):
                plt.annotate(txt, (pca_result[i, 0], pca_result[i, 1]), fontsize=8, alpha=0.7)

            plt.title(f"{study_id}: PCA Plot")
            plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}%)")
            plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1] * 100:.1f}%)")
            plt.tight_layout()
            plt.savefig(f"{prefix}_2_PCA.png")
            plt.close()

        # --- PLOT 3: Correlation Heatmap ---
        plt.figure(figsize=(10, 10))
        corr_matrix = log_counts.corr(method="spearman")
        sns.heatmap(corr_matrix, annot=False, cmap="viridis", vmin=0.8, vmax=1.0)
        plt.title(f"{study_id}: Sample Correlation")
        plt.tight_layout()
        plt.savefig(f"{prefix}_3_corr.png")
        plt.close()

    except Exception as e:
        print(f"    [ERROR] Failed to process {study_id}: {e}")


def run_batch_qc(search_root, output_root):
    """
    Scans search_root for data and saves plots to output_root.
    """
    print("--- Starting Batch QC Scan ---")
    print(f"  Scanning: {search_root}")
    print(f"  Saving to: {output_root}")

    # 1. Create the central output directory if it doesn't exist
    os.makedirs(output_root, exist_ok=True)

    # 2. Find files
    search_pattern = os.path.join(search_root, "*", "kallisto", "kallisto.merged.gene_counts.tsv")
    search_pattern_alt = os.path.join(search_root, "*", "salmon", "salmon.merged.gene_counts.tsv")

    files_found = glob.glob(search_pattern) + glob.glob(search_pattern_alt)

    if not files_found:
        print("No count matrix files found.")
        return

    print(f"Found {len(files_found)} studies.")

    # 3. Process
    for file_path in files_found:
        study_dir = os.path.dirname(os.path.dirname(file_path))
        study_id = os.path.basename(study_dir)
        # Pass the central output folder
        check_rnaseq_quality(file_path, output_root, study_id)


if __name__ == "__main__":
    # --- CONFIGURATION ---
    SEARCH_DIR = "/tudelft.net/staff-umbrella/GeneExpressionStorage/rnaseq_data/processed_rnaseq"  # Where to look for GSE folders
    OUTPUT_DIR = "./outputs/plots/"  # Where to save all images
    run_batch_qc(SEARCH_DIR, OUTPUT_DIR)
