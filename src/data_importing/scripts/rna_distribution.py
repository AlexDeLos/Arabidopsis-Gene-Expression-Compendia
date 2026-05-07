import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

module_dir = "./"
sys.path.append(module_dir)

# Import the new map alongside the others
from src.constants import FIGURES_DIR, STORAGE_DIR  # noqa: E402

# 1. Configuration
INPUT_FILE = f"{STORAGE_DIR}/final_data/rnaseq_processed/filter_norm.csv"

ext = 'TMP_log_norm'
def read_id(path):
    with open(path) as f:
        return f.read().strip()


def evaluate_distributions():
    OUTPUT_PLOT = f"{FIGURES_DIR}distribution_figures/mathematical_study_evaluation_{ext}.svg"
    if not os.path.exists(INPUT_FILE):
        print(f"Error: Could not find {INPUT_FILE}")
        return

    print(f"Loading data from {INPUT_FILE}...")
    df = pd.read_csv(INPUT_FILE, index_col=0)

    try:
        rnaseq_ids = eval("['" + read_id("./study_ids/RNA_seq_ids.txt").replace(",", "','") + "']")
        valid_columns = [col for col in df.columns if col.split("_")[0] in rnaseq_ids]
        df = df[valid_columns]
    except Exception as e:
        print(f"Could not filter by ID list (using all columns): {e}")

    df = df.select_dtypes(include=[np.number])
    print(f"Loaded matrix with {df.shape[0]} genes and {df.shape[1]} samples.")

    # Extract Study IDs from column names (Assuming format GSEXXXXX_SRRYYYYY)
    study_ids = np.array([col.split("_")[0] for col in df.columns])
    unique_studies = np.unique(study_ids)
    print(f"Found {len(unique_studies)} unique studies.")

    # 2. Log Transformation
    print("Applying log1p transformation...")
    df_log = np.log1p(df)

    # 3. Mathematical Evaluation
    print("Calculating quantiles and mathematical distances...")
    percentiles = np.linspace(0, 100, 100)
    quantile_matrix = np.zeros((df_log.shape[1], len(percentiles)))

    # Also prepare the standard histogram for the density plot
    max_val = df_log.max().max()
    bins = np.linspace(0, max_val, 100)
    bin_centers = 0.5 * (bins[1:] + bins[:-1])
    hist_matrix = np.zeros((df_log.shape[1], len(bins) - 1))

    for i, col in enumerate(df_log.columns):
        series = df_log[col].dropna()
        quantile_matrix[i, :] = np.percentile(series, percentiles)
        hist_matrix[i, :], _ = np.histogram(series, bins=bins, density=True)

    # Reference Distribution
    reference_quantiles = np.median(quantile_matrix, axis=0)
    median_hist = np.median(hist_matrix, axis=0)

    p05_quantiles = np.percentile(quantile_matrix, 5, axis=0)
    p95_quantiles = np.percentile(quantile_matrix, 95, axis=0)
    p05_hist = np.percentile(hist_matrix, 5, axis=0)
    p95_hist = np.percentile(hist_matrix, 95, axis=0)

    # --- NEW: Study-Level Aggregation ---
    print("Aggregating deviations by study mean...")
    distances = np.mean(np.abs(quantile_matrix - reference_quantiles), axis=1)

    df_dist = pd.DataFrame({"Sample": df_log.columns, "Distance": distances, "Study": study_ids})

    # Calculate the mean distance for each study
    study_stats = df_dist.groupby("Study")["Distance"].agg(["mean", "count"]).sort_values(by="mean", ascending=False)

    # Identify top 3 worst studies for curve highlighting
    top_worst_studies = study_stats.head(3).index.tolist()

    # Calculate the "Mean Curve" for those worst studies
    study_mean_quantiles = {}
    study_mean_hists = {}
    for study in top_worst_studies:
        s_idx = np.where(study_ids == study)[0]
        study_mean_quantiles[study] = np.mean(quantile_matrix[s_idx, :], axis=0)
        study_mean_hists[study] = np.mean(hist_matrix[s_idx, :], axis=0)

    # 4. Plotting a 4-Panel Dashboard
    print("Generating study-level dashboard...")
    _fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    ax1, ax2, ax3, ax4 = axes.flatten()

    colors = ["#e63946", "#f4a261", "#e9c46a"]  # Distinct colors for the 3 worst studies

    # --- PANEL 1: Q-Q Plot (Highlighting Study Means) ---
    ax1.fill_between(reference_quantiles, p05_quantiles, p95_quantiles, color="#4c72b0", alpha=0.3, label="90% Sample Range")
    ax1.plot(reference_quantiles, reference_quantiles, color="black", linestyle="--", linewidth=2, label="Perfect Match (y=x)")

    for i, study in enumerate(top_worst_studies):
        ax1.plot(reference_quantiles, study_mean_quantiles[study], color=colors[i], linewidth=2.5, label=f"{study} (Avg of {study_stats.loc[study, 'count']} samples)")

    ax1.set_title("Q-Q Plot: Worst Studies vs Reference", fontsize=14)
    ax1.set_xlabel("Reference Quantiles (Median)", fontsize=12)
    ax1.set_ylabel("Study Mean Quantiles", fontsize=12)
    ax1.legend(loc="upper left")
    ax1.grid(True, alpha=0.3)

    # --- PANEL 2: Histogram of All Sample Distances ---
    ax2.hist(distances, bins=50, color="#003366", edgecolor="black", alpha=0.7)
    ax2.set_title("Overall Sample Deviation Distribution", fontsize=14)
    ax2.set_xlabel("Sample Distance from Median", fontsize=12)
    ax2.set_ylabel("Number of Samples", fontsize=12)
    ax2.grid(True, alpha=0.3)

    # --- PANEL 3: NEW! Bar Chart of Worst Studies ---
    top_n_bars = min(15, len(study_stats))
    subset_stats = study_stats.head(top_n_bars)

    # Highlight the top 3 in the same colors as the curves, rest in grey
    bar_colors = [colors[i] if i < 3 else "#b0b0b0" for i in range(top_n_bars)]

    ax3.bar(subset_stats.index, subset_stats["mean"], color=bar_colors, edgecolor="black")
    ax3.set_title(f"Top {top_n_bars} Most Deviant Studies (Batch Effects)", fontsize=14)
    ax3.set_ylabel("Mean Distance from Reference", fontsize=12)
    ax3.tick_params(axis="x", rotation=45)
    ax3.grid(True, axis="y", alpha=0.3)

    # --- PANEL 4: Density Plot (Highlighting Study Means) ---
    ax4.fill_between(bin_centers, p05_hist, p95_hist, color="#4c72b0", alpha=0.3)
    ax4.plot(bin_centers, median_hist, color="#003366", linewidth=2.5, label="Global Median Shape")

    for i, study in enumerate(top_worst_studies):
        ax4.plot(bin_centers, study_mean_hists[study], color=colors[i], linewidth=2.5, label=f"{study} Mean Density")

    ax4.set_title("Expression Density: Worst Studies", fontsize=14)
    ax4.set_xlabel("Expression Level (log1p)", fontsize=12)
    ax4.set_ylabel("Density", fontsize=12)
    ax4.legend()
    ax4.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(OUTPUT_PLOT, dpi=300)
    print(f"Plot successfully saved to {OUTPUT_PLOT}")


def evaluate_distributions_samples():
    OUTPUT_PLOT = f"{FIGURES_DIR}distribution_figures/mathematical_sample_evaluation_{ext}.svg"
    os.makedirs(f"{FIGURES_DIR}distribution_figures/",exist_ok=True)
    if not os.path.exists(INPUT_FILE):
        print(f"Error: Could not find {INPUT_FILE}")
        return

    print(f"Loading data from {INPUT_FILE}...")
    df = pd.read_csv(INPUT_FILE, index_col=0)
    try:
        rnaseq_ids = eval("['" + read_id("./study_ids/RNA_seq_ids.txt").replace(",", "','") + "']")
        valid_columns = [col for col in df.columns if col.split("_")[0] in rnaseq_ids]
        df = df[valid_columns]
    except Exception as e:
        print(f"Could not filter by ID list (using all columns): {e}")

    df = df.select_dtypes(include=[np.number])
    print(f"Loaded matrix with {df.shape[0]} genes and {df.shape[1]} samples.")

    # 2. Log Transformation
    print("Applying log1p transformation...")
    df_log = np.log1p(df)

    # 3. Mathematical Evaluation using Quantiles
    print("Calculating quantiles and mathematical distances...")
    # Calculate percentiles from 0 to 100 for every sample
    percentiles = np.linspace(0, 100, 100)

    # Create an array to hold the quantiles: shape (n_samples, 100)
    quantile_matrix = np.zeros((df_log.shape[1], len(percentiles)))

    for i, col in enumerate(df_log.columns):
        quantile_matrix[i, :] = np.percentile(df_log[col].dropna(), percentiles)

    # The "Reference Distribution" is the median of all sample quantiles
    reference_quantiles = np.median(quantile_matrix, axis=0)

    # Calculate distance: Mean Absolute Error between sample quantiles and reference
    # This is an approximation of the 1-Wasserstein Distance for 1D distributions
    distances = np.mean(np.abs(quantile_matrix - reference_quantiles), axis=1)

    # Identify outliers mathematically
    num_outliers = min(5, df_log.shape[1])
    outlier_indices = np.argsort(distances)[-num_outliers:]

    # Calculate bands for plotting
    p05_quantiles = np.percentile(quantile_matrix, 5, axis=0)
    p95_quantiles = np.percentile(quantile_matrix, 95, axis=0)

    # 4. Plotting a 3-Panel Dashboard
    print("Generating dashboard...")
    _fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # --- PANEL 1: Q-Q Plot ---
    ax = axes[0]
    ax.fill_between(reference_quantiles, p05_quantiles, p95_quantiles, color="#4c72b0", alpha=0.3, label="90% Sample Range")
    ax.plot(reference_quantiles, reference_quantiles, color="black", linestyle="--", linewidth=2, label="Perfect Match (y=x)")

    # Plot mathematical outliers
    for i, idx in enumerate(outlier_indices):
        ax.plot(reference_quantiles, quantile_matrix[idx, :], color="red", alpha=0.8, label=f"Outlier: {df_log.columns[idx]}" if i < 3 else None)

    ax.set_title("Quantile-Quantile (Q-Q) Plot", fontsize=14)
    ax.set_xlabel("Reference Quantiles (Median)", fontsize=12)
    ax.set_ylabel("Sample Quantiles", fontsize=12)
    ax.legend(loc="upper left")
    ax.grid(True, alpha=0.3)

    # --- PANEL 2: Distribution of Distances ---
    ax = axes[1]
    ax.hist(distances, bins=50, color="#003366", edgecolor="black", alpha=0.7)

    # Mark the cutoff for outliers
    outlier_threshold = distances[outlier_indices[0]]  # Lowest distance among the top outliers
    ax.axvline(outlier_threshold, color="red", linestyle="--", label=f"Top {num_outliers} Outlier Threshold")

    ax.set_title("Sample Deviation from Reference", fontsize=14)
    ax.set_xlabel("Distance from Median Distribution", fontsize=12)
    ax.set_ylabel("Number of Samples", fontsize=12)
    ax.legend()
    ax.grid(True, alpha=0.3)

    # --- PANEL 3: Standard Density (for intuition) ---
    ax = axes[2]
    # Re-calculate standard histogram for plotting
    max_val = df_log.max().max()
    bins = np.linspace(0, max_val, 100)
    bin_centers = 0.5 * (bins[1:] + bins[:-1])
    hist_matrix = np.zeros((df_log.shape[1], len(bins) - 1))

    for i, col in enumerate(df_log.columns):
        hist_matrix[i, :], _ = np.histogram(df_log[col].dropna(), bins=bins, density=True)

    median_hist = np.median(hist_matrix, axis=0)
    p05_hist = np.percentile(hist_matrix, 5, axis=0)
    p95_hist = np.percentile(hist_matrix, 95, axis=0)

    ax.fill_between(bin_centers, p05_hist, p95_hist, color="#4c72b0", alpha=0.3)
    ax.plot(bin_centers, median_hist, color="#003366", linewidth=2)

    # Plot same outliers in red for visual connection
    for idx in outlier_indices:
        ax.plot(bin_centers, hist_matrix[idx, :], color="red", alpha=0.8, linewidth=1.5)

    ax.set_title("Expression Density Overlay", fontsize=14)
    ax.set_xlabel("Expression Level (log1p)", fontsize=12)
    ax.set_ylabel("Density", fontsize=12)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(OUTPUT_PLOT, dpi=300)
    print(f"Plot successfully saved to {OUTPUT_PLOT}")


if __name__ == "__main__":
    evaluate_distributions_samples()
    evaluate_distributions()
