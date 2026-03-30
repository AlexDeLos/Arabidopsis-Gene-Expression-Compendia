"""
rnaseq_norm_and_analysis.py

Normalization pipeline for the combined RNA-seq count matrix.

Column format: {GSE_id}_{SRR_id}   (e.g. GSE201685_SRR18937094)
Row format   : gene_id              (e.g. AT1G01010)

Pipeline stages (each cached to disk, skipped on re-run):
  1. filter.csv          — remove low-coverage genes/samples, zero → NaN
  2. combat_seq.csv      — ComBat-seq batch correction (count-space, negative-binomial)
  3. rankin.csv          — Rank-in normalization (cross-platform integration)

ComBat-seq reference:
  Zhang et al. (2020) Genome Biology 21:257
  https://pmc.ncbi.nlm.nih.gov/articles/PMC7518324/
"""

import os
import sys
import numpy as np
import pandas as pd
from collections import Counter
from sklearn.decomposition import TruncatedSVD
import matplotlib.pyplot as plt

# --- R integration ---
import rpy2.robjects as ro
from rpy2.robjects import numpy2ri, pandas2ri
from rpy2.robjects.packages import importr
from rpy2.robjects.conversion import localconverter

module_dir = './'
sys.path.append(module_dir)
from src.constants import *   # provides STORAGE_DIR, COMBINED_DATA_OUTPUT_FILE, etc.

# ---------------------------------------------------------------------------
# Paths  (edit RNASEQ_DATA_DIR / RNASEQ_FIGURES_DIR in src/constants.py
#         or override here)
# ---------------------------------------------------------------------------
RNASEQ_COMBINED   = os.path.join(STORAGE_DIR, "final_data", "Salmon_RNAseq_Combined.csv")
RNASEQ_DATA_DIR   = os.path.join(STORAGE_DIR, "final_data", "rnaseq_processed")
RNASEQ_FIGURES    = os.path.join(STORAGE_DIR, "figures", "rnaseq_filtering")


# ---------------------------------------------------------------------------
# Helpers — batch label extraction
# ---------------------------------------------------------------------------

def get_gse_from_col(col: str) -> str:
    """
    Extracts the GSE study ID from a column named {GSE_id}_{SRR_id}.
    Works whether the column has already been renamed to a GSM or still
    carries the original GSE_SRR format.
    """
    return col.split('_')[0]


def get_batch_labels(columns) -> list[str]:
    """Returns one batch label (GSE ID) per column."""
    return [get_gse_from_col(c) for c in columns]


# ---------------------------------------------------------------------------
# 1. Filtering
# ---------------------------------------------------------------------------

def plot_filtering_summary(df_before: pd.DataFrame, df_after: pd.DataFrame, output_path: str):
    print("  [Plot] Generating filtering summary...")
    fig = plt.figure(figsize=(16, 12))
    fig.suptitle('RNA-seq Data Filtering Summary', fontsize=20, fontweight='bold')

    gs = fig.add_gridspec(2, 2, height_ratios=[1, 2.5])

    ax_bars = fig.add_subplot(gs[0, :])
    categories = ['Genes (Rows)', 'Samples (Columns)']
    before_counts = [df_before.shape[0], df_before.shape[1]]
    after_counts  = [df_after.shape[0],  df_after.shape[1]]
    x, w = np.arange(2), 0.35
    ax_bars.bar(x - w/2, before_counts, w, label='Before', color='#ff9999', edgecolor='black')
    ax_bars.bar(x + w/2, after_counts,  w, label='After',  color='#66b3ff', edgecolor='black')
    ax_bars.set_xticks(x); ax_bars.set_xticklabels(categories, fontsize=12)
    ax_bars.set_ylabel('Count'); ax_bars.legend()
    ax_bars.set_title('Genes and Samples Before vs After Filtering')
    for i, v in enumerate(before_counts):
        ax_bars.text(i - w/2, v * 1.01, f"{v:,}", ha='center', fontweight='bold')
    for i, v in enumerate(after_counts):
        ax_bars.text(i + w/2, v * 1.01, f"{v:,}", ha='center', fontweight='bold')

    for ax, df, label in [
        (fig.add_subplot(gs[1, 0]), df_before, 'Before'),
        (fig.add_subplot(gs[1, 1]), df_after,  'After'),
    ]:
        ax.imshow(df.isna(), aspect='auto', cmap='viridis', interpolation='nearest')
        pct = (df.isna().sum().sum() / df.size) * 100
        ax.set_title(f'Missingness {label}\n({pct:.1f}% missing)')
        ax.set_xlabel('Samples'); ax.set_ylabel('Genes')

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  [Plot] Saved to {output_path}")


def run_filtering(
    raw_df: pd.DataFrame,
    zero_as_nan: bool = True,
    gene_nan_pct: float = 20.0,
    sample_nan_pct: float = 20.0,
    min_count_threshold: int = 10,
    min_samples_expressed: float = 0.1,
) -> pd.DataFrame:
    """
    RNA-seq-aware filtering:

    1. Genes with total count < min_count_threshold across ALL samples → drop.
       (These are almost certainly un-expressed / mapping artefacts.)
    2. Genes expressed (count > 0) in fewer than min_samples_expressed fraction
       of samples → drop.  (Prevents extremely sparse genes from inflating zeros.)
    3. Zero → NaN  (if zero_as_nan=True), then standard NaN-% thresholds.

    Parameters
    ----------
    raw_df               : genes × samples count matrix (integer or float counts)
    zero_as_nan          : treat 0 as missing after the count-based filters above
    gene_nan_pct         : max % NaN per gene   (row) to keep  [default 20]
    sample_nan_pct       : max % NaN per sample (col) to keep  [default 20]
    min_count_threshold  : drop genes whose row-sum < this value [default 10]
    min_samples_expressed: fraction of samples that must have count > 0 [default 0.1]
    """
    print(f"  [Filter] Input: {raw_df.shape[0]} genes × {raw_df.shape[1]} samples")

    # --- a. Low-count gene filter (count domain, before NaN conversion) ---
    row_sums = raw_df.sum(axis=1)
    raw_df = raw_df.loc[row_sums >= min_count_threshold]
    print(f"  [Filter] After total-count filter (≥{min_count_threshold}): {raw_df.shape[0]} genes")

    # --- b. Low-prevalence gene filter ---
    min_samples = int(np.ceil(min_samples_expressed * raw_df.shape[1]))
    expressed   = (raw_df > 0).sum(axis=1)
    raw_df = raw_df.loc[expressed >= min_samples]
    print(f"  [Filter] After prevalence filter (expressed in ≥{min_samples_expressed*100:.0f}% samples): {raw_df.shape[0]} genes")

    # --- c. Zero → NaN ---
    if zero_as_nan:
        raw_df = raw_df.replace(0, np.nan)

    # --- d. NaN-% filters ---
    nan_genes   = raw_df.isna().mean(axis=1) * 100
    raw_df      = raw_df.loc[nan_genes <= gene_nan_pct]
    nan_samples = raw_df.isna().mean() * 100
    raw_df      = raw_df[raw_df.columns[nan_samples <= sample_nan_pct]]
    print(f"  [Filter] After NaN filters: {raw_df.shape[0]} genes × {raw_df.shape[1]} samples")
    return raw_df


# ---------------------------------------------------------------------------
# 2. ComBat-seq  (count-space batch correction via R's sva package)
# ---------------------------------------------------------------------------

def run_combat_seq(
    counts_df: pd.DataFrame,
    batch_labels: list[str],
    covar_df: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """
    Runs R's sva::ComBat_seq on integer count data.

    Unlike the microarray ComBat (Gaussian), ComBat-seq uses a negative-binomial
    model and operates on RAW counts — do NOT log-transform before passing in.

    Parameters
    ----------
    counts_df    : genes × samples integer/float count matrix (no NaN allowed)
    batch_labels : one batch string per column (e.g. GSE study IDs)
    covar_df     : optional DataFrame of biological covariates (samples as rows)
                   passed as the 'covar_mod' argument to ComBat_seq
    """
    print("  [ComBat-seq] Loading R sva package...")
    try:
        sva  = importr('sva')
        base = importr('base')
    except Exception as e:
        raise ImportError(
            "R package 'sva' not found. Install with:\n"
            "  BiocManager::install('sva')"
        ) from e

    # ComBat-seq requires integer counts — coerce
    counts_int = counts_df.fillna(0).astype(int)
    print(f"  [ComBat-seq] Input: {counts_int.shape[0]} genes × {counts_int.shape[1]} samples")

    # Genes × samples matrix → R integer matrix
    counts_r = ro.r.matrix(
        ro.IntVector(counts_int.values.flatten(order='F')),
        nrow=counts_int.shape[0],
        ncol=counts_int.shape[1],
    )

    batch_r = ro.StrVector(batch_labels)

    # Build covariate model matrix if provided
    covar_mod_r = ro.r('NULL')
    if covar_df is not None and not covar_df.empty:
        print(f"  [ComBat-seq] Building covariate model for: {list(covar_df.columns)}")
        covar_df = covar_df.astype(str)
        stats = importr('stats')
        with localconverter(ro.default_converter + pandas2ri.converter):
            covar_r = ro.conversion.py2rpy(covar_df)
        formula = ro.Formula("~ " + " + ".join(covar_df.columns))
        covar_mod_r = stats.model_matrix(formula, data=covar_r)

    print("  [ComBat-seq] Running sva::ComBat_seq (this may take several minutes)...")
    corrected_r = sva.ComBat_seq(
        counts  = counts_r,
        batch   = batch_r,
        covar_mod = covar_mod_r,
    )

    # Free R memory before allocating Python copy
    del counts_r, batch_r, covar_mod_r
    base.gc()

    corrected_np = np.array(corrected_r)
    del corrected_r
    base.gc()

    result = pd.DataFrame(corrected_np, index=counts_df.index, columns=counts_df.columns)
    print("  [ComBat-seq] Done.")
    return result


# ---------------------------------------------------------------------------
# 3. Rank-in normalization  (same algorithm as microarray pipeline)
# ---------------------------------------------------------------------------

def run_rank_in_normalization(
    df: pd.DataFrame,
    n_bins: int = 100,
    variance_threshold: float = 0.95,
    out_path: str | None = None,
) -> pd.DataFrame:
    """
    Rank-in normalization (Tang et al. 2021).
    Operates on log-transformed data or corrected counts — not raw integers.

    Parameters
    ----------
    df                 : genes × samples expression matrix
    n_bins             : rank bins  (default 100, per paper)
    variance_threshold : cumulative SVD variance to retain  (default 0.95)
    out_path           : if given, saves result to this CSV path
    """
    print(f"\n[Rank-in] Step 1: Intra-sample rank transformation ({n_bins} bins)...")
    rank_matrix  = df.rank(pct=True, method='average')
    binned_matrix = pd.DataFrame(
        np.ceil(rank_matrix.values * n_bins),
        index=df.index, columns=df.columns,
    )

    print("[Rank-in] Step 2: Centering for SVD...")
    gene_means      = binned_matrix.mean(axis=1)
    binned_centered = binned_matrix.sub(gene_means, axis=0)
    X               = binned_centered.T   # (samples × genes) for sklearn

    print("[Rank-in] Step 3: SVD to filter platform noise...")
    max_k    = min(X.shape) - 1
    svd_test = TruncatedSVD(n_components=max_k, random_state=42)
    svd_test.fit(X)
    cumvar   = np.cumsum(svd_test.explained_variance_ratio_)
    optimal_k = int(np.argmax(cumvar >= variance_threshold) + 1)
    print(f"  -> Using {optimal_k} SVD components to retain {variance_threshold*100:.0f}% variance.")

    svd_final     = TruncatedSVD(n_components=optimal_k, random_state=42)
    X_reduced     = svd_final.fit_transform(X)
    X_reconstructed = svd_final.inverse_transform(X_reduced)

    adjusted_df = pd.DataFrame(
        X_reconstructed.T,
        index=binned_matrix.index,
        columns=binned_matrix.columns,
    ).add(gene_means, axis=0)

    print("[Rank-in] Done.")
    if out_path:
        adjusted_df.to_csv(out_path)
        print(f"  -> Saved to {out_path}")
    return adjusted_df


# ---------------------------------------------------------------------------
# 4. Main pipeline
# ---------------------------------------------------------------------------

def run_rnaseq_preprocessing():
    os.makedirs(RNASEQ_DATA_DIR, exist_ok=True)
    os.makedirs(RNASEQ_FIGURES,  exist_ok=True)

    filter_path      = os.path.join(RNASEQ_DATA_DIR, "filter.csv")
    combat_seq_path  = os.path.join(RNASEQ_DATA_DIR, "combat_seq.csv")
    rankin_path      = os.path.join(RNASEQ_DATA_DIR, "rankin.csv")

    # ── Stage 1: Filter ──────────────────────────────────────────────────────
    if os.path.exists(filter_path):
        print("Loading cached filter.csv...")
        filtered_df = pd.read_csv(filter_path, index_col=0)
    else:
        print(f"Loading combined count matrix from {RNASEQ_COMBINED}...")
        raw_df = pd.read_csv(RNASEQ_COMBINED, index_col=0)
        print(f"  Raw matrix: {raw_df.shape[0]} genes × {raw_df.shape[1]} samples")

        filtered_df = run_filtering(raw_df)

        plot_filtering_summary(
            raw_df, filtered_df,
            os.path.join(RNASEQ_FIGURES, "filtering_summary.png"),
        )
        filtered_df.to_csv(filter_path)
        print(f"Saved filtered matrix → {filter_path}")

    # ── Stage 2: ComBat-seq ──────────────────────────────────────────────────
    if os.path.exists(combat_seq_path):
        print("Loading cached combat_seq.csv...")
        combat_df = pd.read_csv(combat_seq_path, index_col=0)
    else:
        print("\nRunning ComBat-seq batch correction...")

        batch_labels = get_batch_labels(filtered_df.columns)
        study_counts = Counter(batch_labels)

        # ComBat-seq requires ≥ 2 samples per batch
        single_batches = {s for s, c in study_counts.items() if c < 2}
        if single_batches:
            print(f"  Removing {len(single_batches)} single-sample batches: {single_batches}")

        valid_cols    = [c for c, b in zip(filtered_df.columns, batch_labels) if b not in single_batches]
        valid_batches = [b for b in batch_labels if b not in single_batches]
        df_for_combat = filtered_df[valid_cols]

        # ComBat-seq needs integer counts and no NaN — fill remaining NaNs with 0
        print(f"  NaN in input: {df_for_combat.isna().sum().sum()} → filling with 0 for ComBat-seq")
        df_for_combat = df_for_combat.fillna(0)

        combat_df = run_combat_seq(df_for_combat, valid_batches)
        combat_df.to_csv(combat_seq_path)
        print(f"Saved ComBat-seq result → {combat_seq_path}")

    # ── Stage 3: Rank-in ─────────────────────────────────────────────────────
    if os.path.exists(rankin_path):
        print("Loading cached rankin.csv...")
        rankin_df = pd.read_csv(rankin_path, index_col=0)
    else:
        print("\nRunning Rank-in normalization on ComBat-seq output...")

        # Log1p transform the corrected counts before Rank-in
        # (Rank-in was designed for log-space expression values)
        log_df = np.log1p(combat_df.clip(lower=0))

        rankin_df = run_rank_in_normalization(
            log_df, n_bins=100, variance_threshold=0.95, out_path=rankin_path
        )
        print(f"Saved Rank-in result → {rankin_path}")

    print("\n=== Pipeline Complete ===")
    print(f"  filter.csv    : {filtered_df.shape}")
    print(f"  combat_seq.csv: {combat_df.shape}")
    print(f"  rankin.csv    : {rankin_df.shape}")
    return filtered_df, combat_df, rankin_df


if __name__ == '__main__':
    run_rnaseq_preprocessing()
