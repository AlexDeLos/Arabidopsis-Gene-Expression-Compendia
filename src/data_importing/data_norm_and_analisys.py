"""
microarray_norm_and_analysis.py

Normalization pipeline for the combined Microarray expression matrix.

Input assumption:
  - RMA-processed matrix (Robust Multi-Array Average)
  - RMA output is already log2-transformed and quantile-normalized per array
  - No NaNs in the input matrix
  - Column format: {GSE_id}_{GSM_id}   (e.g. GSE12345_GSM123456)
  - Row format   : gene_id              (e.g. AT1G01010)

Pipeline stages (each cached to disk, skipped on re-run):
  1. filter.csv          — remove low-variance genes/samples
                           NOTE: NO log transform applied here.
                           RMA output is already in log2 space (typically range 2–14).
                           Applying log2 again would be a double-transform error.
  2. filter_norm.csv     — mean-centered version of filter.csv for visualization/reference
                           (mirrors RNA-seq filter_norm.csv for comparability)
  3. combat.csv          — ComBat batch correction (Gaussian model, log2-space)
                           Uses sva::ComBat, NOT ComBat_seq.
                           ComBat expects log2-transformed continuous data — correct for RMA.
                           ComBat output remains in log2 space.
  4. rankin.csv          — Rank-in normalization (cross-platform integration)
                           No log transform before Rank-in: input is already log2 from ComBat.
                           This mirrors the RNA-seq pipeline where log1p is applied to
                           raw ComBat-seq counts before Rank-in — the principle is the same
                           (Rank-in receives log-space data in both pipelines).

ComBat reference:
  Johnson et al. (2007) Biostatistics 8(1):118-127
  https://doi.org/10.1093/biostatistics/kxj037

Rank-in reference:
  Tang et al. (2021) Briefings in Bioinformatics 22(4)
"""

import os
import sys
from collections import Counter
import json
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
import warnings
from scipy.sparse.linalg import svds
from sklearn.decomposition import TruncatedSVD

# --- R integration ---
import rpy2.robjects as ro
from rpy2.robjects import numpy2ri, pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr

module_dir = "./"
sys.path.append(module_dir)
from src.constants import STORAGE_DIR,SAMPLE_STUDY_MAP, LABELS_PATH,RNA_USED # noqa: E402
from src.data_analisys.utils.cluster_exploration_utils_final import get_gsm_id # noqa: E402

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
MICROARRAY_COMBINED = os.path.join(STORAGE_DIR, "final_data", "RMA_Microarray_Combined.csv")
MICROARRAY_DATA_DIR = os.path.join(STORAGE_DIR, "final_data")
MICROARRAY_FIGURES  = os.path.join(STORAGE_DIR, "figures", "microarray_filtering")


# ---------------------------------------------------------------------------
# Helpers — batch label extraction
# ---------------------------------------------------------------------------

def get_gse_from_col(col: str) -> str:
    """
    Looks up the GSE study ID for a GSM sample ID using SAMPLE_STUDY_MAP.
    Falls back to splitting on '_' for any non-GSM formatted columns.
    """
    try:
        sample_key = col.split(".", maxsplit=1)[0]  # strip any suffix
        return SAMPLE_STUDY_MAP.loc[sample_key, "StudyID"]
    except KeyError:
        # Fallback for columns already in GSE_GSM format
        return col.split("_", maxsplit=1)[0]


def get_batch_labels(columns) -> list[str]:
    """Returns one batch label (GSE ID) per column."""
    return [get_gse_from_col(c) for c in columns]

# ---------------------------------------------------------------------------
# 1. Filtering + summary plot
# ---------------------------------------------------------------------------

def plot_filtering_summary(df_before: pd.DataFrame, df_after: pd.DataFrame, output_path: str):
    print("  [Plot] Generating filtering summary...")
    fig = plt.figure(figsize=(16, 12))
    fig.suptitle("Microarray Data Filtering Summary", fontsize=20, fontweight="bold")

    gs = fig.add_gridspec(2, 2, height_ratios=[1, 2.5])

    ax_bars = fig.add_subplot(gs[0, :])
    categories = ["Genes (Rows)", "Samples (Columns)"]
    before_counts = [df_before.shape[0], df_before.shape[1]]
    after_counts  = [df_after.shape[0],  df_after.shape[1]]
    x, w = np.arange(2), 0.35
    ax_bars.bar(x - w / 2, before_counts, w, label="Before", color="#ff9999", edgecolor="black")
    ax_bars.bar(x + w / 2, after_counts,  w, label="After",  color="#66b3ff", edgecolor="black")
    ax_bars.set_xticks(x)
    ax_bars.set_xticklabels(categories, fontsize=12)
    ax_bars.set_ylabel("Count")
    ax_bars.legend()
    ax_bars.set_title("Genes and Samples Before vs After Filtering")
    for i, v in enumerate(before_counts):
        ax_bars.text(i - w / 2, v * 1.01, f"{v:,}", ha="center", fontweight="bold")
    for i, v in enumerate(after_counts):
        ax_bars.text(i + w / 2, v * 1.01, f"{v:,}", ha="center", fontweight="bold")

    for ax, df, label in [
        (fig.add_subplot(gs[1, 0]), df_before, "Before"),
        (fig.add_subplot(gs[1, 1]), df_after,  "After"),
    ]:
        ax.imshow(df.isna(), aspect="auto", cmap="viridis", interpolation="nearest")
        pct = (df.isna().sum().sum() / df.size) * 100
        ax.set_title(f"Missingness {label}\n({pct:.1f}% missing)")
        ax.set_xlabel("Samples")
        ax.set_ylabel("Genes")

    plt.tight_layout(rect=[0, 0, 1, 0.96])  # pyright: ignore[reportArgumentType]
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  [Plot] Saved to {output_path}")


def run_filtering(
    raw_df: pd.DataFrame,
    gene_nan_pct: float = 100.0,
    sample_nan_pct: float = 60.0,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Microarray filtering — mirrors RNA-seq run_filtering() exactly.

    RMA input is already in log2 space with no zeros by design, so the
    zero-or-NaN mask is equivalent to a pure NaN mask here. We keep the
    same (zero | NaN) logic as RNA-seq for consistency.

    Returns
    -------
    filtered_df : log2-space filtered matrix  → saved as filter.csv
    norm_df     : mean-centered version        → saved as filter_norm.csv
                  (analogous to RNA-seq filter_norm.csv; useful for QC plots)
    """
    print(f"  [Filter] Input: {raw_df.shape[0]} genes × {raw_df.shape[1]} samples")
    print(f"  [Filter] Input value range: [{raw_df.values.min():.2f}, {raw_df.values.max():.2f}]")
    print("  [Filter] NOTE: RMA output is already log2-transformed. No additional log applied.")

    # Validate log2 scale — RMA data is typically in range 2–14
    # If max > 20 the input is likely raw intensities, not RMA output → abort early
    if raw_df.values.max() > 20:
        raise ValueError(
            f"Input max value is {raw_df.values.max():.1f}, which suggests raw intensities "
            "rather than RMA log2 output (expected range ~2–14). "
            "Please provide RMA-processed data."
        )

    # STEP 1: Filter Genes (Rows)
    is_zero_or_nan = raw_df.isna() | (raw_df == 0)
    invalid_pct_per_gene = is_zero_or_nan.mean(axis=1) * 100
    raw_df = raw_df.loc[invalid_pct_per_gene <= gene_nan_pct, :]
    print(f"  [Filter] After gene filtering (≤{gene_nan_pct}% 0/NaN): {raw_df.shape[0]} genes × {raw_df.shape[1]} samples")

    # STEP 2: Filter Samples (Columns)
    is_zero_or_nan = raw_df.isna() | (raw_df == 0)
    invalid_pct_per_sample = is_zero_or_nan.mean(axis=0) * 100
    raw_df = raw_df.loc[:, invalid_pct_per_sample <= sample_nan_pct]
    print(f"  [Filter] After sample filtering (≤{sample_nan_pct}% 0/NaN): {raw_df.shape[0]} genes × {raw_df.shape[1]} samples")

    # filter.csv  = log2-space RMA values, filtered only (no further transform)
    filtered_df = pd.DataFrame(raw_df.values, index=raw_df.index, columns=raw_df.columns)

    return filtered_df, filtered_df # due to RMA normal filter matrix is already log normalized


# ---------------------------------------------------------------------------
# 2. ComBat  (Gaussian batch correction, log2-space)
# ---------------------------------------------------------------------------

def run_combat(
    log2_df: pd.DataFrame,
    batch_labels: list[str],
    covar_df: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """
    Runs R's sva::ComBat (Gaussian model) on log2-transformed microarray data.

    Unlike ComBat-seq (negative-binomial, for RNA-seq counts), ComBat uses a
    Gaussian/empirical-Bayes model and requires log2-transformed continuous input.
    RMA output satisfies this requirement directly.

    ComBat output is in the same log2 space as the input — no post-hoc
    log transform is needed before Rank-in (contrast with RNA-seq where
    log1p is applied to raw ComBat-seq counts before Rank-in).

    Parameters
    ----------
    log2_df      : genes × samples log2-space expression matrix (no NaN allowed)
    batch_labels : one batch string per column (GSE study IDs)
    covar_df     : optional DataFrame of biological covariates (samples as rows)
    """
    print("  [ComBat] Loading R sva package...")
    try:
        sva   = importr("sva")
        stats = importr("stats")
        base  = importr("base")
    except Exception as e:
        msg = "R package 'sva' not found. Install with:\n  BiocManager::install('sva')"
        raise ImportError(msg) from e

    df = log2_df.fillna(0)
    if not np.issubdtype(df.values.dtype, np.floating):
        df = df.astype(float)

    print(f"  [ComBat] Input: {df.shape[0]} genes × {df.shape[1]} samples")

    # Transfer data to R using localconverter (rpy2 3.5 compatible)
    with localconverter(ro.default_converter + numpy2ri.converter):
        dat_r = ro.conversion.py2rpy(df.values)

    batch_r = ro.StrVector(batch_labels)

    mod_r = ro.r("NULL")
    if covar_df is not None and not covar_df.empty:
        print(f"  [ComBat] Building covariate model for: {list(covar_df.columns)}")
        covar_df = covar_df.astype(str)
        with localconverter(ro.default_converter + pandas2ri.converter):
            covar_r = ro.conversion.py2rpy(covar_df)
        formula  = ro.Formula("~ " + " + ".join(covar_df.columns))
        mod_r    = stats.model_matrix(formula, data=covar_r)

    print("  [ComBat] Running sva::ComBat (Gaussian model)...")
    combat_r = sva.ComBat(dat=dat_r, batch=batch_r, mod=mod_r, par_prior=True, prior_plots=False)

    del dat_r, batch_r, mod_r
    base.gc()

    combat_np = np.array(combat_r)
    del combat_r
    base.gc()

    result = pd.DataFrame(combat_np, index=log2_df.index, columns=log2_df.columns)
    print("  [ComBat] Done.")
    return result


# ---------------------------------------------------------------------------
# 3. Rank-in normalization
# ---------------------------------------------------------------------------
def run_rank_in_normalization(
    df: pd.DataFrame,
    sample_classes: pd.Series | None,
    n_bins: int = 100,
    k: int | None = None,
    k_max: int = 10,
    out_path: str | None = None,
) -> pd.DataFrame:
    """
    Rank-In normalization per Tang et al. (2021).
    Position-indexed and immune to duplicate column names.
    """
    # =============================================================================
    #   THOROUGH PIPELINE DIAGNOSTIC ENGINE
    # =============================================================================
    print("\n" + "═"*80)
    print(" [RANK-IN DIAGNOSIS] EXAMINING INPUT EXPRESSION MATRIX AND LABELS")
    print("═"*80)
    
    print("─── 1. Core Matrix Geometry ───")
    print(f"  • Matrix Type:              {type(df)}")
    print(f"  • Shape (Genes × Samples):  {df.shape}")
    print(f"  • Row (Gene) Index Type:    {type(df.index)} (Name: '{df.index.name}')")
    print(f"  • Column (Sample) ID Type:  {type(df.columns)} (Name: '{df.columns.name}')")
    print(f"  • Is Row Index MultiIndex?  {isinstance(df.index, pd.MultiIndex)}")
    print(f"  • Is Column MultiIndex?     {isinstance(df.columns, pd.MultiIndex)}")
    
    print("\n─── 2. Class Labels Metadata (sample_classes) ───")
    print(f"  • sample_classes Type:      {type(sample_classes)}")
    if sample_classes is not None:
        print(f"  • Label Dimensions/Length:  {sample_classes.shape}")
    else:
        print("  • WARNING: sample_classes is explicitly 'None'!")
        print("    ↳ Step 3 Centering calculations will use global/fallback column vectors.")
        
    print("\n─── 3. Data Integrity & Column Formats ───")
    print(f"  • Total Null/NaN values inside matrix: {df.isna().sum().sum()}")
    
    print("\n─── 4. Sparsity & Distribution Profiling ───")
    total_cells = df.size
    if total_cells > 0:
        exact_zeros = (df == 0).sum().sum()
        sparsity_pct = (exact_zeros / total_cells) * 100
        print(f"  • Global Sparsity: {exact_zeros} true zeros out of {total_cells} cells ({sparsity_pct:.2f}%)")
    print("═"*80 + "\n")

    df_og_index = df.index.copy()
    df_og_col = df.columns.copy()

    # ------------------------------------------------------------------ #
    # Mechanism to get the sample classes from LABELS_PATH if not provided#
    # ------------------------------------------------------------------ #
    if sample_classes is None:
        print(f"[Rank-In] sample_classes not provided. Loading from {LABELS_PATH}...")
        
        groupings = ['treatment', 'tissue']
        raw_metadata = {}

        for path in sorted(Path(LABELS_PATH).glob("*.json")):
            if path.name.startswith("map_"):
                continue
            try:
                data = json.loads(path.read_text())
            except Exception as e:
                warnings.warn(f"Failed to read metadata file {path.name}: {e}")
                continue
            
            if isinstance(data, dict):
                for sample_id, sample_dict in data.items():
                    if isinstance(sample_dict, dict):
                        raw_metadata[str(sample_id).upper()] = sample_dict
            elif isinstance(data, list):
                for item in data:
                    if isinstance(item, dict):
                        for id_key in ['sample_id', 'id', 'name', 'run', 'sample']:
                            if id_key in item:
                                raw_metadata[str(item[id_key]).upper()] = item
                                break

        if RNA_USED:
            print("[Rank-In] RNA mode: remapping SRR IDs to GSM IDs in the df...")
            df.columns = df.columns.map(lambda x: get_gsm_id(x.split('_')[1]))
            print(f"[Rank-In] RNA mode: remapped form {df_og_col} to {df.columns} temporarily...")

        class_mapping = {}
        for sample_id, sample in raw_metadata.items():
            vals_parts = []
            for label in groupings:
                raw = sample.get(label, [])
                if not isinstance(raw, list):
                    raw = [raw]
                vals = [str(item.get("val") or item.get("value") or "Unspecified") if isinstance(item, dict) else str(item) if item else "Unspecified" for item in raw]
                vals_parts.append("-".join(sorted(vals)))
            class_mapping[sample_id] = "_".join(vals_parts)

        # FIX 1: Use a list to prevent dictionary collisions on duplicate column IDs
        classes_list = []
        for col in df.columns:
            col_str = str(col).upper()
            if col_str in class_mapping:
                classes_list.append(class_mapping[col_str])
            else:
                matched = False
                for map_key, map_val in class_mapping.items():
                    if map_key in col_str or col_str in map_key:
                        classes_list.append(map_val)
                        matched = True
                        break
                if not matched:
                    classes_list.append(np.nan)

        # Construct Series using the original explicit list order
        sample_classes = pd.Series(classes_list, index=df.columns)
        print(f"  -> Successfully generated classes for {sample_classes.notna().sum()} / {len(df.columns)} samples.")
        print(f"len of sample_classes: {len(sample_classes)}")

    # ------------------------------------------------------------------ #
    # Validate inputs                                                      #
    # ------------------------------------------------------------------ #
    if sample_classes.isna().any():
        missing = sample_classes[sample_classes.isna()].index.tolist()
        raise ValueError(f"sample_classes has missing labels for {len(missing)} samples. Examples: {missing[:5]}")
 
    classes = sample_classes.unique()
    if len(classes) < 2:
        raise ValueError("At least two distinct class labels are required.")
 
    # ------------------------------------------------------------------ #
    # Step 1: Binned rank transformation                                   #
    # ------------------------------------------------------------------ #
    print(f"\n[Rank-In] Step 1: Binned rank transformation ({n_bins} bins)...")
    rank_matrix = df.rank(pct=True, method="average")
    binned_matrix = pd.DataFrame(
        np.ceil(rank_matrix.values * n_bins).astype(float),
        index=df.index,
        columns=df.columns,
    )
 
    # ------------------------------------------------------------------ #
    # Step 2: Weight ranks by expression intensity slope                   #
    # ------------------------------------------------------------------ #
    print("[Rank-In] Step 2: Weighting ranks by expression intensity slope...")
    weighted_matrix = np.zeros_like(binned_matrix.values, dtype=float)
 
    # FIX 2: Iterate and slice by numerical column index location (.iloc)
    for j in range(binned_matrix.shape[1]):
        r = binned_matrix.iloc[:, j].values.astype(float)   # Guaranteed 1D array
        e = df.iloc[:, j].values.astype(float)              # Guaranteed 1D array
 
        try:
            coeffs = np.polyfit(r, e, deg=2)
            a, b = coeffs[0], coeffs[1]
        except (np.linalg.LinAlgError, ValueError) as exc:
            warnings.warn(f"Quadratic fit failed for column index {j}. Using identity weighting.", RuntimeWarning)
            a, b = 0.0, 1.0
 
        w = 2 * a * r + b
        weighted_matrix[:, j] = r * w

    # weighted_df = pd.DataFrame(weighted_matrix, index=df_og_index, columns=df_og_col)
 
    # ------------------------------------------------------------------ #
    # Step 3: SVD — estimate and subtract nonbiological effects            #
    # ------------------------------------------------------------------ #
    print("[Rank-In] Step 3: SVD — estimating and subtracting nonbiological effects...")
 
    # FIX 3: Compute class means safely using positional indexes instead of label matching
    group_mean_matrix = np.zeros_like(weighted_matrix, dtype=float)
    for cls in classes:
        class_indices = np.where(sample_classes.values == cls)[0]
        if len(class_indices) > 0:
            class_mean = np.mean(weighted_matrix[:, class_indices], axis=1)
            for idx in class_indices:
                group_mean_matrix[:, idx] = class_mean
 
    variance_matrix = weighted_matrix - group_mean_matrix
 
    n_genes, n_samples = variance_matrix.shape
    k_max_safe = min(k_max, n_samples - 1, n_genes - 1)
 
    if k is not None:
        k_use = min(int(k), k_max_safe)
        print(f"  -> Using user-supplied k={k_use} SVD components.")
    else:
        _, singular_values, _ = svds(variance_matrix, k=k_max_safe)
        singular_values = singular_values[::-1]
        ratios = singular_values[:-1] / (singular_values[1:] + 1e-12)
        k_use = int(np.argmax(ratios)) + 1
        print(f"  -> Auto-selected k={k_use} dominant SVD components.")
 
    U, s, Vt = np.linalg.svd(variance_matrix, full_matrices=False)
 
    U_k  = U[:, :k_use]
    s_k  = np.diag(s[:k_use])
    Vt_k = Vt[:k_use, :]
    nonbio_effects = U_k @ s_k @ Vt_k
 
    adjusted_values = weighted_matrix - nonbio_effects
    adjusted_df = pd.DataFrame(adjusted_values, index=df_og_index, columns=df_og_col)
 
    print("[Rank-In] Normalization complete.")
    if out_path is not None:
        adjusted_df.to_csv(out_path)
        print(f"  -> Saved adjusted matrix to: {out_path}")
 
    return adjusted_df

def run_rank_in_normalization_old(
    df: pd.DataFrame,
    n_bins: int = 100,
    variance_threshold: float = 0.95,
    out_path: str | None = None,
) -> pd.DataFrame:
    """
    Rank-in normalization for cross-platform integration (Tang et al., 2021).

    Identical implementation to the RNA-seq pipeline — the algorithm is
    platform-agnostic and operates on any log-space expression matrix.
    """
    print(f"\n[Rank-in] Step 1: Rank transformation ({n_bins} bins)...")
    rank_matrix  = df.rank(pct=True, method="average")
    binned_matrix = pd.DataFrame(
        np.ceil(rank_matrix.values * n_bins),
        index=df.index,
        columns=df.columns,
    )

    print("[Rank-in] Step 2: Centering for SVD...")
    gene_means     = binned_matrix.mean(axis=1)
    binned_centered = binned_matrix.sub(gene_means, axis=0)
    X = binned_centered.T  # sklearn SVD expects (n_samples, n_features)

    print("[Rank-in] Step 3: SVD adjustment...")
    max_components = min(X.shape) - 1
    svd_test = TruncatedSVD(n_components=max_components, random_state=42)
    svd_test.fit(X)

    cumulative_variance = np.cumsum(svd_test.explained_variance_ratio_)
    optimal_k = np.argmax(cumulative_variance >= variance_threshold) + 1
    print(f"  -> Selected top {optimal_k} SVD components ({variance_threshold * 100:.0f}% variance retained).")

    svd_final    = TruncatedSVD(n_components=optimal_k, random_state=42)
    X_reduced    = svd_final.fit_transform(X)
    X_reconstructed = svd_final.inverse_transform(X_reduced)

    adjusted_df = pd.DataFrame(
        X_reconstructed.T,
        index=binned_matrix.index,
        columns=binned_matrix.columns,
    ).add(gene_means, axis=0)

    print("[Rank-in] Normalization complete.")
    save_to = out_path if out_path is not None else os.path.join(MICROARRAY_DATA_DIR, "rankin.csv")
    adjusted_df.to_csv(save_to)
    return adjusted_df


# ---------------------------------------------------------------------------
# 4. Main pipeline
# ---------------------------------------------------------------------------

def run_microarray_preprocessing():
    os.makedirs(MICROARRAY_DATA_DIR, exist_ok=True)
    os.makedirs(MICROARRAY_FIGURES,  exist_ok=True)

    filter_path      = os.path.join(MICROARRAY_DATA_DIR, "filter.csv")
    filter_norm_path = os.path.join(MICROARRAY_DATA_DIR, "filter_norm.csv")
    combat_path      = os.path.join(MICROARRAY_DATA_DIR, "combat_seq_norm.csv")
    rankin_path      = os.path.join(MICROARRAY_DATA_DIR, "rankin.csv")

    # ── Stage 1: Filter ──────────────────────────────────────────────────────
    if os.path.exists(filter_path):
        print("Loading cached filter.csv...")
        filtered_df = pd.read_csv(filter_path, index_col=0)
        norm_df = pd.read_csv(filter_norm_path, index_col=0)
    else:
        print(f"Loading RMA matrix from {MICROARRAY_COMBINED}...")
        raw_df = pd.read_csv(MICROARRAY_COMBINED, index_col=0)
        print(f"  Raw matrix: {raw_df.shape[0]} genes × {raw_df.shape[1]} samples")

        filtered_df, norm_df = run_filtering(raw_df)

        plot_filtering_summary(
            raw_df,
            filtered_df,
            os.path.join(MICROARRAY_FIGURES, "filtering_summary.svg"),
        )
        filtered_df.to_csv(filter_path)
        norm_df.to_csv(filter_norm_path)
        print(f"Saved filtered matrix  → {filter_path}")
        print(f"Saved norm matrix      → {filter_norm_path}")

    # ── Stage 2: ComBat ───────────────────────────────────────────────────────
    if os.path.exists(combat_path):
        print("Loading cached combat.csv...")
        combat_df = pd.read_csv(combat_path, index_col=0)
    else:
        print("\nRunning ComBat batch correction...")

        batch_labels  = get_batch_labels(filtered_df.columns)
        study_counts  = Counter(batch_labels)

        # ComBat requires ≥ 2 samples per batch
        single_batches = {s for s, c in study_counts.items() if c < 2}
        if single_batches:
            print(f"  Removing {len(single_batches)} single-sample batches: {single_batches}")

        valid_cols    = [c for c, b in zip(filtered_df.columns, batch_labels, strict=False) if b not in single_batches]
        valid_batches = [b for b in batch_labels if b not in single_batches]
        df_for_combat = filtered_df[valid_cols]

        print(f"  NaN in input: {df_for_combat.isna().sum().sum()} → filling with 0 for ComBat")
        df_for_combat = df_for_combat.fillna(0)

        combat_df = run_combat(df_for_combat, valid_batches)
        combat_df.to_csv(combat_path)
        print(f"Saved ComBat result    → {combat_path}")

    # ── Stage 3: Rank-in ──────────────────────────────────────────────────────
    if os.path.exists(rankin_path):
        print("Loading cached rankin.csv...")
        rankin_df = pd.read_csv(rankin_path, index_col=0)
    else:
        print("\nRunning Rank-in normalization on normalized filter output...")
        # ComBat output is already in log2 space (same space as its input).
        # No additional log transform is applied here — contrast with RNA-seq
        # where log1p(combat_seq_counts) is needed because ComBat-seq outputs
        # raw counts. Here the data is already log-transformed.
        rankin_df = run_rank_in_normalization(
            df=norm_df,
            sample_classes=None,
            out_path=rankin_path,
        )
        print(f"Saved Rank-in result   → {rankin_path}")

    print("\n=== Pipeline Complete ===")
    print(f"  filter.csv      : {filtered_df.shape}  (log2 RMA, filtered)")
    print(f"  filter_norm.csv : {filtered_df.shape}  (mean-centered, for QC)")
    print(f"  combat.csv      : {combat_df.shape}   (log2, batch corrected)")
    print(f"  rankin.csv      : {rankin_df.shape}   (rank normalized)")
    return filtered_df, combat_df, rankin_df


if __name__ == "__main__":
    run_microarray_preprocessing()