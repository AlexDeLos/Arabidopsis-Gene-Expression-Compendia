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

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
import warnings
# --- R integration ---
import rpy2.robjects as ro
from rpy2.robjects import numpy2ri, pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr
from sklearn.decomposition import TruncatedSVD

module_dir = "./"
sys.path.append(module_dir)
from src.constants import STORAGE_DIR,SAMPLE_STUDY_MAP, LABELS_PATH  # noqa: E402

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
 
    Parameters
    ----------
    df : pd.DataFrame
        Log2-transformed expression matrix, shape (n_genes, n_samples).
        Genes as rows, samples as columns.
    sample_classes : pd.Series
        Class label for each sample (index must match df.columns).
        E.g. pd.Series({"S1": "cancer", "S2": "normal", ...})
        Used to compute per-group means for the centering step (Step 3),
        matching the paper's use of experimental vs. control group means.
    n_bins : int
        Number of rank bins (default 100, as in the paper).
    k : int or None
        Number of SVD components to treat as nonbiological effects.
        If None (default), k is selected automatically by finding the
        elbow / dominant drop in singular values (small k, typically 1–5).
        The paper says k should capture the dominant nonbiological variance,
        NOT explain a large fraction of total variance.
    k_max : int
        Upper bound when auto-selecting k (default 10).
    out_path : str or None
        If provided, the adjusted matrix is saved to this CSV path.
 
    Returns
    -------
    pd.DataFrame
        Adjusted ranking matrix (genes × samples), nonbiological effects removed.
    """
# ------------------------------------------------------------------ #
    # Mechanism to get the sample classes from LABELS_PATH if not provided#
    # ------------------------------------------------------------------ #
    if sample_classes is None:
        print(f"[Rank-In] sample_classes not provided. Loading from {LABELS_PATH}...")
        
        # Define groupings to build the class signature (matching your script)
        groupings = ['treatment', 'tissue']
        raw_metadata = {}

        # 1. Parse all files and map Sample IDs to their raw dictionary entries
        for path in sorted(Path(LABELS_PATH).glob("*.json")):
            if path.name.startswith("map_"):
                continue
            try:
                data = json.loads(path.read_text())
            except Exception as e:
                warnings.warn(f"Failed to read metadata file {path.name}: {e}")
                continue
            
            if isinstance(data, dict):
                # Assuming top-level keys are sample identifiers (e.g., GSM/SRR IDs)
                for sample_id, sample_dict in data.items():
                    if isinstance(sample_dict, dict):
                        raw_metadata[str(sample_id)] = sample_dict
            elif isinstance(data, list):
                # If it's a list, look inside the sample dictionary for common ID keys
                for item in data:
                    if isinstance(item, dict):
                        for id_key in ['sample_id', 'id', 'name', 'run', 'sample']:
                            if id_key in item:
                                raw_metadata[str(item[id_key])] = item
                                break

        # 2. Extract standardized combinations to formulate class strings
        class_mapping = {}
        for sample_id, sample in raw_metadata.items():
            vals_parts = []
            for label in groupings:
                raw = sample.get(label, [])
                if not isinstance(raw, list):
                    raw = [raw]
                
                vals = []
                for item in raw:
                    if isinstance(item, dict):
                        vals.append(item.get("val") or item.get("value") or "Unspecified")
                    else:
                        vals.append(str(item) if item else "Unspecified")
                
                # Combine multiple values within a single factor (e.g. ['leaf', 'stem'] -> 'leaf-stem')
                vals_parts.append("-".join(sorted(vals)))
            
            # Combine across multiple groupings (e.g., 'drought' + 'leaf' -> 'drought_leaf')
            class_mapping[sample_id] = "_".join(vals_parts)

        # 3. Align class labels strictly with the columns present in df
        series_dict = {}
        for col in df.columns:
            col_str = str(col)
            if col_str in class_mapping:
                series_dict[col] = class_mapping[col_str]
            else:
                # Soft fallback matching: checks if the dict key is embedded inside the column name
                matched = False
                for k, v in class_mapping.items():
                    if k in col_str or col_str in k:
                        series_dict[col] = v
                        matched = True
                        break
                if not matched:
                    series_dict[col] = np.nan

        sample_classes = pd.Series(series_dict)
        print(f"  -> Successfully generated classes for {sample_classes.notna().sum()} / {len(df.columns)} samples.")    # ------------------------------------------------------------------ #
    # Validate inputs                                                      #
    # ------------------------------------------------------------------ #
    if not df.columns.equals(sample_classes.index):
        # Align just in case columns and series index are in different order
        sample_classes = sample_classes.reindex(df.columns)
    if sample_classes.isna().any():
        raise ValueError(
            "sample_classes has missing labels for some samples. "
            "Every sample column must have a class label."
        )
 
    classes = sample_classes.unique()
    if len(classes) < 2:
        raise ValueError(
            "At least two distinct class labels are required "
            "(e.g. 'cancer' and 'normal')."
        )
 
    # ------------------------------------------------------------------ #
    # Step 1: Binned rank transformation                                   #
    # ------------------------------------------------------------------ #
    print(f"\n[Rank-In] Step 1: Binned rank transformation ({n_bins} bins)...")
 
    # Rank genes within each sample, low → high.
    # pct=True gives fractional rank in (0, 1]; ceil * n_bins → [1, n_bins].
    rank_matrix = df.rank(pct=True, method="average")
    binned_matrix = pd.DataFrame(
        np.ceil(rank_matrix.values * n_bins).astype(float),
        index=df.index,
        columns=df.columns,
    )  # values in [1, n_bins]
 
    # ------------------------------------------------------------------ #
    # Step 2: Weight ranks by expression intensity slope                   #
    # ------------------------------------------------------------------ #
    print("[Rank-In] Step 2: Weighting ranks by expression intensity slope...")
 
    weighted_matrix = np.zeros_like(binned_matrix.values, dtype=float)
 
    for j, col in enumerate(binned_matrix.columns):
        r = binned_matrix[col].values.astype(float)   # bin ranks [1, n_bins]
        e = df[col].values.astype(float)               # log2 expression
 
        # Fit quadratic:  e = a*r^2 + b*r + c  (least-squares, per paper)
        try:
            coeffs = np.polyfit(r, e, deg=2)           # [a, b, c]
            a, b = coeffs[0], coeffs[1]
        except (np.linalg.LinAlgError, ValueError) as exc:
            warnings.warn(
                f"Quadratic fit failed for sample '{col}' ({exc}). "
                "Falling back to identity weighting (w=1) for this sample.",
                RuntimeWarning,
                stacklevel=2,
            )
            a, b = 0.0, 1.0                             # w = 1 → R_ij = r_ij
 
        # Weight:  w_ij = 2a*r_ij + b  (derivative of the quadratic)
        w = 2 * a * r + b
        weighted_matrix[:, j] = r * w                  # R_ij = r_ij * w_ij
 
    weighted_df = pd.DataFrame(
        weighted_matrix,
        index=binned_matrix.index,
        columns=binned_matrix.columns,
    )
 
    # ------------------------------------------------------------------ #
    # Step 3: SVD — estimate and subtract nonbiological effects            #
    # ------------------------------------------------------------------ #
    print("[Rank-In] Step 3: SVD — estimating and subtracting nonbiological effects...")
 
    # --- FIX 1: center using per-class group means, not a grand mean ----
    #
    # The paper models:  R = x + y + α_j + ε
    # and approximates the "real" signal as:  R_real ≈ Me_ij
    # where Me_ij is the mean of gene i within each experimental group.
    # Subtracting these group means leaves the variance matrix R̃ that
    # is dominated by nonbiological (batch/platform) effects.
    #
    group_mean_matrix = pd.DataFrame(
        np.zeros_like(weighted_df.values, dtype=float),
        index=weighted_df.index,
        columns=weighted_df.columns,
    )
    for cls in classes:
        cols_in_class = sample_classes[sample_classes == cls].index.tolist()
        class_mean = weighted_df[cols_in_class].mean(axis=1).values  # (n_genes,)
        for col in cols_in_class:
            group_mean_matrix[col] = class_mean
 
    # Variance matrix:  R̃ = R - Me_ij  (genes × samples)
    variance_matrix = weighted_df - group_mean_matrix   # shape: (n_genes, n_samples)
 
    # --- FIX 2: choose k as the number of DOMINANT singular values ------
    #
    # The paper says k captures the main nonbiological variables.
    # This should be a small number (typically 1–5), NOT chosen to explain
    # a large fraction of variance (which would remove biological signal).
    #
    # Auto-selection: compute up to k_max singular values and find the
    # elbow — the point where the ratio between consecutive singular
    # values drops most sharply, indicating the end of the dominant
    # nonbiological structure.
    #
    n_genes, n_samples = variance_matrix.shape
    k_max_safe = min(k_max, n_samples - 1, n_genes - 1)
 
    if k is not None:
        # User-supplied k: trust it but warn if it looks large
        k_use = int(k)
        if k_use > k_max_safe:
            warnings.warn(
                f"Requested k={k_use} exceeds k_max_safe={k_max_safe}. "
                f"Clamping to {k_max_safe}.",
                RuntimeWarning,
                stacklevel=2,
            )
            k_use = k_max_safe
        print(f"  -> Using user-supplied k={k_use} SVD components.")
    else:
        # Compute top-(k_max_safe) singular values via sparse SVD
        # svds returns singular values in *ascending* order → reverse
        _, singular_values, _ = svds(
            variance_matrix.values.astype(float), k=k_max_safe
        )
        singular_values = singular_values[::-1]  # descending
 
        # Elbow detection: largest *ratio drop* between consecutive values
        # s[i] / s[i+1] — biggest ratio means sharpest drop after component i
        ratios = singular_values[:-1] / (singular_values[1:] + 1e-12)
        k_use = int(np.argmax(ratios)) + 1   # +1 because k counts from 1
 
        print(
            f"  -> Auto-selected k={k_use} dominant SVD components "
            f"(singular values: {np.round(singular_values[:k_use+2], 2)})."
        )
 
    # --- Full SVD for reconstruction (economy / thin SVD) ---------------
    #
    # We need U, Σ, Vᵀ of the variance matrix (genes × samples).
    # Using numpy's full SVD with full_matrices=False for efficiency.
    #
    U, s, Vt = np.linalg.svd(variance_matrix.values.astype(float), full_matrices=False)
 
    # Nonbiological effect matrix:  y ≈ U_k * Σ_k * V_kᵀ
    U_k  = U[:, :k_use]                        # (n_genes,  k)
    s_k  = np.diag(s[:k_use])                  # (k, k)
    Vt_k = Vt[:k_use, :]                       # (k, n_samples)
    nonbio_effects = U_k @ s_k @ Vt_k          # (n_genes, n_samples)
 
    # Adjusted matrix:  R_adj = R - y  (subtract nonbio from *weighted* matrix,
    # not the variance matrix — we want to keep the biological signal Me_ij)
    adjusted_values = weighted_df.values - nonbio_effects
 
    adjusted_df = pd.DataFrame(
        adjusted_values,
        index=weighted_df.index,
        columns=weighted_df.columns,
    )
 
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