"""
rnaseq_norm_and_analysis.py

Normalization pipeline for the combined RNA-seq count matrix.

Column format: {GSE_id}_{SRR_id}   (e.g. GSE201685_SRR18937094)
Row format   : gene_id              (e.g. AT1G01010)

Pipeline stages (each cached to disk, skipped on re-run):
  1. filter.csv          — remove low-coverage genes/samples, zero → NaN
  2. combat.csv      — ComBat-seq batch correction (count-space, negative-binomial)
  3. rankin.csv          — Rank-in normalization (cross-platform integration)

ComBat-seq reference:
  Zhang et al. (2020) Genome Biology 21:257
  https://pmc.ncbi.nlm.nih.gov/articles/PMC7518324/
"""

import os
import sys
from collections import Counter

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

module_dir = "./"
sys.path.append(module_dir)
from src.constants import STORAGE_DIR  # provides STORAGE_DIR, COMBINED_DATA_OUTPUT_FILE, etc.  # noqa: E402
from src.data_importing.data_norm_and_analisys import run_rank_in_normalization,run_combat  # noqa: E402

# ---------------------------------------------------------------------------
# Paths  (edit RNASEQ_DATA_DIR / RNASEQ_FIGURES_DIR in src/constants.py
#         or override here)
# ---------------------------------------------------------------------------
RNASEQ_COMBINED = os.path.join(STORAGE_DIR, "final_data/rnaseq_processed", "Salmon_RNAseq_Combined_TPM.csv")
RNASEQ_DATA_DIR = os.path.join(STORAGE_DIR, "final_data", "rnaseq_processed")
RNASEQ_FIGURES = os.path.join(STORAGE_DIR, "figures", "rnaseq_filtering")


# ---------------------------------------------------------------------------
# Helpers — batch label extraction
# ---------------------------------------------------------------------------


def get_gse_from_col(col: str) -> str:
    """
    Extracts the GSE study ID from a column named {GSE_id}_{SRR_id}.
    Works whether the column has already been renamed to a GSM or still
    carries the original GSE_SRR format.
    """
    return col.split("_", maxsplit=1)[0]


def get_batch_labels(columns) -> list[str]:
    """Returns one batch label (GSE ID) per column."""
    return [get_gse_from_col(c) for c in columns]


# ---------------------------------------------------------------------------
# 1. Filtering
# ---------------------------------------------------------------------------


def plot_filtering_summary(df_before: pd.DataFrame, df_after: pd.DataFrame, output_path: str):
    print("  [Plot] Generating filtering summary...")
    fig = plt.figure(figsize=(16, 12))
    fig.suptitle("RNA-seq Data Filtering Summary", fontsize=20, fontweight="bold")

    gs = fig.add_gridspec(2, 2, height_ratios=[1, 2.5])

    ax_bars = fig.add_subplot(gs[0, :])
    categories = ["Genes (Rows)", "Samples (Columns)"]
    before_counts = [df_before.shape[0], df_before.shape[1]]
    after_counts = [df_after.shape[0], df_after.shape[1]]
    x, w = np.arange(2), 0.35
    ax_bars.bar(x - w / 2, before_counts, w, label="Before", color="#ff9999", edgecolor="black")
    ax_bars.bar(x + w / 2, after_counts, w, label="After", color="#66b3ff", edgecolor="black")
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
        (fig.add_subplot(gs[1, 1]), df_after, "After"),
    ]:
        ax.imshow((df.isna() | (df == 0)), aspect="auto", cmap="viridis", interpolation="nearest")
        pct = ((df.isna() | (df == 0)).sum().sum() / df.size) * 100
        ax.set_title(f"Missingness {label}\n({pct:.1f}% missing)")
        ax.set_xlabel("Samples")
        ax.set_ylabel("Genes")

    plt.tight_layout(rect=[0, 0, 1, 0.96])  # pyright: ignore[reportArgumentType]
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  [Plot] Saved to {output_path}")


def run_filtering(raw_df: pd.DataFrame, gene_nan_pct: float = 100.0, sample_nan_pct: float = 60.0) -> tuple[pd.DataFrame,pd.DataFrame]:
    """
    RNA-seq sample and gene filtering:

    1. Removes genes (rows) if they have a value of NaN or 0 in more than
       `gene_nan_pct`% of the samples.
    2. Removes samples (columns) if they have a value of NaN or 0 in more than
       `sample_nan_pct`% of the remaining genes.

    Parameters
    ----------
    raw_df         : genes × samples count matrix (integer or fl
# ---------------------------------------------------------------------------
# 4. Main pipeline
# ---------------------------------------------------------------------------


oat counts)
    gene_nan_pct   : max % of (0 or NaN) per gene (row) to keep [default 20.0]
    sample_nan_pct : max % of (0 or NaN) per sample (col) to keep [default 20.0]
    """
    print(f"  [Filter] Input: {raw_df.shape[0]} genes × {raw_df.shape[1]} samples")

    # ==========================================
    # STEP 1: Filter Genes (Rows)
    # ==========================================
    is_zero_or_nan = raw_df.isna() | (raw_df == 0)

    # Calculate the percentage of 0/NaN values per gene (across columns -> axis=1)
    invalid_pct_per_gene = is_zero_or_nan.mean(axis=1) * 100

    # Keep only the genes (rows) that fall below or equal to the threshold
    raw_df = raw_df.loc[invalid_pct_per_gene <= gene_nan_pct, :]

    print(f"  [Filter] After gene filtering (≤{gene_nan_pct}% 0/NaN): {raw_df.shape[0]} genes × {raw_df.shape[1]} samples")

    # ==========================================
    # STEP 2: Filter Samples (Columns)
    # ==========================================
    # Recalculate 0/NaN mask for the newly reduced dataframe
    is_zero_or_nan = raw_df.isna() | (raw_df == 0)

    # Calculate the percentage of 0/NaN values per sample (across rows -> axis=0)
    invalid_pct_per_sample = is_zero_or_nan.mean(axis=0) * 100

    # Keep only the samples (columns) that fall below or equal to the threshold
    raw_df = raw_df.loc[:, invalid_pct_per_sample <= sample_nan_pct]

    # if log_norm:
    print("  [Norm] Applying log2(x + 1) normalization...")
    # We add 1 to avoid log(0)
    norm_df = pd.DataFrame(
        np.log2(raw_df + 1), 
        index=raw_df.index, 
        columns=raw_df.columns
    )
    raw_df = pd.DataFrame(
        raw_df, 
        index=raw_df.index, 
        columns=raw_df.columns
    )
    print(f"  [Filter] After sample filtering (≤{sample_nan_pct}% 0/NaN): {raw_df.shape[0]} genes × {raw_df.shape[1]} samples")

    return raw_df,norm_df

# ---------------------------------------------------------------------------
# 4. Main pipeline
# ---------------------------------------------------------------------------


def run_rnaseq_preprocessing():
    os.makedirs(RNASEQ_DATA_DIR, exist_ok=True)
    os.makedirs(RNASEQ_FIGURES, exist_ok=True)

    filter_path = os.path.join(RNASEQ_DATA_DIR, "filter.csv")
    filter_norm_path = os.path.join(RNASEQ_DATA_DIR, "filter_norm.csv")
    combat_path = os.path.join(RNASEQ_DATA_DIR, "combat.csv")
    combat_norm_path = os.path.join(RNASEQ_DATA_DIR, "combat_norm.csv")
    rankin_path = os.path.join(RNASEQ_DATA_DIR, "rankin.csv")

    # ── Stage 1: Filter ──────────────────────────────────────────────────────
    if os.path.exists(filter_path):
        print("Loading cached filter.csv...")
        filtered_df = pd.read_csv(filter_path, index_col=0)
        norm_df = pd.read_csv(filter_norm_path, index_col=0)
    else:
        print(f"Loading combined count matrix from {RNASEQ_COMBINED}...")
        raw_df = pd.read_csv(RNASEQ_COMBINED, index_col=0)
        print(f"  Raw matrix: {raw_df.shape[0]} genes × {raw_df.shape[1]} samples")

        filtered_df,norm_df = run_filtering(raw_df)

        plot_filtering_summary(
            raw_df,
            filtered_df,
            os.path.join(RNASEQ_FIGURES, "filtering_summary.svg"),
        )
        norm_df.to_csv(filter_norm_path)
        filtered_df.to_csv(filter_path)
        print(f"Saved filtered matrix → {filter_path}")
    # return
    # ── Stage 2: ComBat ──────────────────────────────────────────────────
    if os.path.exists(combat_path):
        print("Loading cached combat.csv...")
        combat_df = pd.read_csv(combat_path, index_col=0)
    else:
        print("\nRunning ComBat batch correction (Gaussian, log2 input)...")
        batch_labels = get_batch_labels(norm_df.columns)
        study_counts = Counter(batch_labels)

        # ComBat requires ≥ 2 samples per batch
        single_batches = {s for s, c in study_counts.items() if c < 2}
        if single_batches:
            print(f"  Removing {len(single_batches)} single-sample batches: {single_batches}")

        valid_cols = [c for c, b in zip(norm_df.columns, batch_labels, strict=False)
                    if b not in single_batches]
        valid_batches = [b for b in batch_labels if b not in single_batches]
        df_for_combat = norm_df[valid_cols]

        assert df_for_combat.isna().sum().sum() == 0, \
            f"Unexpected NaNs in ComBat input: {df_for_combat.isna().sum().sum()}"

        combat_df = run_combat(df_for_combat, valid_batches)
        combat_df.to_csv(combat_path)
        print(f"Saved ComBat result → {combat_path}")
    # ── Stage 2.5: ComBat-seq log norm ──────────────────────────────────────────────────
    if os.path.exists(combat_norm_path):
        print("pre existing combat_norm.csv...")
    else:
        print("shifting matrix to avoid clipping...")
        shifted = combat_df - combat_df.values.min()
        log_df = pd.DataFrame(
            np.log1p(shifted.values),
            index=combat_df.index,
            columns=combat_df.columns,
        )
        log_df.to_csv(combat_norm_path)
    # ── Stage 3: Rank-in ─────────────────────────────────────────────────────
    if os.path.exists(rankin_path):
        print("Loading cached rankin.csv...")
        rankin_df = pd.read_csv(rankin_path, index_col=0)
    else:
        print("\nRunning Rank-in normalization on log normalized filter output...")

        # Log1p-transform corrected counts before Rank-in.
        # np.log1p returns an ndarray so we reconstruct the DataFrame explicitly.
        # log_df = pd.DataFrame(
        #     # np.log1p(combat_df.clip(lower=0).values),
        #     np.log1p(norm_df.values),
        #     index=norm_df.index,
        #     columns=norm_df.columns,
        # )

        rankin_df = run_rank_in_normalization(df=norm_df,sample_classes=None, out_path=rankin_path)
        print(f"Saved Rank-in result → {rankin_path}")

    print("\n=== Pipeline Complete ===")
    print(f"  filter.csv    : {filtered_df.shape}")
    print(f"  combat.csv: {combat_df.shape}")
    print(f"  rankin.csv    : {rankin_df.shape}")
    return filtered_df, combat_df, rankin_df


if __name__ == "__main__":
    run_rnaseq_preprocessing()
