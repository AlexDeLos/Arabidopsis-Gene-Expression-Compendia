"""
plot_pathway_recovery_overview.py
-----------------------------------
Compendium-wide overview plots answering:

    "Is the biology preserved (can I find the expected enriched pathways),
     and how does this preservation change across the 3 normalized data
     matrices and hyperparameters?"

Two public plotting functions, both covering every stress treatment in
STRESS_GO_ROOTS (src/constants_labeling.py) simultaneously:

  plot_pathway_recovery_percentile(...)
      - Heatmap of how strongly each treatment's expected pathway is
        enriched, expressed as a RANK PERCENTILE among all GO terms tested
        for that treatment/matrix. Percentile rank is gene-set-size- and
        term-count-invariant by construction (it's a position, 0-100%,
        not a magnitude), so it is comparable across treatments and across
        normalization methods even when the number of tested terms or
        their gene-set sizes differ.

  plot_pathway_recovery_binary(...)
      - Simple treatment x matrix grid marking PASS/FAIL against a
        significance threshold (FDR q-val by default, configurable),
        with the threshold value always printed on the plot.

Both functions take the SAME long-format recovery table (one row per
treatment x normalization-method-or-hyperparameter-combo), produced by
collect_pathway_recovery_table(...) in this module, which itself reads the
same GSEA result CSVs and the same STRESS_GO_ROOTS mapping used everywhere
else in this pipeline (diff_and_GSEA_pipeline.py, pathway_recovery_analysis.py).

Typical usage, mirroring plot_enrichment_scatter_interactive's calling style:

    recovery_df = collect_pathway_recovery_table(
        base_dir=f"{FIGURES_DIR}GSEA_enrichment_results/",
        experiment_version="matched_control_v1",
        config="All_tissues_full_mixed_min_group_0",
        normalizations=data_types,
    )

    plot_pathway_recovery_percentile(
        recovery_df,
        save_path=f"{FIGURES_DIR}biological_overview_plots/pathway_recovery/recovery_percentile.png",
        title="Pathway recovery (percentile rank)",
        treatments=TREATMENTS,
        normalizations=data_types,
        filter_options=filter_low_combination,
    )

    plot_pathway_recovery_binary(
        recovery_df,
        save_path=f"{FIGURES_DIR}biological_overview_plots/pathway_recovery/recovery_binary.png",
        title="Pathway recovery (pass/fail)",
        treatments=TREATMENTS,
        normalizations=data_types,
        filter_options=filter_low_combination,
        sig_threshold=0.01,
    )
"""

import sys
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

module_dir = "./"
sys.path.append(module_dir)

matplotlib.rc("font", **{"size": 12})

# =============================================================================
# 2. Plot: percentile-rank heatmap
# =============================================================================

def plot_pathway_recovery_percentile(
    recovery_df: pd.DataFrame,
    save_path: str,
    title: str = "Pathway recovery (percentile rank)",
    treatments: list | None = None,
    normalizations: list | None = None,
    filter_options: list[int] | None = None,
    sig_threshold: float = 0.01,
) -> None:
    """
    Heatmap: treatment (rows) x normalization method (columns), colored by
    the expected GO term's RANK PERCENTILE among all terms tested for that
    treatment/matrix (0% = top hit, 100% = last place; lower/darker is
    better - colormap is reversed so strong recovery reads as the most
    visually prominent color).

    Percentile rank is used instead of NES or raw FDR q-val specifically
    because it is invariant to gene-set size and to how many GO terms were
    tested in a given run - a treatment compared against 50 candidate
    terms and one compared against 300 are still both on the same 0-100%
    scale, which is what makes the heatmap meaningfully comparable cell to
    cell. Raw ES/NES do not have this property: their scale depends on
    gene-set size and the shape of the ranked gene list, so two NES values
    of "1.5" are not guaranteed to mean the same thing for two different
    treatments or matrices.

    Cells additionally get a star marker wherever FDR q-val < sig_threshold
    (the conventional significance call), so percentile (continuous
    strength) and significance (discrete pass/fail) are both visible at
    once without needing a second plot for that information.

    Parameters
    ----------
    recovery_df : pd.DataFrame
        Output of collect_pathway_recovery_table(...). Must contain
        columns: treatment, norm_method, percentile, FDR q-val.
    save_path : str
        Output PNG path. Parent directories created as needed.
    title : str
        Plot title.
    treatments, normalizations, filter_options : list or None
        Accepted for calling-convention consistency with
        plot_enrichment_scatter_interactive(...); used here only to fix
        row/column ORDER (rather than relying on whatever order appears
        in recovery_df) so repeated calls with different data still
        produce visually consistent plots. Safe to omit - falls back to
        sorted-unique values from recovery_df itself.
    sig_threshold : float
        FDR q-val threshold for the significance star marker. ALWAYS
        printed in the plot title/subtitle so the reader never has to
        guess what threshold produced the stars.
    """
    row_order = treatments if treatments else sorted(recovery_df["treatment"].unique())
    col_order = normalizations if normalizations else sorted(recovery_df["norm_method"].unique())
    row_order = [t for t in row_order if t in recovery_df["treatment"].unique()]
    col_order = [n for n in col_order if n in recovery_df["norm_method"].unique()]

    pct_matrix = np.full((len(row_order), len(col_order)), np.nan)
    sig_matrix = np.zeros((len(row_order), len(col_order)), dtype=bool)

    for i, t in enumerate(row_order):
        for j, n in enumerate(col_order):
            sub = recovery_df[(recovery_df["treatment"] == t) & (recovery_df["norm_method"] == n)]
            if sub.empty or sub["percentile"].isna().all():
                continue
            pct_matrix[i, j] = sub["percentile"].iloc[0]
            fdr = sub["FDR q-val"].iloc[0]
            sig_matrix[i, j] = bool(pd.notna(fdr) and fdr < sig_threshold)

    fig, ax = plt.subplots(figsize=(max(6, len(col_order) * 1.3), max(5, len(row_order) * 0.45)))
    # Reversed colormap: low percentile (= strong recovery, near rank 1) is the most intense color.
    im = ax.imshow(pct_matrix, aspect="auto", cmap="viridis_r", vmin=0, vmax=100)

    ax.set_xticks(range(len(col_order)))
    ax.set_xticklabels(col_order, rotation=30, ha="right")
    ax.set_yticks(range(len(row_order)))
    ax.set_yticklabels(row_order)

    for i in range(len(row_order)):
        for j in range(len(col_order)):
            if np.isnan(pct_matrix[i, j]):
                ax.text(j, i, "n/a", ha="center", va="center", color="gray", fontsize=8)
                continue
            label = f"{pct_matrix[i, j]:.0f}%"
            if sig_matrix[i, j]:
                label = "\u2605 " + label
            ax.text(j, i, label, ha="center", va="center",
                     color="white" if pct_matrix[i, j] < 50 else "black", fontsize=9)

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Percentile rank of expected pathway\n(0% = top hit, 100% = last place)")

    filt_note = f" | filter_low_combination in {filter_options}" if filter_options else ""
    ax.set_title(
        f"{title}\n"
        f"(\u2605 = FDR q-val < {sig_threshold}{filt_note})",
        fontsize=10,
    )
    plt.tight_layout()

    dirpath = Path(save_path).parent
    if str(dirpath):
        dirpath.mkdir(parents=True, exist_ok=True)
    fig.savefig(save_path, dpi=150)
    plt.close(fig)
    print(f"Pathway recovery percentile heatmap saved -> {save_path}")


# =============================================================================
# 3. Plot: binary pass/fail heatmap
# =============================================================================

def plot_pathway_recovery_binary(
    recovery_df: pd.DataFrame,
    save_path: str,
    title: str = "Pathway recovery (pass/fail)",
    treatments: list | None = None,
    normalizations: list | None = None,
    filter_options: list[int] | None = None,
    sig_threshold: float = 0.01,
    sig_col: str = "FDR q-val",
    require_positive_nes: bool = True,
) -> None:
    """
    Simple treatment x normalization-method grid marking PASS (expected
    pathway significantly enriched, correct direction) vs. FAIL, against
    `sig_threshold` on `sig_col`. The threshold value is always rendered
    directly in the plot title - never implicit - per the requirement
    that any plotted significance call must show its own threshold.

    This plot intentionally throws away the continuous strength
    information that plot_pathway_recovery_percentile shows - it answers
    only "did this matrix recover this pathway, yes or no," which is the
    simplest possible summary for a reader who wants the headline result
    before looking at magnitude/rank detail.

    Parameters
    ----------
    recovery_df : pd.DataFrame
        Output of collect_pathway_recovery_table(...). Must contain
        columns: treatment, norm_method, NES, and `sig_col`.
    save_path : str
        Output PNG path.
    title : str
        Plot title.
    treatments, normalizations, filter_options : list or None
        Same role as in plot_pathway_recovery_percentile - fixes
        row/column order and is echoed into the subtitle.
    sig_threshold : float
        Significance threshold on `sig_col`. ALWAYS shown in the title.
    sig_col : str
        Which column to threshold against - "FDR q-val" by default. Pass
        "NOM p-val" or "FWER p-val" if you want a different significance
        criterion; the threshold and column name are both rendered in the
        title regardless of which is chosen, so the reader always knows
        exactly what test produced each PASS/FAIL cell.
    require_positive_nes : bool
        If True (default), a cell only counts as PASS when, in addition
        to clearing `sig_threshold`, NES > 0 (enrichment in the expected
        direction, not depletion). Set False if direction shouldn't matter
        for your use case (e.g. you intentionally expect some treatments,
        such as dark/etiolated "Low Light" samples, to show depletion
        rather than enrichment of the mapped term).
    """
    row_order = treatments if treatments else sorted(recovery_df["treatment"].unique())
    col_order = normalizations if normalizations else sorted(recovery_df["norm_method"].unique())
    row_order = [t for t in row_order if t in recovery_df["treatment"].unique()]
    col_order = [n for n in col_order if n in recovery_df["norm_method"].unique()]

    pass_matrix = np.zeros((len(row_order), len(col_order)), dtype=int)  # 0=n/a, 1=fail, 2=pass

    for i, t in enumerate(row_order):
        for j, n in enumerate(col_order):
            sub = recovery_df[(recovery_df["treatment"] == t) & (recovery_df["norm_method"] == n)]
            if sub.empty or sub[sig_col].isna().all():
                pass_matrix[i, j] = 0
                continue
            sig_val = sub[sig_col].iloc[0]
            nes_val = sub["NES"].iloc[0] if "NES" in sub.columns else np.nan
            passes_sig = pd.notna(sig_val) and sig_val < sig_threshold
            passes_dir = (not require_positive_nes) or (pd.notna(nes_val) and nes_val > 0)
            pass_matrix[i, j] = 2 if (passes_sig and passes_dir) else 1

    fig, ax = plt.subplots(figsize=(max(6, len(col_order) * 1.3), max(5, len(row_order) * 0.45)))

    cmap = matplotlib.colors.ListedColormap(["#dddddd", "#c0392b", "#27ae60"])
    bounds = [-0.5, 0.5, 1.5, 2.5]
    bnorm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
    ax.imshow(pass_matrix, aspect="auto", cmap=cmap, norm=bnorm)

    ax.set_xticks(range(len(col_order)))
    ax.set_xticklabels(col_order, rotation=30, ha="right")
    ax.set_yticks(range(len(row_order)))
    ax.set_yticklabels(row_order)

    for i in range(len(row_order)):
        for j in range(len(col_order)):
            label = {0: "n/a", 1: "FAIL", 2: "PASS"}[pass_matrix[i, j]]
            ax.text(j, i, label, ha="center", va="center",
                     color="black" if pass_matrix[i, j] == 0 else "white", fontsize=9, fontweight="bold")

    dir_note = " & NES > 0" if require_positive_nes else ""
    filt_note = f" | filter_low_combination in {filter_options}" if filter_options else ""
    ax.set_title(
        f"{title}\n"
        f"(PASS: {sig_col} < {sig_threshold}{dir_note}{filt_note})",
        fontsize=10,
    )
    plt.tight_layout()

    dirpath = Path(save_path).parent
    if str(dirpath):
        dirpath.mkdir(parents=True, exist_ok=True)
    fig.savefig(save_path, dpi=150)
    plt.close(fig)
    print(f"Pathway recovery binary heatmap saved -> {save_path}")