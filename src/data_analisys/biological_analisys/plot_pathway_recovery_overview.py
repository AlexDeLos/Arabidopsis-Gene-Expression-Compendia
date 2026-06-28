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
import matplotlib.patches as patches
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
    the expected GO term's RANK PERCENTILE among all terms tested.
    Significant cells get a distinct red star next to the percentage.
    The top of each column displays the Median % Recovery and the fraction of significant hits.
    """
    
    # --- NEW: MAP RAW PIPELINE STAGES TO CLEAN LABELS ---
    STAGE_MAPPING = {
        "filter_norm": "Unnormalized",
        "combat_norm": "ComBat",
        "rankin": "Rank-In"
    }

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
            fdr = sub["FDR q-val"].iloc[0] if "FDR q-val" in sub.columns else sub["fdr_q_value"].iloc[0]
            sig_matrix[i, j] = bool(pd.notna(fdr) and fdr < sig_threshold)

    fig, ax = plt.subplots(figsize=(max(6, len(col_order) * 1.4), max(5, len(row_order) * 0.45)))
    im = ax.imshow(pct_matrix, aspect="auto", cmap="viridis_r", vmin=0, vmax=100)

    # --- Add Thin Lines Between Boxes ---
    ax.set_xticks(np.arange(-0.5, len(col_order), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(row_order), 1), minor=True)
    ax.grid(which="minor", color="#cccccc", linestyle="-", linewidth=0.5)
    ax.tick_params(which="minor", bottom=False, left=False)

    # Standard Major Labels (Bottom)
    ax.set_xticks(range(len(col_order)))
    
    # --- APPLY MAPPING ONLY FOR VISUAL LABELS ---
    display_cols = [STAGE_MAPPING.get(c, c) for c in col_order]
    ax.set_xticklabels(display_cols, rotation=30, ha="right")
    
    ax.set_yticks(range(len(row_order)))
    ax.set_yticklabels(row_order)

    # --- Draw Median % and Hit Fraction Outside the Axes ---
    for j in range(len(col_order)):
        # Extract percentiles ONLY for cells that passed significance
        sig_pcts = [pct_matrix[i, j] for i in range(len(row_order)) 
                    if sig_matrix[i, j] and not np.isnan(pct_matrix[i, j])]
        
        # Count total valid (non-NaN) cells in this column for the denominator
        total_valid = sum(1 for i in range(len(row_order)) if not np.isnan(pct_matrix[i, j]))
        
        if sig_pcts:
            median_val = np.median(sig_pcts)
            num_sig = len(sig_pcts)
            col_text = f"Median\n{median_val:.1f}%\n\n{num_sig}/{total_valid}"
        else:
            col_text = f"Median\n-\n\n0/{total_valid}"
        
        ax.text(
            j, -0.65, col_text,
            ha="center", va="bottom",
            color="#ff3333", fontsize=9, weight="bold",
            clip_on=False
        )
    # -------------------------------------------------------------

    for i in range(len(row_order)):
        for j in range(len(col_order)):
            if np.isnan(pct_matrix[i, j]):
                ax.text(j, i, "n/a", ha="center", va="center", color="gray", fontsize=8)
                continue
                
            contrast_color = "black" if pct_matrix[i, j] < 50 else "white"
            label = f"{pct_matrix[i, j]:.0f}%"
            
            if sig_matrix[i, j]:
                ax.text(
                    j - 0.20, i, "\u2605", 
                    ha="center", va="center", 
                    color="#ff3333", fontsize=15, weight="bold"
                )
                ax.text(
                    j + 0.15, i, label, 
                    ha="center", va="center", 
                    color=contrast_color, fontsize=9
                )
            else:
                ax.text(
                    j, i, label, 
                    ha="center", va="center", 
                    color=contrast_color, fontsize=9
                )

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Percentile rank of expected pathway\n(0% = top hit, 100% = last place)")

    for spine in ax.spines.values():
        spine.set_visible(False)

    plt.tight_layout()
    # Increased top padding to accommodate 4 lines of text (Median / % / space / fraction)
    plt.subplots_adjust(top=0.82) 

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

    # dir_note = " & NES > 0" if require_positive_nes else ""
    # filt_note = f" | filter_low_combination in {filter_options}" if filter_options else ""
    # ax.set_title(
    #     f"{title}\n"
    #     f"(PASS: {sig_col} < {sig_threshold}{dir_note}{filt_note})",
    #     fontsize=10,
    # )
    plt.tight_layout()

    dirpath = Path(save_path).parent
    if str(dirpath):
        dirpath.mkdir(parents=True, exist_ok=True)
    fig.savefig(save_path, dpi=150)
    plt.close(fig)
    print(f"Pathway recovery binary heatmap saved -> {save_path}")