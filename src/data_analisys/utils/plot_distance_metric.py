"""
plot_distance_metrics
======================
Visualizes the output of run_distance_evaluation() collected across all
normalization stages into `all_dist_metrics`.

Usage
-----
    from src.data_analisys.plot_distance_metrics import plot_distance_metrics

    plot_distance_metrics(
        all_dist_metrics=all_dist_metrics,     # {stage_name: pd.DataFrame | dict}
        output_folder=comparison_output_dir,
        experiment_name="Distance_Metrics",
        show=True,
    )
"""

from __future__ import annotations

import os
from typing import Dict, Union

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Color palette (matches existing exploration figures)
# ---------------------------------------------------------------------------
_COL_INTER  = "#378ADD"   # blue  — inter-study
_COL_INTRA  = "#1D9E75"   # teal  — intra-study
_COL_GOOD   = "#1D9E75"   # green — positive separation
_COL_BAD    = "#E24B4A"   # red   — negative separation
_COL_REF    = "#BA7517"   # amber — zero reference
_ALPHA      = 0.75


def _extract_row(record: Union[pd.DataFrame, dict]) -> dict:
    """Accept either a single-row DataFrame or a plain dict of metric values."""
    if isinstance(record, pd.DataFrame):
        return record.iloc[0].to_dict()
    return record


def plot_distance_metrics(
    all_dist_metrics: Dict[str, Union[pd.DataFrame, dict]],
    output_folder: str | None = None,
    experiment_name: str = "Distance_Metrics",
    show: bool = False,
    figsize: tuple[float, float] = (16, 14),
    plot_ratio: bool = False,
) -> plt.Figure:

    plt.rcParams.update(
        {
            "font.size": 28,
            "axes.titlesize": 30,
            "axes.labelsize": 28,
            "xtick.labelsize": 24,
            "ytick.labelsize": 24,
            "legend.fontsize": 22,
        }
    )

    if not all_dist_metrics:
        raise ValueError("all_dist_metrics is empty.")

    # --------------------------------------------------
    # Extract data
    # --------------------------------------------------

    stages = list(all_dist_metrics.keys())
    rows = [_extract_row(all_dist_metrics[s]) for s in stages]

    inter = np.array(
        [r["G_d_inter"] for r in rows],
        dtype=float,
    )

    intra = np.array(
        [r["G_d_intra"] for r in rows],
        dtype=float,
    )

    sep = np.array(
        [r["BiologicalSeparation"] for r in rows],
        dtype=float,
    )

    spearman = np.array(
        [
            r.get(
                "SimilarityDistanceSpearman",
                np.nan,
            )
            for r in rows
        ],
        dtype=float,
    )

    spearman_p = np.array(
        [
            r.get(
                "SimilarityDistanceSpearmanP",
                np.nan,
            )
            for r in rows
        ],
        dtype=float,
    )

    x = np.arange(len(stages))
    width = 0.35

    sep_colors = [
        _COL_GOOD if s > 0 else _COL_BAD
        for s in sep
    ]

    # --------------------------------------------------
    # Layout
    # --------------------------------------------------

    if plot_ratio:

        fig, axes = plt.subplots(
            3,
            1,
            figsize=(figsize[0], figsize[1] * 1.5),
            gridspec_kw={
                "height_ratios": [3, 2, 2]
            },
        )

        ax_dist = axes[0]
        ax_sep = axes[1]
        ax_corr = axes[2]

    else:

        fig, axes = plt.subplots(
            2,
            1,
            figsize=(figsize[0], figsize[1]),
            gridspec_kw={
                "height_ratios": [3, 2]
            },
        )

        ax_dist = axes[0]
        ax_corr = axes[1]

    fig.patch.set_facecolor("white")

    # ==================================================
    # DISTANCE METRICS
    # ==================================================

    ax_dist.set_facecolor("#fafafa")

    bars_inter = ax_dist.bar(
        x - width / 2,
        inter,
        width,
        label=r"$G_{d,\mathrm{inter}}$",
        color=_COL_INTER,
        alpha=_ALPHA,
        edgecolor=_COL_INTER,
        linewidth=1.5,
    )

    bars_intra = ax_dist.bar(
        x + width / 2,
        intra,
        width,
        label=r"$G_{d,\mathrm{intra}}$",
        color=_COL_INTRA,
        alpha=_ALPHA,
        edgecolor=_COL_INTRA,
        linewidth=1.5,
        hatch="//",
    )

    ymax = np.nanmax(
        np.concatenate([inter, intra])
    )

    for bars in (bars_inter, bars_intra):

        for bar in bars:

            h = bar.get_height()

            if np.isfinite(h):

                ax_dist.text(
                    bar.get_x() + bar.get_width() / 2,
                    h + ymax * 0.02,
                    f"{h:.3f}",
                    ha="center",
                    va="bottom",
                    fontsize=18,
                )

    ax_dist.set_xticks(x)

    ax_dist.set_xticklabels(
        stages,
        rotation=15,
        ha="right",
    )

    ax_dist.set_ylabel(
        "Normalized weighted distance"
    )

    ax_dist.set_title(
        "Inter-study and intra-study biological distances",
        pad=20,
    )

    ax_dist.grid(
        axis="y",
        linestyle=":",
        linewidth=1,
        color="#cccccc",
    )

    ax_dist.spines[
        ["top", "right"]
    ].set_visible(False)

    ax_dist.legend(
        framealpha=0.8,
    )

    # ==================================================
    # SEPARATION SCORE
    # ==================================================

    if plot_ratio:

        ax_sep.set_facecolor("#fafafa")

        bars_sep = ax_sep.bar(
            x,
            sep,
            width=0.7,
            color=[c + "bb" for c in sep_colors],
            edgecolor=sep_colors,
            linewidth=1.5,
        )

        ax_sep.axhline(
            0,
            color=_COL_REF,
            linestyle="--",
            linewidth=2,
        )

        for bar, value in zip(
            bars_sep,
            sep,
        ):

            if np.isfinite(value):

                ax_sep.text(
                    bar.get_x()
                    + bar.get_width() / 2,
                    value,
                    f"{value:.3f}",
                    ha="center",
                    va="bottom"
                    if value >= 0
                    else "top",
                    fontsize=20,
                    fontweight="bold",
                )

        ax_sep.set_xticks(x)

        ax_sep.set_xticklabels(
            stages,
            rotation=15,
            ha="right",
        )

        ax_sep.set_ylabel(
            "Separation score"
        )

        ax_sep.set_title(
            "Biological separation metric"
        )

        ax_sep.grid(
            axis="y",
            linestyle=":",
            linewidth=1,
            color="#cccccc",
        )

        ax_sep.spines[
            ["top", "right"]
        ].set_visible(False)

    # ==================================================
    # SIMILARITY-DISTANCE CORRELATION
    # ==================================================

    ax_corr.set_facecolor("#fafafa")

    corr_colors = [
        _COL_GOOD if c < 0 else _COL_BAD
        for c in spearman
    ]

    bars_corr = ax_corr.bar(
        x,
        spearman,
        width=0.7,
        color=[c + "bb" for c in corr_colors],
        edgecolor=corr_colors,
        linewidth=1.5,
    )

    ax_corr.axhline(
        0,
        color=_COL_REF,
        linestyle="--",
        linewidth=2,
    )

    for bar, corr, pval in zip(
        bars_corr,
        spearman,
        spearman_p,
    ):

        if np.isfinite(corr):

            label = (
                f"ρ={corr:.3f}\n"
                f"p={pval:.1e}"
            )

            ax_corr.text(
                bar.get_x()
                + bar.get_width() / 2,
                corr,
                label,
                ha="center",
                va="bottom"
                if corr >= 0
                else "top",
                fontsize=18,
            )

    ax_corr.set_xticks(x)

    ax_corr.set_xticklabels(
        stages,
        rotation=15,
        ha="right",
    )

    ax_corr.set_ylabel(
        "Spearman correlation"
    )

    ax_corr.set_title(
        "Similarity vs distance relationship"
    )

    ax_corr.grid(
        axis="y",
        linestyle=":",
        linewidth=1,
        color="#cccccc",
    )

    ax_corr.spines[
        ["top", "right"]
    ].set_visible(False)

    # ==================================================
    # FINALIZE
    # ==================================================

    plt.tight_layout(
        h_pad=3.0
    )

    if output_folder is not None:

        os.makedirs(
            output_folder,
            exist_ok=True,
        )

        out_path = os.path.join(
            output_folder,
            f"{experiment_name}.pdf",
        )

        fig.savefig(
            out_path,
            dpi=300,
            bbox_inches="tight",
        )

        print(
            f"[DistMetrics] Figure saved → {out_path}"
        )

    if show:
        plt.show()

    return fig

def plot_similarity_distance_scatter(
    pairwise_df: pd.DataFrame,
    output_folder: str | None = None,
    experiment_name: str = "SimilarityDistance",
    show: bool = False,
    max_points: int = 100_000,
) -> plt.Figure:

    from scipy.stats import spearmanr

    df = pairwise_df.copy()

    if len(df) > max_points:
        df = df.sample(
            max_points,
            random_state=42,
        )

    corr, pval = spearmanr(
        df["Distance"],
        df["Similarity"],
    )

    fig, ax = plt.subplots(
        figsize=(10, 8)
    )

    ax.scatter(
        df["Distance"],
        df["Similarity"],
        alpha=0.05,
        s=4,
    )

    ax.set_xlabel(
        "PCA distance"
    )

    ax.set_ylabel(
        "Biological similarity"
    )

    ax.set_title(
        (
            "Similarity vs Distance\n"
            f"Spearman = {corr:.3f} "
            f"(p={pval:.2e})"
        )
    )

    ax.grid(
        linestyle=":",
        alpha=0.5,
    )

    if output_folder is not None:

        os.makedirs(
            output_folder,
            exist_ok=True,
        )

        path = os.path.join(
            output_folder,
            f"{experiment_name}.pdf",
        )

        fig.savefig(
            path,
            dpi=300,
            bbox_inches="tight",
        )

        print(
            f"[DistMetrics] Scatter saved → {path}"
        )

    if show:
        plt.show()

    return fig