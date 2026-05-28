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
_COL_GOOD   = "#1D9E75"   # green — ratio > 1
_COL_BAD    = "#E24B4A"   # red   — ratio ≤ 1
_COL_REF    = "#BA7517"   # amber — ratio = 1 reference line
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
    figsize: tuple[float, float] = (10, 8),
    plot_ratio: bool = False,
) -> plt.Figure:
    """
    Plot G_d_inter, G_d_intra, and optionally Ratio_global for every evaluated stage.

    Parameters
    ----------
    all_dist_metrics : dict
        {stage_name: single-row DataFrame or dict} as produced by
        run_distance_evaluation() collected in a loop.
    output_folder : str, optional
        Directory to save the figure.  Skipped if None.
    experiment_name : str
        Base filename (no extension) for the saved figure.
    show : bool
        Call plt.show() after plotting.
    figsize : tuple
        Overall figure size in inches.
    plot_ratio : bool
        If True (default), include the Ratio_global panel below the distance
        bars.  If False, only the distance bar chart is shown.

    Returns
    -------
    matplotlib.figure.Figure
    """
    if not all_dist_metrics:
        raise ValueError("all_dist_metrics is empty.")

    # ------------------------------------------------------------------ #
    # Unpack data
    # ------------------------------------------------------------------ #
    stages = list(all_dist_metrics.keys())
    rows   = [_extract_row(all_dist_metrics[s]) for s in stages]

    inter  = np.array([r["G_d_inter"]    for r in rows], dtype=float)
    intra  = np.array([r["G_d_intra"]    for r in rows], dtype=float)
    ratio  = np.array([r["Ratio_global"] for r in rows], dtype=float)

    x      = np.arange(len(stages))
    width  = 0.35

    ratio_colors = [_COL_GOOD if r > 1 else _COL_BAD for r in ratio]

    # ------------------------------------------------------------------ #
    # Layout: 1 or 2 rows depending on plot_ratio
    # ------------------------------------------------------------------ #
    if plot_ratio:
        fig, axes = plt.subplots(
            2, 1,
            figsize=figsize,
            gridspec_kw={"height_ratios": [3, 2]},
        )
        ax_dist = axes[0]
        ax_rat  = axes[1]
    else:
        fig, ax_dist = plt.subplots(1, 1, figsize=figsize)
    fig.patch.set_facecolor("white")

    # ── Top panel: grouped bar chart ────────────────────────────────── #
    ax_dist = axes[0]
    ax_dist.set_facecolor("#fafafa")

    bars_inter = ax_dist.bar(
        x - width / 2, inter, width,
        label="$G_{d,\\mathrm{inter}}$  (inter-study)",
        color=_COL_INTER, alpha=_ALPHA, edgecolor=_COL_INTER, linewidth=1,
    )
    bars_intra = ax_dist.bar(
        x + width / 2, intra, width,
        label="$G_{d,\\mathrm{intra}}$  (intra-study)",
        color=_COL_INTRA, alpha=_ALPHA, edgecolor=_COL_INTRA, linewidth=1,
        hatch="//",
    )

    # value labels on bars
    for bars in (bars_inter, bars_intra):
        for bar in bars:
            h = bar.get_height()
            if np.isfinite(h):
                ax_dist.text(
                    bar.get_x() + bar.get_width() / 2,
                    h + 0.005 * max(inter.max(), intra.max()),
                    f"{h:.3f}",
                    ha="center", va="bottom", fontsize=8, color="#444",
                )

    ax_dist.set_xticks(x)
    ax_dist.set_xticklabels(stages, fontsize=10, rotation=15, ha="right")
    ax_dist.set_ylabel("Weighted distance (PCA space)", fontsize=10)
    ax_dist.set_title("Normalization stage evaluation — distance metrics", fontsize=12, pad=10)
    ax_dist.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.3f"))
    ax_dist.legend(fontsize=9, framealpha=0.6)
    ax_dist.spines[["top", "right"]].set_visible(False)
    ax_dist.grid(axis="y", linestyle=":", linewidth=0.6, color="#ccc")

    # ── Bottom panel: ratio bar chart ───────────────────────────────── #
    if plot_ratio:
        ax_rat.set_facecolor("#fafafa")

        bars_ratio = ax_rat.bar(
            x, ratio, width * 1.6,
            color=[c + "bb" for c in ratio_colors],   # hex + alpha suffix
            edgecolor=ratio_colors, linewidth=1.2,
        )

        # reference line at 1.0
        ax_rat.axhline(1.0, color=_COL_REF, linewidth=1.5, linestyle="--", zorder=3)
        ax_rat.text(
            len(stages) - 0.5, 1.0,
            " ratio = 1.0",
            va="bottom", ha="right", fontsize=8, color=_COL_REF,
        )

        # value labels
        for bar, r in zip(bars_ratio, ratio):
            if np.isfinite(r):
                ax_rat.text(
                    bar.get_x() + bar.get_width() / 2,
                    r + 0.01,
                    f"{r:.3f}",
                    ha="center", va="bottom", fontsize=9, fontweight="bold",
                    color=_COL_GOOD if r > 1 else _COL_BAD,
                )

        ax_rat.set_xticks(x)
        ax_rat.set_xticklabels(stages, fontsize=10, rotation=15, ha="right")
        ax_rat.set_ylabel("$\\mathrm{Ratio_{global}}$ = intra / inter", fontsize=10)
        ax_rat.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.2f"))
        ax_rat.spines[["top", "right"]].set_visible(False)
        ax_rat.grid(axis="y", linestyle=":", linewidth=0.6, color="#ccc")

        # coloured legend for pass/fail
        from matplotlib.patches import Patch
        ax_rat.legend(
            handles=[
                Patch(color=_COL_GOOD, alpha=0.7, label="Ratio > 1  ✓  biology > batch"),
                Patch(color=_COL_BAD,  alpha=0.7, label="Ratio ≤ 1  ✗  batch dominates"),
            ],
            fontsize=8, framealpha=0.6, loc="upper left",
        )

    plt.tight_layout(h_pad=2.5)

    # ── Save ───────────────────────────────────────────────────────── #
    if output_folder is not None:
        os.makedirs(output_folder, exist_ok=True)
        out_path = os.path.join(output_folder, f"{experiment_name}.png")
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        print(f"[DistMetrics] Figure saved → {out_path}")

    if show:
        plt.show()

    return fig