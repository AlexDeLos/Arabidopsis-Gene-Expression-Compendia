"""
gsea_normalization_impact.py
-----------------------------
Quantifies how much normalization method (combat_norm, rankin vs.
filter_norm baseline) and tissue scope (All_tissues vs. leaf) shift each
GSEA stat (ES, NES, NOM p-val, FDR q-val, FWER p-val), expressed as
percent change relative to the filter_norm / All_tissues baseline.

For each treatment and each stat, only GO terms present in BOTH the
baseline and the comparison run are used (paired comparison) — this
avoids comparing apples to oranges when a term is significant in one
run but wasn't even tested/returned in the other.

Produces one boxplot per stat, showing the distribution of percent
changes across all GO terms within each treatment (not collapsed to a
single average), so you can see both the typical shift AND how much it
varies treatment to treatment.

Usage:
    python gsea_normalization_impact.py
    python gsea_normalization_impact.py --config-b leaf_full_mixed_min_group_0 --axis tissue
    python gsea_normalization_impact.py --out-dir ./plots/
"""

import argparse
from pathlib import Path
import sys
import matplotlib
matplotlib.use("Agg")  # headless — no display on the cluster
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# =============================================================================
# Paths
# =============================================================================

BASE = Path(
    "/tudelft.net/staff-umbrella/GeneExpressionStorage/outputs/5.0/GSEA_enrichment_results"
)
EXPERIMENT_NAME = "5.0"
PERMUTATIONS = 1000

STAT_COLS = ["ES", "NES", "NOM p-val", "FDR q-val", "FWER p-val"]

# Stats that can cross zero — percent change is unstable near zero
# denominators for these, so we flag (not silently drop) extreme outliers.
SIGNED_STATS = {"ES", "NES"}
EXTREME_PCT_THRESHOLD = 500.0  # percent change beyond this is flagged as unstable

module_dir = "./"
sys.path.append(module_dir)
from src.constants_labeling import TreatmentEnum  # noqa: E402



NON_STRESS_LABELS = {
    TreatmentEnum.OTHER,
    TreatmentEnum.CONTROL,
    TreatmentEnum.UNKNOWN,
    TreatmentEnum.UNSPECIFIED,
}
STRESS_TREATMENTS = [t for t in TreatmentEnum if t not in NON_STRESS_LABELS]


def build_csv_path(norm: str, config: str, treatment: str) -> Path:
    folder = f"GSEA_enrichment_{EXPERIMENT_NAME}_{norm}_{config}"
    filename = f"{treatment}_gsea_go_enrichment_results_{PERMUTATIONS}.csv"
    return BASE / folder / filename


def load_csv(norm: str, config: str, treatment: str) -> pd.DataFrame | None:
    path = build_csv_path(norm, config, treatment)
    if not path.exists():
        return None
    df = pd.read_csv(path)
    missing = [c for c in STAT_COLS + ["go_id"] if c not in df.columns]
    if missing:
        print(f"    [warn] {path} missing columns {missing}, skipping")
        return None
    return df


def paired_percent_change(
    baseline_df: pd.DataFrame, compare_df: pd.DataFrame, treatment: str, comparison_label: str
) -> pd.DataFrame:
    """
    Build a long-format table of percent change per matched go_id, per stat,
    for one treatment and one baseline-vs-comparison pairing.
    """
    merged = baseline_df[["go_id"] + STAT_COLS].merge(
        compare_df[["go_id"] + STAT_COLS],
        on="go_id",
        suffixes=("_base", "_cmp"),
    )

    rows = []
    n_extreme = 0
    n_undefined = 0
    for stat in STAT_COLS:
        base_vals = merged[f"{stat}_base"]
        cmp_vals = merged[f"{stat}_cmp"]

        with np.errstate(divide="ignore", invalid="ignore"):
            pct_change = (cmp_vals - base_vals) / base_vals * 100.0

        undefined = ~np.isfinite(pct_change)
        n_undefined += int(undefined.sum())
        pct_change = pct_change.where(~undefined, np.nan)

        if stat in SIGNED_STATS:
            extreme = pct_change.abs() > EXTREME_PCT_THRESHOLD
            n_extreme += int(extreme.sum())

        for go_id, val in zip(merged["go_id"], pct_change):
            rows.append(
                {
                    "treatment": treatment,
                    "comparison": comparison_label,
                    "stat": stat,
                    "go_id": go_id,
                    "pct_change": val,
                }
            )

    if n_undefined:
        print(
            f"    [note] {treatment} / {comparison_label}: {n_undefined} "
            f"undefined pct-change values (zero baseline) excluded from plot"
        )
    if n_extreme:
        print(
            f"    [note] {treatment} / {comparison_label}: {n_extreme} "
            f"extreme pct-change values (>|{EXTREME_PCT_THRESHOLD}%|) in ES/NES "
            f"— likely near-zero baseline denominators, interpret with caution"
        )

    return pd.DataFrame(rows)


def collect_all(config_a: str, config_b: str, axis_b_label: str) -> pd.DataFrame:
    """
    For every treatment, build percent-change rows for:
      - combat_norm vs filter_norm   (same config_a, normalization axis)
      - rankin      vs filter_norm   (same config_a, normalization axis)
      - filter_norm config_b vs filter_norm config_a   (tissue-scope axis,
        using filter_norm on both sides so normalization isn't a confound)
    """
    all_rows = []

    for treatment in STRESS_TREATMENTS:
        t = treatment.value
        print(f"--- {t} ---")

        filt_a = load_csv("filter_norm", config_a, t)
        if filt_a is None:
            print(f"    [skip] no filter_norm baseline for {t} / {config_a}")
            continue

        combat_a = load_csv("combat_norm", config_a, t)
        rankin_a = load_csv("rankin", config_a, t)
        filt_b = load_csv("filter_norm", config_b, t)

        if combat_a is not None:
            all_rows.append(paired_percent_change(filt_a, combat_a, t, "combat_norm vs filter_norm"))
        else:
            print(f"    [skip] no combat_norm for {t} / {config_a}")

        if rankin_a is not None:
            all_rows.append(paired_percent_change(filt_a, rankin_a, t, "rankin vs filter_norm"))
        else:
            print(f"    [skip] no rankin for {t} / {config_a}")

        if filt_b is not None:
            all_rows.append(paired_percent_change(filt_a, filt_b, t, f"{axis_b_label} vs filter_norm"))
        else:
            print(f"    [skip] no filter_norm for {t} / {config_b}")

    if not all_rows:
        return pd.DataFrame()
    return pd.concat(all_rows, ignore_index=True)


def plot_stat(df: pd.DataFrame, stat: str, comparison: str, out_dir: Path) -> list[Path]:
    """
    One boxplot: x-axis = treatment, y-axis = % change for this stat/comparison.
    For signed stats (ES, NES), also writes a y-axis-clipped variant, since a
    handful of near-zero-baseline outliers can otherwise crush the visible
    box into an unreadable sliver near zero.
    """
    sub = df[(df["stat"] == stat) & (df["comparison"] == comparison)].dropna(subset=["pct_change"])
    if sub.empty:
        return []

    treatments = [t.value for t in STRESS_TREATMENTS if t.value in sub["treatment"].unique()]
    data = [sub.loc[sub["treatment"] == t, "pct_change"].values for t in treatments]

    safe_stat = stat.replace(" ", "_").replace(".", "")
    safe_cmp = comparison.replace(" ", "_")
    written = []

    def _make_plot(ylim: tuple[float, float] | None, suffix: str) -> Path:
        fig, ax = plt.subplots(figsize=(max(8, len(treatments) * 0.8), 5))
        ax.boxplot(data, labels=treatments, showfliers=True)
        ax.axhline(0, color="gray", linestyle="--", linewidth=1)
        if ylim:
            ax.set_ylim(*ylim)
        ax.set_ylabel(f"% change in {stat}\n({comparison})")
        title_note = " (y-axis clipped — outliers exist beyond range)" if ylim else ""
        ax.set_title(f"{stat}: {comparison}{title_note}\n"
                      f"(per-term % change, paired by go_id, filter_norm baseline)")
        ax.tick_params(axis="x", rotation=45)
        plt.tight_layout()
        out_path = out_dir / f"pct_change_{safe_stat}_{safe_cmp}{suffix}.pdf"
        fig.savefig(out_path, dpi=150)
        plt.close(fig)
        return out_path

    written.append(_make_plot(ylim=None, suffix=""))

    if stat in SIGNED_STATS:
        # Clip to ~1.5x the IQR-based whisker range across all treatments,
        # so the readable box-and-whisker structure isn't dwarfed by a few
        # near-zero-denominator outliers.
        all_vals = np.concatenate(data) if data else np.array([])
        finite_vals = all_vals[np.isfinite(all_vals)]
        if finite_vals.size:
            q1, q3 = np.percentile(finite_vals, [25, 75])
            iqr = q3 - q1
            pad = max(iqr * 3, 50)  # floor so near-zero IQR doesn't over-clip
            lo, hi = q1 - pad, q3 + pad
            written.append(_make_plot(ylim=(lo, hi), suffix="_clipped"))

    return written


def main():
    parser = argparse.ArgumentParser(
        description="Plot per-term percent change in GSEA stats vs. filter_norm baseline, "
                     "across normalization methods and tissue scope."
    )
    parser.add_argument("--config-a", default="All_tissues_full_mixed_min_group_0",
                         help="Baseline config (filter_norm here = the reference point)")
    parser.add_argument("--config-b", default="leaf_full_mixed_min_group_0",
                         help="Comparison config for the tissue-scope axis")
    parser.add_argument("--axis-b-label", default="leaf",
                         help="Label for the tissue-scope comparison (e.g. 'leaf')")
    parser.add_argument("--out-dir", default="./gsea_norm_impact_plots",
                         help="Directory to write PNG plots and the combined CSV into")
    parser.add_argument("--out-csv", default="gsea_pct_change.csv",
                         help="Filename (within --out-dir) for the combined long-format CSV")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Baseline config: filter_norm / {args.config_a}")
    print(f"Tissue-scope comparison: filter_norm / {args.config_b} ({args.axis_b_label})\n")

    df = collect_all(args.config_a, args.config_b, args.axis_b_label)

    if df.empty:
        print("\nNo data collected — nothing to plot.")
        return

    csv_path = out_dir / args.out_csv
    df.to_csv(csv_path, index=False)
    print(f"\nWrote combined long-format percent-change table ({len(df)} rows) → {csv_path}")

    comparisons = df["comparison"].unique()
    print(f"\nGenerating boxplots for {len(STAT_COLS)} stats × {len(comparisons)} comparisons...")
    written = []
    for stat in STAT_COLS:
        for comparison in comparisons:
            paths = plot_stat(df, stat, comparison, out_dir)
            for path in paths:
                written.append(path)
                print(f"  wrote {path}")

    print(f"\nDone. {len(written)} plots written to {out_dir}")


if __name__ == "__main__":
    main()
