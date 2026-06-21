"""
pathway_recovery_analysis.py
------------------------------
Biological validation analysis (Results section, Priority 1):

  "Does the compendium recover known Arabidopsis stress-response pathways?"

For each stress treatment, this script checks whether the literature-expected
GO term (e.g. Cold -> "response to cold", GO:0009409) is recovered by GSEA —
i.e. shows up as significantly enriched, in the correct direction — and how
that recovery compares across normalization methods (filter_norm, combat_norm,
rankin).

This reuses the SAME GO root mapping, fallback logic, and output paths as the
rest of the pipeline (diff_and_GSEA_pipeline.py / pr_rank_gene_enrich.py /
comapre_gsea.py), so results here are directly consistent with everything
already generated — no separate GO resolution step, no separate path scheme.

Outputs
-------
1. A combined long-format CSV: one row per (treatment, norm_method), with
   the expected GO term's ES/NES/FDR q-val/rank-among-tested-terms/hit flag.
2. A recovery heatmap: treatment (rows) x normalization method (columns),
   colored by -log10(FDR q-val) of the expected term, annotated with a hit
   marker when recovery is confirmed (FDR < threshold AND NES > 0).
3. A per-treatment rank strip-plot showing where the expected term landed
   among all tested GO terms, per normalization method.
4. Console-printed provenance for every value used, so a reviewer (or future
   you) can trace exactly which file/row produced each number in the figures.

Usage:
    python pathway_recovery_analysis.py
    python pathway_recovery_analysis.py --config All_tissues_full_mixed_min_group_0
    python pathway_recovery_analysis.py --fdr 0.1 --out-dir ./pathway_recovery_plots
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # headless — no display on the cluster
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
module_dir = "./"
sys.path.append(module_dir)
from src.constants_labeling import TreatmentEnum,STRESS_GO_ROOTS  # noqa: E402
# =============================================================================
# Paths — same scheme used throughout the pipeline (comapre_gsea.py,
# gsea_normalization_impact.py). Change BASE here once if the output root
# ever moves; everything else in this script derives from it.
# =============================================================================

BASE = Path(
    "/tudelft.net/staff-umbrella/GeneExpressionStorage/outputs/5.0/GSEA_enrichment_results"
)
EXPERIMENT_NAME = "5.0"
PERMUTATIONS = 1000
NORM_METHODS = ["filter_norm", "combat_norm", "rankin"]
STAT_COLS = ["ES", "NES", "NOM p-val", "FDR q-val", "FWER p-val"]


# =============================================================================
# Treatment enum — identical to the rest of the pipeline / comparison scripts
# =============================================================================

NON_STRESS_LABELS = {
    TreatmentEnum.OTHER,
    TreatmentEnum.CONTROL,
    TreatmentEnum.UNKNOWN,
    TreatmentEnum.UNSPECIFIED,
}
STRESS_TREATMENTS = [t for t in TreatmentEnum if t not in NON_STRESS_LABELS]


# =============================================================================
# GO root mapping — finalized version from the GO_NAME_OVERRIDES debugging
# session: Drought/Flood corrected away from bad fuzzy-matcher hits
# (Drought was originally mismatched to "response to DDT"; Flood's literal
# "response to flooding" term has zero annotated genes, so it is anchored to
# the broader, populated "response to decreased oxygen levels" parent term
# instead, since flooding stress in plants is mechanistically a hypoxia
# response). All other entries are the original fuzzy-matched names that
# were confirmed correct (score >= 0.847) earlier in this project.
#
# IMPORTANT: keep this dict in sync with GO_NAME_OVERRIDES in
# diff_and_GSEA_pipeline.py if that ever changes — this script intentionally
# does NOT import it directly, since the cluster entrypoint for this script
# may run in a lighter-weight environment without the full pipeline's
# dependencies (goatools/gseapy) loaded at import time. See the GO ID
# resolution note printed at runtime for how to re-verify these.
# =============================================================================



# Flagged separately so it's impossible to miss in the console output —
# this mapping was identified as biologically wrong (a cholesterol-drug
# response term, not a salt-stress term) but never corrected with a
# verified replacement during the original debugging session.
# KNOWN_BAD_GO_MAPPINGS = {"Salinity"}


def build_csv_path(norm: str, config: str, treatment: str) -> Path:
    """Reconstruct the GSEA results CSV path — same scheme as comapre_gsea.py."""
    folder = f"GSEA_enrichment_{EXPERIMENT_NAME}_{norm}_{config}"
    filename = f"{treatment}_gsea_go_enrichment_results_{PERMUTATIONS}.csv"
    return BASE / folder / filename


def load_csv(norm: str, config: str, treatment: str) -> pd.DataFrame | None:
    path = build_csv_path(norm, config, treatment)
    print(f"    [path] {norm:12s} -> {path}")
    if not path.exists():
        print("    [skip] file does not exist")
        return None
    df = pd.read_csv(path)
    missing = [c for c in STAT_COLS + ["go_id", "Term"] if c not in df.columns]
    if missing:
        print(f"    [warn] missing columns {missing}, skipping")
        return None
    print(f"    [ok]   loaded {len(df)} GO terms")
    return df


def get_term_row_with_rank(df: pd.DataFrame, go_id: str) -> dict | None:
    """
    Find the expected GO term's row in a GSEA results dataframe, AND its
    rank by FDR q-val among all tested terms (rank 1 = most significant).

    Mirrors the proxy-fallback logic used in get_spider_plots /
    _get_term_row (diff_and_GSEA_pipeline.py): if the exact root term has
    no row (common for broad GO terms with no direct gene annotations),
    falls back to nothing here deliberately — pathway RECOVERY should be
    evaluated against the literal expected term, not a descendant proxy,
    since "did we recover the descendant instead of the root" is a
    different (weaker) claim than "did we recover the expected term."
    Treat a missing exact match as a non-recovery (rank = NaN), and print
    why, so it's visible in the console output rather than silently
    becoming an empty cell in the figure.
    """
    ranked = df.sort_values("FDR q-val", ascending=True).reset_index(drop=True)
    match = ranked[ranked["go_id"] == go_id]

    if match.empty:
        print(f"      [recovery] go_id {go_id} not found among {len(df)} tested terms "
              f"(likely no direct gene annotations at this root — see proxy-fallback "
              f"note in pr_rank_gene_enrich.py if you want descendant-level recovery instead)")
        return None

    row = match.iloc[0]
    rank = int(match.index[0]) + 1  # 1-indexed rank by FDR q-val ascending
    print(f"      [recovery] go_id {go_id} found at rank {rank}/{len(ranked)} "
          f"(FDR q-val={row['FDR q-val']:.4g}, NES={row['NES']:.4g})")

    return {
        "go_id": go_id,
        "Term": row["Term"],
        "ES": row["ES"],
        "NES": row["NES"],
        "NOM p-val": row["NOM p-val"],
        "FDR q-val": row["FDR q-val"],
        "FWER p-val": row["FWER p-val"],
        "rank": rank,
        "n_terms_tested": len(ranked),
    }


def collect_recovery_table(config: str, fdr_threshold: float) -> pd.DataFrame:
    """
    For every stress treatment and every normalization method, load the
    GSEA results, locate the expected GO term, and record its stats +
    rank + hit/miss flag. Returns one long-format dataframe.
    """
    rows = []

    for treatment in STRESS_TREATMENTS:
        t = treatment.value
        print(f"\n--- {t} ---")

        if t not in STRESS_GO_ROOTS:
            print(f"    [warn] no GO root mapping defined for '{t}' — skipping")
            continue

        go_id, go_name = STRESS_GO_ROOTS[t]
        # bad_flag = " [** KNOWN-BAD MAPPING — see KNOWN_BAD_GO_MAPPINGS **]" if t in KNOWN_BAD_GO_MAPPINGS else ""
        # print(f"    Expected pathway: {go_id} ({go_name}){bad_flag}")

        for norm in NORM_METHODS:
            df = load_csv(norm, config, t)
            if df is None:
                continue

            result = get_term_row_with_rank(df, go_id)
            if result is None:
                rows.append({
                    "treatment": t,
                    "norm_method": norm,
                    "go_id": go_id,
                    "go_name": go_name,
                    "ES": np.nan, "NES": np.nan, "NOM p-val": np.nan,
                    "FDR q-val": np.nan, "FWER p-val": np.nan,
                    "rank": np.nan, "n_terms_tested": np.nan,
                    "hit": False,
                    # "known_bad_mapping": t in KNOWN_BAD_GO_MAPPINGS,
                })
                continue

            hit = bool(result["FDR q-val"] < fdr_threshold) and bool(result["NES"] > 0)
            print(f"      [recovery] hit={hit} (threshold: FDR<{fdr_threshold} AND NES>0)")

            rows.append({
                "treatment": t,
                "norm_method": norm,
                "go_id": go_id,
                "go_name": go_name,
                **{k: v for k, v in result.items() if k not in ("go_id", "Term")},
                "hit": hit,
                # "known_bad_mapping": t in KNOWN_BAD_GO_MAPPINGS,
            })

    return pd.DataFrame(rows)


def plot_recovery_heatmap(recovery_df: pd.DataFrame, fdr_threshold: float, out_dir: Path) -> Path:
    """
    Treatment (rows) x normalization method (columns), colored by
    -log10(FDR q-val) of the expected term. Cells with a confirmed hit
    (FDR < threshold AND NES > 0) get a star marker.
    """
    treatments = [t.value for t in STRESS_TREATMENTS if t.value in recovery_df["treatment"].unique()]
    norms = NORM_METHODS

    score_matrix = np.full((len(treatments), len(norms)), np.nan)
    hit_matrix = np.zeros((len(treatments), len(norms)), dtype=bool)

    for i, t in enumerate(treatments):
        for j, n in enumerate(norms):
            sub = recovery_df[(recovery_df["treatment"] == t) & (recovery_df["norm_method"] == n)]
            if sub.empty or sub["FDR q-val"].isna().all():
                continue
            fdr = sub["FDR q-val"].iloc[0]
            score_matrix[i, j] = -np.log10(max(fdr, 1e-300))
            hit_matrix[i, j] = bool(sub["hit"].iloc[0])

    fig, ax = plt.subplots(figsize=(7.5, max(5, len(treatments) * 0.45)))
    im = ax.imshow(score_matrix, aspect="auto", cmap="viridis")
    ax.set_xticks(range(len(norms)))
    ax.set_xticklabels(norms, rotation=30, ha="right")
    ax.set_yticks(range(len(treatments)))

    ytick_labels = []
    # for t in treatments:
    #     # flag = " *" if t in KNOWN_BAD_GO_MAPPINGS else ""
    #     ytick_labels.append(f"{t}{flag}")
    ax.set_yticklabels(ytick_labels)

    for i in range(len(treatments)):
        for j in range(len(norms)):
            if hit_matrix[i, j]:
                ax.text(j, i, "\u2605", ha="center", va="center", color="white", fontsize=12)

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("-log10(FDR q-val) of expected pathway")

    ax.set_title(
        "Pathway recovery per treatment\n"
        f"(\u2605 hit: FDR<{fdr_threshold} & NES>0  |  * label = known-bad mapping)",
        fontsize=10,
    )
    plt.tight_layout()

    out_path = out_dir / "pathway_recovery_heatmap.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


def plot_rank_strip(recovery_df: pd.DataFrame, out_dir: Path) -> Path:
    """
    Per treatment, where does the expected term rank among all tested
    GO terms, for each normalization method? Lower rank = more strongly
    recovered. One strip of points per treatment.
    """
    treatments = [t.value for t in STRESS_TREATMENTS if t.value in recovery_df["treatment"].unique()]
    norm_colors = {"filter_norm": "tab:blue", "combat_norm": "tab:orange", "rankin": "tab:green"}

    fig, ax = plt.subplots(figsize=(max(8, len(treatments) * 0.7), 5))

    for x, t in enumerate(treatments):
        for norm, color in norm_colors.items():
            sub = recovery_df[(recovery_df["treatment"] == t) & (recovery_df["norm_method"] == norm)]
            if sub.empty or sub["rank"].isna().all():
                continue
            rank = sub["rank"].iloc[0]
            ax.scatter(x, rank, color=color, label=norm, zorder=3)

    ax.set_xticks(range(len(treatments)))
    ax.set_xticklabels(treatments, rotation=45, ha="right")
    ax.set_ylabel("Rank of expected term among tested GO terms\n(1 = most significant)")
    ax.invert_yaxis()  # rank 1 at top — visually "better" is higher
    ax.set_title("Where does the expected stress-response pathway rank?\n(lower = stronger recovery)")

    handles, labels = ax.get_legend_handles_labels()
    unique = dict(zip(labels, handles))
    ax.legend(unique.values(), unique.keys(), title="normalization")

    plt.tight_layout()
    out_path = out_dir / "pathway_recovery_rank_strip.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


def main():
    parser = argparse.ArgumentParser(
        description="Pathway recovery analysis: does GSEA recover the literature-"
                     "expected stress-response GO term for each treatment, and how "
                     "does that compare across normalization methods?"
    )
    parser.add_argument("--config", default="All_tissues_full_mixed_min_group_0",
                         help="Folder suffix for the run config to evaluate")
    parser.add_argument("--fdr", type=float, default=0.25,
                         help="FDR q-val threshold for calling recovery a 'hit' (default: 0.25, "
                              "matching standard GSEA convention)")
    parser.add_argument("--out-dir", default="./pathway_recovery_plots",
                         help="Directory to write plots and the combined CSV into")
    parser.add_argument("--out-csv", default="pathway_recovery_table.csv",
                         help="Filename (within --out-dir) for the combined long-format CSV")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("PATHWAY RECOVERY ANALYSIS")
    print("=" * 70)
    print(f"Config:        {args.config}")
    print(f"FDR threshold: {args.fdr}")
    print(f"GSEA base dir: {BASE}")
    print(f"Treatments:    {[t.value for t in STRESS_TREATMENTS]}")
    # if KNOWN_BAD_GO_MAPPINGS:
    #     print(f"\n[!] WARNING: known-bad GO mappings present for: {sorted(KNOWN_BAD_GO_MAPPINGS)}")
    #     print("    These results will be computed and plotted, but flagged — do not")
    #     print("    report them as biological validation without fixing the mapping first.")
    print("=" * 70)

    recovery_df = collect_recovery_table(args.config, args.fdr)

    if recovery_df.empty:
        print("\nNo recovery data collected — nothing to plot. Check paths/config above.")
        return

    csv_path = out_dir / args.out_csv
    recovery_df.to_csv(csv_path, index=False)
    print(f"\nWrote combined recovery table ({len(recovery_df)} rows) -> {csv_path}")

    print("\n--- Recovery summary (hits by normalization method) ---")
    summary = recovery_df.groupby("norm_method")["hit"].sum()
    totals = recovery_df.groupby("norm_method")["hit"].count()
    for norm in NORM_METHODS:
        if norm in summary.index:
            print(f"  {norm:12s}: {int(summary[norm])} / {int(totals[norm])} treatments recovered")

    print("\n--- Generating figures ---")
    heatmap_path = plot_recovery_heatmap(recovery_df, args.fdr, out_dir)
    print(f"  wrote {heatmap_path}")
    rank_path = plot_rank_strip(recovery_df, out_dir)
    print(f"  wrote {rank_path}")

    print(f"\nDone. Outputs in {out_dir}")


if __name__ == "__main__":
    main()
