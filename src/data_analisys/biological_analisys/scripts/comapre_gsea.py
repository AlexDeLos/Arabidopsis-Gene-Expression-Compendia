"""
comapre_gsea.py
----------------
Compares GSEA GO-term enrichment results across normalization methods
(filter_norm, combat_norm, rankin) for ALL stress treatments, surfacing
both:

  1. The raw stat values (ES, NES, NOM p-val, FDR q-val, FWER p-val) for
     each GO term, side by side across normalization methods — so you can
     see directly how normalization shifts the numbers for the same term.
  2. Significant-term overlap (Jaccard) per treatment, as before.

Outputs one combined long-format CSV (one row per treatment × go_id × norm
method) plus a per-treatment console summary.

Usage:
    python comapre_gsea.py
    python comapre_gsea.py --config All_tissues_full_mixed_min_group_0
    python comapre_gsea.py --fdr 0.1 --out-csv my_comparison.csv
"""

import argparse
from enum import Enum
from pathlib import Path

import pandas as pd

# =============================================================================
# Paths — adjust BASE if your output root ever moves
# =============================================================================

BASE = Path(
    "/tudelft.net/staff-umbrella/GeneExpressionStorage/outputs/5.0/GSEA_enrichment_results"
)

EXPERIMENT_NAME = "5.0"
PERMUTATIONS = 1000

NORM_METHODS = ["filter_norm", "combat_norm", "rankin"]

DEFAULT_CONFIG = "All_tissues_full_mixed_min_group_0"

STAT_COLS = ["ES", "NES", "NOM p-val", "FDR q-val", "FWER p-val"]


# =============================================================================
# Treatments — excludes non-stress / structural labels (Control, Other,
# unknown, unspecified) since no GSEA contrast was ever run "treatment vs
# control" for those — they ARE the control / catch-all side of contrasts.
# =============================================================================

class TreatmentEnum(str, Enum):
    DROUGHT = "Drought"
    FLOOD = "Flood"
    DEHYDRATION = "Dehydration"
    SALINITY = "Salinity"
    HEAT = "Heat"
    COLD = "Cold"
    CHEMICAL = "Chemical"
    BIOTIC = "Biotic"
    ABIOTIC = "Abiotic"
    LOW_LIGHT = "Low Light"
    HIGH_LIGHT = "High Light"
    OTHER_LIGHT = "Other Light"
    CUT = "Cut"
    NUTRIENT_DEFICIENCY = "Nutrient Deficiency"
    OTHER = "Other"
    CONTROL = "Control"
    UNKNOWN = "unknown"
    UNSPECIFIED = "unspecified"


NON_STRESS_LABELS = {
    TreatmentEnum.OTHER,
    TreatmentEnum.CONTROL,
    TreatmentEnum.UNKNOWN,
    TreatmentEnum.UNSPECIFIED,
}

STRESS_TREATMENTS = [t for t in TreatmentEnum if t not in NON_STRESS_LABELS]


def build_csv_path(norm: str, config: str, treatment: str) -> Path:
    """Reconstruct the GSEA results CSV path for a given normalization method."""
    folder = f"GSEA_enrichment_{EXPERIMENT_NAME}_{norm}_{config}"
    filename = f"{treatment}_gsea_go_enrichment_results_{PERMUTATIONS}.csv"
    return BASE / folder / filename


def load_results(treatment: str, config: str) -> dict[str, pd.DataFrame]:
    """Load available GSEA result CSVs per normalization method for one treatment."""
    dfs = {}
    for norm in NORM_METHODS:
        path = build_csv_path(norm, config, treatment)
        if not path.exists():
            print(f"    [skip] {norm}: not found")
            continue
        df = pd.read_csv(path)
        missing_cols = [c for c in STAT_COLS + ["go_id", "Term"] if c not in df.columns]
        if missing_cols:
            print(f"    [warn] {norm}: missing columns {missing_cols}, skipping")
            continue
        dfs[norm] = df
        print(f"    [ok]   {norm}: {len(df)} terms")
    return dfs


def build_stat_comparison(treatment: str, dfs: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """
    Long-format table: one row per (go_id, norm_method), with the five
    stat columns, for a single treatment. Makes it easy to pivot/filter
    later (e.g. df.pivot(index="go_id", columns="norm_method", values="NES")).
    """
    rows = []
    for norm, df in dfs.items():
        sub = df[["go_id", "Term"] + STAT_COLS].copy()
        sub["treatment"] = treatment
        sub["norm_method"] = norm
        rows.append(sub)
    if not rows:
        return pd.DataFrame()
    return pd.concat(rows, ignore_index=True)


def summarize_overlap(treatment: str, dfs: dict[str, pd.DataFrame], fdr_threshold: float) -> None:
    """Print Jaccard overlap of significant terms across normalization methods."""
    sig_sets = {
        norm: set(df.loc[df["FDR q-val"] < fdr_threshold, "go_id"])
        for norm, df in dfs.items()
    }
    names = list(sig_sets.keys())
    if len(names) < 2:
        print("    (fewer than 2 methods loaded — skipping overlap)")
        return

    for i, a in enumerate(names):
        for b in names[i + 1:]:
            set_a, set_b = sig_sets[a], sig_sets[b]
            overlap = set_a & set_b
            union = set_a | set_b
            jaccard = len(overlap) / len(union) if union else float("nan")
            print(
                f"    {a:11s} ∩ {b:11s}: {len(overlap):3d} shared "
                f"(of {len(set_a)}, {len(set_b)})  Jaccard={jaccard:.1%}"
            )

    common_to_all = set.intersection(*sig_sets.values()) if len(names) >= 2 else set()
    print(f"    Significant in ALL {len(names)} methods: {len(common_to_all)} "
          f"{sorted(common_to_all) if common_to_all else ''}")


def main():
    parser = argparse.ArgumentParser(
        description="Compare GSEA stat values + significant-term overlap "
                     "across normalization methods, for all stress treatments."
    )
    parser.add_argument(
        "--config", default=DEFAULT_CONFIG,
        help=f"Folder suffix for the run config (default: {DEFAULT_CONFIG})",
    )
    parser.add_argument(
        "--fdr", type=float, default=0.25,
        help="FDR q-val threshold for 'significant' (default: 0.25)",
    )
    parser.add_argument(
        "--out-csv", default="gsea_norm_comparison.csv",
        help="Path to write the combined long-format stat comparison CSV",
    )
    args = parser.parse_args()

    print(f"Config: {args.config}")
    print(f"Treatments to compare: {[t.value for t in STRESS_TREATMENTS]}\n")

    all_comparisons = []

    for treatment in STRESS_TREATMENTS:
        print(f"--- {treatment.value} ---")
        dfs = load_results(treatment.value, args.config)

        if len(dfs) < 1:
            print("    No files found for any normalization method. Skipping.\n")
            continue

        comparison = build_stat_comparison(treatment.value, dfs)
        if not comparison.empty:
            all_comparisons.append(comparison)

        if len(dfs) >= 2:
            print(f"  Applying FDR < {args.fdr} threshold:")
            summarize_overlap(treatment.value, dfs, args.fdr)
        else:
            print("    Only one normalization method available — no overlap to compute.")
        print()

    if not all_comparisons:
        print("No comparable results found for any treatment. Nothing written.")
        return

    combined = pd.concat(all_comparisons, ignore_index=True)
    combined.to_csv(args.out_csv, index=False)
    print(f"Wrote combined stat comparison ({len(combined)} rows) → {args.out_csv}")
    print("\nColumns: treatment, go_id, Term, norm_method, "
          + ", ".join(STAT_COLS))
    print("\nTip: pivot for a specific stat to see normalization side by side, e.g.:")
    print("  df = pd.read_csv(out_csv)")
    print("  df[df['treatment']=='Cold'].pivot(index='go_id', columns='norm_method', values='NES')")


if __name__ == "__main__":
    main()