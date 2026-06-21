"""
compare_gsea_overlap.py
------------------------
Compares GSEA GO-term enrichment results across normalization methods
(filter_norm, combat_norm, rankin) for a given treatment, to quantify
how consistent "significant" terms are between methods.

Usage:
    python compare_gsea_overlap.py --treatment Cold
    python compare_gsea_overlap.py --treatment Cold --fdr 0.1
    python compare_gsea_overlap.py --treatment Cold --config leaf_full_mixed_min_group_0
"""

import argparse
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

# Folder suffix matching your run config, e.g. "All_tissues_full_mixed_min_group_0"
DEFAULT_CONFIG = "All_tissues_full_mixed_min_group_0"


def build_csv_path(norm: str, config: str, treatment: str) -> Path:
    """
    Reconstruct the GSEA results CSV path for a given normalization method.

    Matches the pattern:
      {BASE}/GSEA_enrichment_{EXPERIMENT_NAME}_{norm}_{config}/
          {treatment}_gsea_go_enrichment_results_{PERMUTATIONS}.csv
    """
    folder = f"GSEA_enrichment_{EXPERIMENT_NAME}_{norm}_{config}"
    filename = f"{treatment}_gsea_go_enrichment_results_{PERMUTATIONS}.csv"
    return BASE / folder / filename


def load_results(treatment: str, config: str) -> dict[str, pd.DataFrame]:
    """Load available GSEA result CSVs per normalization method."""
    dfs = {}
    for norm in NORM_METHODS:
        path = build_csv_path(norm, config, treatment)
        if not path.exists():
            print(f"  [skip] {norm}: file not found at {path}")
            continue
        dfs[norm] = pd.read_csv(path)
        print(f"  [ok]   {norm}: loaded {len(dfs[norm])} terms from {path}")
    return dfs


def compare_significant_sets(
    dfs: dict[str, pd.DataFrame], fdr_threshold: float
) -> dict[str, set[str]]:
    """Build per-method sets of go_id passing the FDR threshold."""
    sig_sets = {}
    for name, df in dfs.items():
        if "go_id" not in df.columns or "FDR q-val" not in df.columns:
            print(f"  [warn] {name}: missing expected columns, skipping")
            continue
        sig = set(df.loc[df["FDR q-val"] < fdr_threshold, "go_id"])
        sig_sets[name] = sig
        print(f"  {name}: {len(sig)} significant terms (FDR < {fdr_threshold})")
    return sig_sets


def print_pairwise_overlap(sig_sets: dict[str, set[str]]) -> None:
    """Print pairwise overlap counts and Jaccard similarity between methods."""
    names = list(sig_sets.keys())
    if len(names) < 2:
        print("\nNeed at least two successfully loaded methods to compare overlap.")
        return

    print("\nPairwise overlap of significant terms:")
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            set_a, set_b = sig_sets[a], sig_sets[b]
            overlap = set_a & set_b
            union = set_a | set_b
            jaccard = len(overlap) / len(union) if union else float("nan")
            print(
                f"  {a:10s} ∩ {b:10s}: {len(overlap):3d} shared "
                f"(out of {len(set_a)}, {len(set_b)})  "
                f"Jaccard={jaccard:.2%}"
            )

    if len(names) >= 2:
        common_to_all = set.intersection(*sig_sets.values())
        print(f"\nTerms significant in ALL {len(names)} methods: {len(common_to_all)}")
        if common_to_all:
            for go_id in sorted(common_to_all):
                print(f"    {go_id}")


def main():
    parser = argparse.ArgumentParser(
        description="Compare GSEA significant-term overlap across normalization methods."
    )
    parser.add_argument(
        "--treatment", required=True,
        help="Treatment name as used in filenames, e.g. Cold, Salinity, 'Low Light'",
    )
    parser.add_argument(
        "--config", default=DEFAULT_CONFIG,
        help=f"Folder suffix for the run config (default: {DEFAULT_CONFIG})",
    )
    parser.add_argument(
        "--fdr", type=float, default=0.25,
        help="FDR q-val threshold for 'significant' (default: 0.25, matching GSEA convention)",
    )
    args = parser.parse_args()

    print(f"Comparing GSEA results for treatment='{args.treatment}', config='{args.config}'\n")

    print("Loading result files:")
    dfs = load_results(args.treatment, args.config)

    if len(dfs) < 2:
        print("\nFewer than 2 normalization methods loaded — cannot compare. Exiting.")
        return

    print(f"\nApplying FDR < {args.fdr} threshold:")
    sig_sets = compare_significant_sets(dfs, args.fdr)

    print_pairwise_overlap(sig_sets)


if __name__ == "__main__":
    main()