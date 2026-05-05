"""
pr_rank_gene_enrich.py
----------------------
Provides two public functions:

  get_go_data(...)          – Parse GO OBO + TAIR GAF annotation files,
                              optionally filtered to stress-related GO terms.

  perform_gsea_enrichment(...)  – Run GSEApy prerank on a ranked gene list
                                  and return a tidy results DataFrame.
"""

import gzip
import os
import re

import Bio.UniProt.GOA as GOA
import gseapy
import pandas as pd
from goatools.base import download_go_basic_obo
from goatools.obo_parser import GODag


# =============================================================================
# 1. GAF Parsing
# =============================================================================

def parse_gaf(annotation_file: str) -> dict[str, set[str]]:
    """
    Parse a TAIR GAF file and return a gene→GO mapping keyed by AGI codes.

    AGI codes (e.g. AT1G01010) are extracted from the 'Synonym' field of
    each GAF record.  Only canonical AGI patterns are accepted
    (AT[1-5MC]G followed by 5 digits).

    Parameters
    ----------
    annotation_file : str
        Path to a gzip-compressed GAF file (e.g. tair.gaf.gz).

    Returns
    -------
    dict[str, set[str]]
        Mapping of uppercase AGI code → set of GO term IDs.
    """
    geneid2gos: dict[str, set[str]] = {}
    agi_pattern = re.compile(r"AT[1-5MC]G\d{5}", re.IGNORECASE)

    with gzip.open(annotation_file, "rt") as fp:
        for entry in GOA.gafiterator(fp):
            go_id = entry["GO_ID"]
            synonyms = entry.get("Synonym", "")   # already a list in Biopython GAF

            for synonym in synonyms:
                match = agi_pattern.search(synonym)
                if match:
                    agi_id = match.group(0).upper()
                    geneid2gos.setdefault(agi_id, set()).add(go_id)

    return geneid2gos


# =============================================================================
# 2. GO Data Loading with Optional Filtering
# =============================================================================

def get_go_data(
    go_obo_file: str,
    annotation_file: str,
    namespaces: set[str] | None = None,
    stress_root_go_ids: set[str] | None = None,
) -> tuple[GODag, dict[str, set[str]]]:
    """
    Load and optionally filter GO ontology and TAIR annotation data.

    Downloads missing files automatically.  Filtering is cumulative: a GO
    term must pass *all* supplied filters to be retained.

    Parameters
    ----------
    go_obo_file : str
        Path to go-basic.obo.  Downloaded if absent.
    annotation_file : str
        Path to tair.gaf.gz.  Downloaded if absent.
    namespaces : set of str or None
        If given, keep only GO terms in these namespaces
        (e.g. {'biological_process'}).
    stress_root_go_ids : set of str or None
        If given, keep only GO terms that are descendants (or the root itself)
        of these GO IDs.  Used to restrict GSEA to stress-relevant terms.

    Returns
    -------
    obodag : GODag
    geneid2gos : dict[str, set[str]]
        Gene → filtered set of GO term IDs.
    """
    # --- OBO ---
    if not os.path.exists(go_obo_file):
        print(f"Downloading GO OBO to '{go_obo_file}'...")
        download_go_basic_obo(go_obo_file)

    print("Parsing GO OBO file...")
    obodag = GODag(go_obo_file)

    # --- Expand stress root IDs to full descendant set ---
    all_stress_go_ids: set[str] = set()
    if stress_root_go_ids:
        print(f"Expanding {len(stress_root_go_ids)} root stress GO terms to descendants...")
        for go_id in stress_root_go_ids:
            if go_id in obodag:
                go_term = obodag[go_id]
                all_stress_go_ids.add(go_id)
                all_stress_go_ids.update(go_term.get_all_children())
            else:
                print(f"  Warning: '{go_id}' not found in OBO DAG.")
        print(f"  Total stress-related GO terms (with descendants): {len(all_stress_go_ids)}")

    # --- GAF ---
    if not os.path.exists(annotation_file):
        print(f"Downloading TAIR annotation to '{annotation_file}'...")
        import requests
        r = requests.get("http://current.geneontology.org/annotations/tair.gaf.gz")
        with open(annotation_file, "wb") as f:
            f.write(r.content)

    print("Parsing TAIR annotation file...")
    geneid2gos = parse_gaf(annotation_file)

    # --- Apply filters ---
    if not namespaces and not all_stress_go_ids:
        print("No filters applied — returning full annotation set.")
        return obodag, geneid2gos

    print("Filtering annotations...")
    filtered: dict[str, set[str]] = {}
    for gene_id, go_ids in geneid2gos.items():
        kept = set()
        for go_id in go_ids:
            if go_id not in obodag:
                continue
            passes_ns = not namespaces or obodag[go_id].namespace in namespaces
            passes_stress = not all_stress_go_ids or go_id in all_stress_go_ids
            if passes_ns and passes_stress:
                kept.add(go_id)
        if kept:
            filtered[gene_id] = kept

    print(f"  {len(filtered)} genes retained after filtering.")
    return obodag, filtered


# =============================================================================
# 3. GSEA Prerank
# =============================================================================

def perform_gsea_enrichment(
    ranked_gene_df: pd.DataFrame,
    gene_col: str,
    rank_col: str,
    obodag: GODag,
    geneid2gos: dict[str, set[str]],
    keys: list | None,
    stress: str,
    out_path: str,
    permutations: int,
) -> pd.DataFrame:
    """
    Run GSEApy prerank and return a tidy results DataFrame.

    The GO gene sets are built by inverting `geneid2gos` and annotating each
    set name with its human-readable GO term name
    (format: ``"GO:XXXXXXX (term name)"``).

    After running, zero p-values are replaced with ``1 / permutations`` so
    that downstream -log10 transforms are finite.

    Parameters
    ----------
    ranked_gene_df : pd.DataFrame
        Must contain at least `gene_col` (AGI codes) and `rank_col`
        (numeric ranking metric, e.g. logFC × -log10(adj.P.Val)).
    gene_col : str
        Column name for gene identifiers.
    rank_col : str
        Column name for the ranking metric.  Genes are ranked descending.
    obodag : GODag
        Parsed GO DAG (from :func:`get_go_data`).
    geneid2gos : dict[str, set[str]]
        Filtered gene → GO mapping (from :func:`get_go_data`).
    keys : list or None
        If given, restrict gene sets to those whose GO ID starts with one of
        these keys (legacy filter; set to None to keep all).
    stress : str
        Stress label used to name the GSEApy output subdirectory.
    out_path : str
        Parent directory for GSEApy prerank output files.
    permutations : int
        Number of permutations for p-value estimation.

    Returns
    -------
    pd.DataFrame
        Sorted by descending ES.  Contains columns:
        Term, go_id, ES, NES, NOM p-val, FDR q-val, FWER p-val,
        enrichment_score, normalized_enrichment_score, leading_edge_genes.
        Zero p-values are set to ``1 / permutations``.
    """
    # ------------------------------------------------------------------
    # 1. Build GO gene-set dictionary  {GO:XXXXXXX (name): [gene, ...]}
    # ------------------------------------------------------------------
    print("Building GO gene sets...")
    go_gene_sets: dict[str, list[str]] = {}
    for gene, go_ids in geneid2gos.items():
        for go_id in go_ids:
            label = f"{go_id} ({obodag[go_id].name if go_id in obodag else 'Unknown'})"
            go_gene_sets.setdefault(label, []).append(gene.upper())

    # Optional key-based filter (legacy)
    if keys is not None:
        keys_set = set(keys)
        go_gene_sets = {k: v for k, v in go_gene_sets.items() if k.split(" ")[0] in keys_set}

    print(f"  {len(go_gene_sets)} GO term gene sets prepared.")

    # ------------------------------------------------------------------
    # 2. Prepare ranked gene list
    # ------------------------------------------------------------------
    print("Preparing ranked gene list...")
    rnk = ranked_gene_df[[gene_col, rank_col]].copy()
    rnk[gene_col] = rnk[gene_col].str.upper()
    rnk = rnk.loc[rnk.groupby(gene_col)[rank_col].idxmax()]
    rnk = rnk.set_index(gene_col).sort_values(by=rank_col, ascending=False)
    print(f"  {len(rnk)} unique genes in ranked list.")

    if rnk.empty:
        print("Error: Ranked gene list is empty after deduplication.")
        return pd.DataFrame()

    # ------------------------------------------------------------------
    # 3. Run GSEApy prerank
    # ------------------------------------------------------------------
    print(f"Running GSEA prerank ({permutations} permutations)...")
    pre_res = gseapy.prerank(
        rnk=rnk,
        gene_sets=go_gene_sets,
        permutation_num=permutations,
        outdir=f"{out_path}{stress}_gsea_prerank_results",
        ascending=False,
        verbose=True,
    )

    # ------------------------------------------------------------------
    # 4. Tidy output
    # ------------------------------------------------------------------
    print("Formatting results...")
    results_df = pre_res.res2d.copy()

    # Replace zero p-values with the minimum representable value
    p_floor = 1.0 / permutations
    for col in ["NOM p-val", "FDR q-val", "FWER p-val"]:
        results_df[col] = results_df[col].astype(float).replace(0.0, p_floor)

    # Friendly column aliases (keep originals too for compatibility)
    results_df.rename(
        columns={
            "es": "enrichment_score",
            "nes": "normalized_enrichment_score",
            "pval": "p_value",
            "fdr": "fdr_q_value",
            "genes": "leading_edge_genes",
        },
        inplace=True,
    )

    # Extract GO ID from the Term string
    results_df["go_id"] = results_df["Term"].str.split(" ").str[0]

    return results_df.sort_values(by="ES", ascending=False)
