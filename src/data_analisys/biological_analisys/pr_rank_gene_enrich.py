"""
pr_rank_gene_enrich.py
----------------------
Provides two public functions:

  get_go_data(...)		  – Parse GO OBO + TAIR GAF annotation files,
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

    # --- Apply filters and propagate lineages ---
    print("Filtering and propagating annotations...")
    filtered: dict[str, set[str]] = {}
    for gene_id, go_ids in geneid2gos.items():
        kept = set()
        
        # 1. Expand direct annotations to include all parents (True Path Rule)
        expanded_go_ids = set()
        for go_id in go_ids:
            if go_id in obodag:
                expanded_go_ids.add(go_id)
                expanded_go_ids.update(obodag[go_id].get_all_parents())
                
        # 2. Evaluate the expanded set against your filters
        for go_id in expanded_go_ids:
            passes_ns = not namespaces or obodag[go_id].namespace in namespaces
            passes_stress = not all_stress_go_ids or go_id in all_stress_go_ids
            
            if passes_ns and passes_stress:
                kept.add(go_id)
                
        if kept:
            filtered[gene_id] = kept

    print(f"  {len(filtered)} genes retained after filtering and propagation.")
    return obodag, filtered

def get_go_data_old(
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
# 3. Shared GO Gene-Set Construction
# =============================================================================

def build_go_gene_sets(
	obodag: GODag,
	geneid2gos: dict[str, set[str]],
	keys: list | None = None,
) -> dict[str, list[str]]:
	"""
	Build a {GO:XXXXXXX (term name): [gene, ...]} dict by inverting
	`geneid2gos`. Shared by both GSEA prerank and ORA, so both methods
	are tested against exactly the same gene-set universe.

	Parameters
	----------
	obodag : GODag
		Parsed GO DAG (from :func:`get_go_data`).
	geneid2gos : dict[str, set[str]]
		Filtered gene → GO mapping (from :func:`get_go_data`).
	keys : list or None
		If given, restrict gene sets to those whose GO ID starts with one
		of these keys (legacy filter; set to None to keep all).

	Returns
	-------
	dict[str, list[str]]
	"""
	go_gene_sets: dict[str, list[str]] = {}
	for gene, go_ids in geneid2gos.items():
		for go_id in go_ids:
			label = f"{go_id} ({obodag[go_id].name if go_id in obodag else 'Unknown'})"
			go_gene_sets.setdefault(label, []).append(gene.upper())

	if keys is not None:
		keys_set = set(keys)
		go_gene_sets = {k: v for k, v in go_gene_sets.items() if k.split(" ")[0] in keys_set}

	return go_gene_sets


# =============================================================================
# 4. GSEA Prerank
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
	go_gene_sets = build_go_gene_sets(obodag, geneid2gos, keys=keys)
	print(f"  {len(go_gene_sets)} GO term gene sets prepared.")
	# print(f"{go_gene_sets}")

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


# =============================================================================
# 5. Overrepresentation Analysis (ORA) — complementary to GSEA prerank
# =============================================================================

def perform_ora_enrichment(
	diff_results: pd.DataFrame,
	gene_col: str,
	obodag: GODag,
	geneid2gos: dict[str, set[str]],
	keys: list | None,
	background_genes: list[str],
	adj_p_threshold: float = 0.05,
	logfc_threshold: float = 1.0,
	out_path: str | None = None,
) -> pd.DataFrame:
	"""
	Run a hypergeometric overrepresentation test (ORA) on a thresholded
	significant-gene list, using gseapy.enrich. Complementary to
	:func:`perform_gsea_enrichment`: ORA tests a hard-thresholded gene
	list rather than a continuous ranking, so it is far less sensitive
	to how the ranking metric was constructed (ties, near-zero p-values,
	clipping choices, etc.). Useful for checking whether GSEA results
	that look unstable across normalization methods reflect genuinely
	unstable differential expression calls, or just instability in the
	continuous ranking metric itself.

	Parameters
	----------
	diff_results : pd.DataFrame
		Differential expression results. Must contain `gene_col`,
		"adj.P.Val", and "logFC" columns.
	gene_col : str
		Column name for gene identifiers (AGI codes).
	obodag : GODag
		Parsed GO DAG (from :func:`get_go_data`).
	geneid2gos : dict[str, set[str]]
		Filtered gene → GO mapping (from :func:`get_go_data`).
	keys : list or None
		Optional legacy GO-ID-prefix filter, same as in
		:func:`perform_gsea_enrichment`.
	background_genes : list[str]
		The full gene universe actually tested in this comparison
		(e.g. the index of the ranked gene list used for GSEA on the
		same contrast). Using the wrong background — e.g. the entire
		TAIR annotation set instead of just the genes that passed
		expression filtering for this run — will bias ORA p-values,
		since gseapy.enrich treats `background` as "all genes that
		could have been selected."
	adj_p_threshold : float
		Adjusted p-value cutoff for calling a gene significant.
	logfc_threshold : float
		Absolute logFC cutoff for calling a gene significant.
	out_path : str or None
		If given, passed to gseapy.enrich as `outdir` so it also
		writes its own output files; if None, results are only
		returned as a DataFrame.

	Returns
	-------
	pd.DataFrame
		gseapy ORA results (Term, go_id, Overlap, P-value, Adjusted P-value,
		Genes, ...), sorted by ascending Adjusted P-value. Empty DataFrame
		if no genes pass the significance threshold.
	"""
	print("Building GO gene sets (ORA)...")
	go_gene_sets = build_go_gene_sets(obodag, geneid2gos, keys=keys)
	print(f"  {len(go_gene_sets)} GO term gene sets prepared.")

	sig_mask = (
		(diff_results["adj.P.Val"] < adj_p_threshold)
		& (diff_results["logFC"].abs() > logfc_threshold)
	)
	sig_genes = diff_results.loc[sig_mask, gene_col].str.upper().unique().tolist()
	print(
		f"  {len(sig_genes)} genes pass adj.P.Val < {adj_p_threshold} "
		f"& |logFC| > {logfc_threshold}."
	)

	if not sig_genes:
		print("Error: no genes pass the significance threshold for ORA.")
		return pd.DataFrame(
			columns=["Term", "go_id", "Overlap", "P-value", "Adjusted P-value", "Genes"]
		)

	fallback_cols = ["Term", "go_id", "Overlap", "P-value", "Adjusted P-value", "Genes"]

	if not sig_genes:
		print("Error: no genes pass the significance threshold for ORA.")
		return pd.DataFrame(columns=fallback_cols)

	background = [g.upper() for g in background_genes]
	print(f"  Background gene universe: {len(background)} genes.")

	print("Running ORA (gseapy.enrich)...")
	try:
		enr = gseapy.enrich(
			gene_list=sig_genes,
			gene_sets=go_gene_sets,
			background=background,
			outdir=out_path,
		)
	except Exception as e:
		print(f"  ! Runtime exception in gseapy.enrich: {e}")
		return pd.DataFrame(columns=fallback_cols)
	if not isinstance(enr.results, pd.DataFrame) or enr.results.empty:
		print("  ! Warning: No valid enrichment results returned (0 overlapping hits).")
		return pd.DataFrame(columns=fallback_cols)
	results_df = enr.results.copy()
	results_df["go_id"] = results_df["Term"].str.split(" ").str[0]

	sort_col = "Adjusted P-value" if "Adjusted P-value" in results_df.columns else "P-value"
	return results_df.sort_values(by=sort_col, ascending=True)


# =============================================================================
# 6. ORA term lookup — shared primitive for pathway-recovery checks
# =============================================================================

def get_ora_term_row(
	ora_results: pd.DataFrame,
	go_id: str,
) -> dict | None:
	"""
	Find one GO term's row in an ORA results DataFrame (from
	:func:`perform_ora_enrichment`), along with its rank by Adjusted
	P-value among every term ORA actually tested (rank 1 = most
	significant).

	This is the single-call primitive behind pathway-recovery checks: do we
	look for the literal expected term here, by design, rather than falling
	back to a GO descendant — "recovered the literal expected term" and
	"recovered a related descendant" are different-strength claims, and
	conflating them silently would overstate recovery. Callers that want
	descendant-level fallback (e.g. matching the spider-plot behavior in
	diff_and_GSEA_pipeline.py's `_get_term_row`) should implement that
	separately using `obodag[go_id].get_all_children()`, since it needs
	the GODag object that this function deliberately does not require.

	No assumption is made about a treatment/stress label — this function
	is treatment-agnostic; callers building treatment-aware tables (e.g.
	a recovery table across many treatments and normalization methods)
	should loop over their own treatment -> go_id mapping and call this
	once per (treatment, normalization method) pair, exactly as
	perform_gsea_enrichment's callers loop per stress treatment.

	Parameters
	----------
	ora_results : pd.DataFrame
		Output of :func:`perform_ora_enrichment`. May be empty (no genes
		passed the significance threshold) — handled gracefully, returns
		None rather than raising.
	go_id : str
		The literal GO ID to look up (e.g. "GO:0009409").

	Returns
	-------
	dict or None
		None if `ora_results` is empty or `go_id` was not among the terms
		ORA actually tested. Otherwise a dict with keys:
		go_id, Term, Overlap, P-value, Adjusted P-value, rank,
		n_terms_tested.
	"""
	if ora_results is None or ora_results.empty:
		return None

	sort_col = "Adjusted P-value" if "Adjusted P-value" in ora_results.columns else "P-value"
	ranked = ora_results.sort_values(sort_col, ascending=True).reset_index(drop=True)
	match = ranked[ranked["go_id"] == go_id]

	if match.empty:
		return None

	row = match.iloc[0]
	rank = int(match.index[0]) + 1

	return {
		"go_id": go_id,
		"Term": row["Term"],
		"Overlap": row["Overlap"],
		"P-value": row["P-value"],
		"Adjusted P-value": row.get("Adjusted P-value", float("nan")),
		"rank": rank,
		"n_terms_tested": len(ranked),
	}


def count_ora_significant_genes(ora_results: pd.DataFrame) -> int:
	"""
	Lower-bound count of genes that contributed to at least one ORA hit,
	by taking the union of every term's "Genes" field (gseapy's
	semicolon-separated convention).

	This UNDERCOUNTS the true significant-gene-list size whenever a
	significant gene didn't land in any tested GO term at all (e.g. an
	unannotated gene, or a gene whose only annotations were filtered out
	upstream by `stress_root_go_ids` in :func:`get_go_data`) — there is no
	way to recover the true count from `ora_results` alone, since
	:func:`perform_ora_enrichment` does not currently return the original
	`sig_genes` list alongside its results. Treat this as a diagnostic
	lower bound, useful for distinguishing "ORA found nothing because
	there's no signal" from "ORA found nothing because almost no genes
	cleared the hard threshold at this sample size" — not as an exact
	significant-gene count.
	"""
	if ora_results is None or ora_results.empty or "Genes" not in ora_results.columns:
		return 0
	all_genes: set[str] = set()
	for genes_str in ora_results["Genes"].dropna():
		all_genes.update(g.strip() for g in str(genes_str).split(";") if g.strip())
	return len(all_genes)