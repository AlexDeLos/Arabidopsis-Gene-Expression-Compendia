import re
import gzip
import Bio.UniProt.GOA as GOA

import os
import pandas as pd
import gseapy # New library for GSEA
from goatools.base import download_go_basic_obo # For GO data
from goatools.obo_parser import GODag
from typing import Dict, Set, Optional

def parse_gaf(annotation_file: str) -> Dict[str, Set[str]]:
    """
    Parses a GAF file to extract GO terms, mapping them to AGI codes.
    AGI codes (e.g., AT1G01010) are extracted from the synonym column.
    """
    geneid2gos = {}
    # Regex to find AGI codes like AT1G01234 or ATMG01234
    agi_pattern = re.compile(r'AT[1-5MC]G\d{5}', re.IGNORECASE)

    with gzip.open(annotation_file, 'rt') as arab_gaf_fp:
        # The gafiterator skips comment lines automatically
        for entry in GOA.gafiterator(arab_gaf_fp):
            go_id = entry['GO_ID']
            
            # CORRECTED LINE: The key for synonyms is 'Synonym', not 'DB_Object_Synonym'
            synonyms_str = entry.get('Synonym', '')
            synonyms = synonyms_str#.split('|')
            
            for synonym in synonyms:
                match = agi_pattern.search(synonym)
                if match:
                    # Standardize to uppercase (e.g., at1g... -> AT1G...)
                    agi_id = match.group(0).upper()
                    if agi_id in geneid2gos:
                        geneid2gos[agi_id].add(go_id)
                    else:
                        geneid2gos[agi_id] = {go_id}
    return geneid2gos


def get_go_data(
    go_obo_file: str, 
    annotation_file: str,
    namespaces: Optional[Set[str]] = None,
    stress_root_go_ids: Optional[Set[str]] = None  # <-- NEW ARGUMENT
):
    """
    Downloads and parses GO and annotation files.
    This helper function separates the data loading from the analysis.
    
    Can optionally filter by namespace AND/OR a set of "root" stress GO IDs
    (which will be expanded to include all their descendants).
    """
    # Download GO OBO file if it doesn't exist
    if not os.path.exists(go_obo_file):
        print(f"GO OBO file not found. Downloading to '{go_obo_file}'...")
        download_go_basic_obo(go_obo_file)

    print("Parsing GO OBO file...")
    obodag = GODag(go_obo_file)

    # --- NEW: Calculate the complete set of stress-related GO IDs ---
    all_stress_go_ids = set()
    if stress_root_go_ids:
        print(f"Finding all descendant GO terms for {len(stress_root_go_ids)} root stress terms...")
        for go_id in stress_root_go_ids:
            if go_id in obodag:
                all_stress_go_ids.add(go_id) # Add the root term itself
                all_stress_go_ids.update(obodag[go_id].children)
            else:
                print(f"Warning: Root stress GO ID '{go_id}' not found in OBO DAG.")
        print(f"Total stress-related GO terms (including descendants): {len(all_stress_go_ids)}")
    # --- END NEW BLOCK ---

    # Download and parse Arabidopsis annotation file (GAF format)
    if not os.path.exists(annotation_file):
        print(f"Annotation file not found. Downloading to '{annotation_file}'...")
        import requests
        url = "http://current.geneontology.org/annotations/tair.gaf.gz"
        r = requests.get(url)
        with open(annotation_file, 'wb') as f:
            f.write(r.content)

    print("Parsing Arabidopsis annotation file...")
    # This returns a dictionary mapping gene ID -> set of GO IDs
    geneid2gos = parse_gaf(annotation_file)
    
    # Check if *any* filtering is needed
    if namespaces or all_stress_go_ids:
        if namespaces:
            print(f"Filtering annotations for namespaces: {namespaces}")
        if all_stress_go_ids:
             print(f"Filtering annotations for {len(all_stress_go_ids)} stress-related GO terms.")

        filtered_geneid2gos = {}
        for gene_id, go_ids in geneid2gos.items():
            
            filtered_set = set()
            for go_id in go_ids:
                # Basic check: must be in the DAG
                if go_id not in obodag:
                    continue
                    
                # Condition 1: Namespace filter
                # If namespaces are given, the GO term must be in them.
                passes_namespace = (
                    not namespaces or 
                    obodag[go_id].namespace in namespaces
                )
                
                # Condition 2: Stress filter
                # If stress_go_ids are given (i.e., all_stress_go_ids is not empty), 
                # the GO term must be in that set.
                passes_stress = (
                    not all_stress_go_ids or 
                    go_id in all_stress_go_ids
                )
                
                # Only keep the GO ID if it passes *both* filters
                if passes_namespace and passes_stress:
                    filtered_set.add(go_id)
            
            # If the gene still has annotations after filtering, add it to the results
            if filtered_set:
                filtered_geneid2gos[gene_id] = filtered_set
        
        print(f"Filtering complete. {len(filtered_geneid2gos)} genes have annotations matching all criteria.")
        return obodag, filtered_geneid2gos
        
    # If no filters were provided, return the full, unfiltered data
    print("No filters applied. Returning all parsed annotations.")
    return obodag, geneid2gos

def perform_gsea_enrichment(
    ranked_gene_df: pd.DataFrame,
    gene_col: str,
    rank_col: str,
    obodag: GODag,
    geneid2gos: Dict[str, Set[str]],
    keys: Optional[list],
    stress:str,
    out_path: str,
    permutations: int
) -> pd.DataFrame:
    """
    Performs Gene Set Enrichment Analysis (GSEA) for a ranked list of Arabidopsis genes.

    Args:
        ranked_gene_df (pd.DataFrame): DataFrame with at least two columns:
                                      one for gene identifiers (e.g., AGI codes) and
                                      one for the ranking metric (e.g., log2FoldChange, stat).
        gene_col (str): The name of the column containing gene identifiers.
        rank_col (str): The name of the column to rank genes by.
        obodag (GODag): A parsed GO DAG object from goatools.
        geneid2gos (Dict[str, Set[str]]): A dictionary mapping gene IDs to a set of GO IDs.

    Returns:
        pd.DataFrame: A DataFrame containing the GSEA results, sorted by significance.
    """
    if keys is not None:
        list_of_keys = list(obodag.keys())
        list_of_keys_start = map(lambda x: x.split(' ')[0],list_of_keys)
        for el,dic_el in zip(list_of_keys,list_of_keys_start):
            if dic_el not in keys:
                del obodag[el]
    # 1. Prepare the Gene Sets (GO terms) for GSEApy
    # GSEApy requires a dictionary mapping: {go_term: [gene1, gene2, ...]}
    print("Creating GO gene sets for GSEA...")
    go_gene_sets = {}
    # Invert the geneid2gos dictionary
    for gene, go_ids in geneid2gos.items():
        for go_id in go_ids:
            if go_id not in go_gene_sets:
                go_gene_sets[go_id] = []
            # Ensure gene IDs are uppercase to match common conventions
            go_gene_sets[go_id].append(gene.upper())

    # For better readability, map GO IDs to their names
    # THIS MAKES ME LOOSE SOME TREATMENTS
    go_gene_sets_named = {
        f"{go_id} ({obodag[go_id].name if go_id in obodag else 'Unknown'})": genes
        for go_id, genes in go_gene_sets.items()
    }
    print(f"Created {len(go_gene_sets_named)} GO term gene sets.")

    # 2. Prepare the ranked gene list for GSEApy
    print("Preparing ranked gene list...")
    rnk = ranked_gene_df[[gene_col, rank_col]].copy()
    rnk[gene_col] = rnk[gene_col].str.upper()
    # Remove any duplicates, keeping the one with the max rank value
    rnk = rnk.loc[rnk.groupby(gene_col)[rank_col].idxmax()]
    rnk = rnk.set_index(gene_col)
    rnk = rnk.sort_values(by=rank_col, ascending=False)
    print(f"Ranked list contains {len(rnk)} unique genes.")

    if rnk.empty:
        print("Error: The ranked gene list is empty after processing.")
        return pd.DataFrame()
    # 3. Run GSEA Prerank
    if keys is not None:
        list_of_keys = list(go_gene_sets_named.keys())
        list_of_keys_start = map(lambda x: x.split(' ')[0],list_of_keys)
        for el,dic_el in zip(list_of_keys,list_of_keys_start):
            if dic_el not in keys:
                del go_gene_sets_named[el]
    print("\nRunning GSEA Prerank analysis...")
    permutation_num = permutations
    pre_res = gseapy.prerank(
        rnk=rnk,
        gene_sets=go_gene_sets_named,
        # min_size=,              # Min size of a gene set to be considered
        # max_size=len(rnk),        # Max size of a gene set to be considered
        permutation_num=permutation_num,# 00,#2000,     # Number of permutations for p-value calculation
        outdir=f'{out_path}{stress}_gsea_prerank_results', # Directory to save plots and results
        ascending=False,
        verbose=True,
    )
    
    print("GSEA complete. Formatting results...")
    results_df = pre_res.res2d
    neg = 0
    results_df['NOM p-val'] = results_df['NOM p-val'].astype(float).replace(0, 1/permutation_num)
    results_df['FDR q-val'] = results_df['FDR q-val'].astype(float).replace(0, 1/permutation_num)
    results_df['FWER p-val'] = results_df['FWER p-val'].astype(float).replace(0, 1/permutation_num)
    # A standard FDR cutoff for GSEA is often more lenient (e.g., < 0.25)
    significant_results:pd.D = results_df #[results_df['FDR q-val'] < 0.05].copy()

    # Make the output more user-friendly
    significant_results.rename(columns={
        'es': 'enrichment_score',
        'nes': 'normalized_enrichment_score',
        'pval': 'p_value',
        'fdr': 'fdr_q_value',
        'genes': 'leading_edge_genes'
    }, inplace=True)
    
    # Extract GO ID and Name from the index (Term)
    significant_results['go_id'] = significant_results['Term'].str.split(' ').str[0]
    # significant_results = significant_results[[
    #     'go_id','Term', 'normalized_enrichment_score', 'p_value', 'fdr_q_value', 'leading_edge_genes'
    # ]]
    
    return significant_results.sort_values(by='ES', ascending=False)

if __name__ == "__main__":
    # 1. Define file paths
    constant_data_path = 'data/'
    out_path = 'GSEA_enrichment_corrected_0.sanity/'
    GO_OBO_FILE = f'{constant_data_path}/go-basic.obo'
    ANNOTATION_FILE = f'{constant_data_path}/goa_arabidopsis.gaf.gz'
    treatments = [
        # "Drought Stress",
        "Salinity Stress",
        "Heat Stress",
        "Cold Stress",
        # "Chemical Stress",
        # "Pathogen Attack",
        # "Low Light Stress",
        "High Light Stress",
        # "Red Light Stress",
        # "Other Light Stress"
        ]
    # 2. Load the necessary GO data first
    # This only needs to be done once

    go_terms = pd.read_excel(f'{constant_data_path}GO_terms_to_enrich_2.xlsx')

    for stress in treatments:
        try:
            obodag, geneid2gos = get_go_data(GO_OBO_FILE, ANNOTATION_FILE)

            diff_exp_results = pd.read_csv(f'diff_exp/0.sanity_2_way_norm/All-Tissues_{stress}_genes.csv')
            # diff_exp_results['rank'] = diff_exp_results.index
            # 3. Now run the GSEA enrichment
            # Let's say your differential expression results are in a DataFrame `diff_exp_results`
            # with columns 'ID' for gene names and 'logFC' for ranking.
            
            # Load your actual differential expression results
            # diff_exp_results = pd.read_excel('your_limma_output.xlsx') 
            keys = None#go_terms['GO Term Name']
            ids =  go_terms['GO Term ID']
            gsea_results_df = perform_gsea_enrichment(
                ranked_gene_df=diff_exp_results,
                gene_col='ID',         # The column with AGI codes in your data
                rank_col='t',      # The column to rank by (logFC is perfect for this)
                obodag=obodag,
                geneid2gos=geneid2gos,
                keys = list(ids),
                stress = stress,
                out_path = out_path
            )

            if not gsea_results_df.empty:
                print("\n--- GSEA Enrichment Results (FDR < 0.1) ---")
                print(gsea_results_df.head())
                
                # Save results to a file
                gsea_results_df.to_csv(f'{out_path}{stress}_gsea_go_enrichment_results.csv', index=False)
                print("\nResults saved to 'gsea_go_enrichment_results.csv'")
            else:
                print("\nNo significant GO terms found with FDR < 0.25.")

        except Exception as e:
            print(f"An error occurred during the analysis: {e}")