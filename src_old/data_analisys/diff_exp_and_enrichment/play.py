import gseapy
import pandas as pd
from goatools.obo_parser import GODag
from typing import Dict, Set, Union
from src.diff_exp_and_enrichment.pr_rank_gene_enrich import get_go_data
import os

def create_cls_file(labels: list, output_dir: str, pheno_a: str, pheno_b: str) -> str:
    """
    Generates a GSEA-compatible .cls file from a list of labels.
    """
    os.makedirs(output_dir, exist_ok=True)
    cls_path = os.path.join(output_dir, "temp_phenotypes.cls")
    
    # Line 1: sample_count class_count 1
    line1 = f"{len(labels)} 2 1"
    
    # Line 2: # phenotype_A phenotype_B
    # Clean phenotypes to be single words
    pheno_a_clean = pheno_a.replace(' ', '_')
    pheno_b_clean = pheno_b.replace(' ', '_')
    line2 = f"# {pheno_a_clean} {pheno_b_clean}"
    
    # Line 3: a a a b b b ...
    line3_labels = [label.replace(' ', '_') for label in labels]
    line3 = " ".join(line3_labels)
    
    with open(cls_path, 'w') as f:
        f.write(f"{line1}\n")
        f.write(f"{line2}\n")
        f.write(f"{line3}\n")
        
    return cls_path

def perform_gsea_analysis(
    expression_df: pd.DataFrame,
    class_vector: Union[str, list],
    obodag: GODag,
    geneid2gos: Dict[str, Set[str]],
    stress: str,
    out_path: str
) -> pd.DataFrame:
    """
    Performs Gene Set Enrichment Analysis (GSEA) using raw expression data.
    ... (docstring remains the same) ...
    """
    # 1. Prepare the Gene Sets (GO terms) for GSEApy
    print("Creating GO gene sets for GSEA...")
    go_gene_sets = {}
    for gene, go_ids in geneid2gos.items():
        for go_id in go_ids:
            if go_id not in go_gene_sets:
                go_gene_sets[go_id] = []
            go_gene_sets[go_id].append(gene.upper())

    go_gene_sets_named = {
        f"{go_id} ({obodag[go_id].name if go_id in obodag else 'Unknown'})": genes
        for go_id, genes in go_gene_sets.items()
    }
    print(f"Created {len(go_gene_sets_named)} GO term gene sets.")

    # 2. Run GSEA
    print("\nRunning GSEA analysis...")
    
    # Ensure gene names in expression data are uppercase to match gene sets
    expression_df.index = expression_df.index.str.upper()

    # --- FIX STARTS HERE ---
    # Prepare the dataframe in the format gseapy expects:
    # Reset index to move gene names from the index to a column.
    # The GSEA function internally handles the first column as gene names.
    data_for_gsea = expression_df.reset_index()
    # --- FIX ENDS HERE ---

    gs_res = gseapy.gsea(
        data=data_for_gsea,  # Pass the modified dataframe
        gene_sets=go_gene_sets_named,
        cls=class_vector,
        method='log2_ratio_of_classes',
        permutation_type='phenotype',
        min_size=15,
        max_size=500,
        outdir=f'{out_path}/{stress}_gsea_results',
        verbose=True
    )
    
    # ... (rest of the function remains the same) ...
    print("GSEA complete. Formatting results...")
    results_df = gs_res.res2d

    significant_results = results_df[results_df['fdr'] < 0.25].copy()

    if significant_results.empty:
        print("No significant GO terms found with FDR < 0.25.")
        return pd.DataFrame()

    significant_results.rename(columns={
        'es': 'enrichment_score',
        'nes': 'normalized_enrichment_score',
        'pval': 'p_value',
        'fdr': 'fdr_q_value',
        'ledge_genes': 'leading_edge_genes'
    }, inplace=True)
    
    significant_results['go_id'] = significant_results.index.str.split(' ').str[0]
    
    final_df = significant_results[[
        'go_id',
        'normalized_enrichment_score',
        'p_value',
        'fdr_q_value',
        'leading_edge_genes'
    ]]
    
    return final_df.sort_values(by='normalized_enrichment_score', ascending=False)



if __name__ == "__main__":
    # --- Setup ---
    constant_data_path = 'data/'
    GO_OBO_FILE = f'{constant_data_path}/go-basic.obo'
    ANNOTATION_FILE = f'{constant_data_path}/goa_arabidopsis.gaf.gz'
    save_dir:str = '/home/alex/Documents/GitHub/Data_collection/df_final'
    data_type = '2_way_norm'
    
    # --- Load and Pre-process Data (once) ---
    norm_exp_mat = pd.read_csv(f'{save_dir}/{data_type}.csv', index_col='ID_REF')
    norm_exp_mat.columns = [col.split('_')[0] for col in norm_exp_mat.columns]
    
    labels = pd.read_csv('outputs/data_outputs/1.1-complete-2_way_norm/labels.csv')
    labels['TREATMENT_CLEAN'] = labels['TREATMENT'].str.replace(r"[\[\]']", "", regex=True)

    # Load GO data
    obodag, geneid2gos = get_go_data(GO_OBO_FILE, ANNOTATION_FILE)

    # --- Loop Through Each Treatment and Run GSEA vs. Control ---
    control_group = 'No stress'
    all_treatments = labels['TREATMENT_CLEAN'].unique()
    stress_treatments = [t for t in all_treatments if t != control_group]

    for treatment in stress_treatments:
        print(f"\n{'='*20}\nRunning GSEA for: {treatment} vs. {control_group}\n{'='*20}")

        try:
            # 1. Filter labels for the current comparison
            current_labels = labels[labels['TREATMENT_CLEAN'].isin([treatment, control_group])].copy()

            # 2. Find and filter for common samples
            common_samples = list(set(norm_exp_mat.columns) & set(current_labels['sample_id']))
            
            if len(common_samples) < 6: # GSEA requires at least 3 per group
                print(f"Skipping {treatment}: Not enough common samples found (found {len(common_samples)}).")
                continue

            current_expression_data = norm_exp_mat[common_samples]
            current_labels = current_labels[current_labels['sample_id'].isin(common_samples)]
            
            # Data Validation
            if current_expression_data.index.duplicated().any():
                num_dups = current_expression_data.index.duplicated().sum()
                print(f"WARNING: Found {num_dups} duplicate gene IDs. Keeping first instance.")
                current_expression_data = current_expression_data[~current_expression_data.index.duplicated(keep='first')]
            if current_expression_data.isnull().values.any():
                num_nans = current_expression_data.isnull().values.sum()
                print(f"WARNING: Found {num_nans} NaN values. Dropping affected genes.")
                current_expression_data.dropna(inplace=True)
            if len(current_expression_data) < 15:
                print(f"Skipping {treatment}: Not enough valid genes remaining after cleaning.")
                continue
            
            # 3. Sort to ensure perfect alignment
            current_expression_data = current_expression_data.sort_index(axis=1)
            current_labels = current_labels.sort_values(by='sample_id').reset_index(drop=True)

            # 4. ⭐ ROBUST FIX: Create a .cls file for GSEA ⭐
            # The order of phenotypes matters for comparison, get it from the labels
            phenotypes = current_labels['TREATMENT_CLEAN'].unique().tolist()
            # Ensure the treatment is the first phenotype for consistent comparison
            if phenotypes[0] == control_group:
                phenotypes = phenotypes[::-1]
            
            output_directory = f'play_test/{treatment.replace(" ", "_")}_gsea_results'
            
            cls_file_path = create_cls_file(
                labels=current_labels['TREATMENT_CLEAN'].tolist(),
                output_dir=output_directory,
                pheno_a=phenotypes[0],
                pheno_b=phenotypes[1]
            )

            # 5. Run the GSEA function using the path to the .cls file
            perform_gsea_analysis(
                expression_df=current_expression_data,
                class_vector=cls_file_path,  # Pass the file path string
                obodag=obodag,
                geneid2gos=geneid2gos,
                stress=treatment.replace(" ", "_"),
                out_path='play_test/'
            )
            print(f"✅ Successfully completed GSEA for {treatment}.")

        except Exception as e:
            print(f"❌ ERROR processing {treatment}: {e}")