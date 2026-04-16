import pandas as pd
from goatools.obo_parser import GODag
from goatools.go_enrichment import GOEnrichmentStudy
from typing import List, Dict, Set
import gzip
import Bio.UniProt.GOA as GOA
import os
from ftplib import FTP
import re
import matplotlib.pyplot as plt
import numpy as np

def parse_gaf_for_agi_codes(annotation_file: str) -> Dict[str, Set[str]]:
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

def perform_arabidopsis_enrichment(gene_list: List[str], go_obo_file: str, annotation_file: str) -> pd.DataFrame:
    """
    Performs gene enrichment analysis for a list of Arabidopsis thaliana genes.
    """
    print("Parsing GO OBO file...")
    try:
        obodag = GODag(go_obo_file, optional_attrs={'relationship'})
    except Exception as e:
        print(f"Error parsing OBO file: {e}")
        return pd.DataFrame()

    print("Parsing Arabidopsis annotation file...")
    try:
        if not os.path.isfile(annotation_file):
            print(f"Annotation file '{annotation_file}' not found. Downloading...")
            with FTP('ftp.ebi.ac.uk') as ebi_ftp:
                ebi_ftp.login()
                arab_uri = '/pub/databases/GO/goa/ARABIDOIS/goa_arabidopsis.gaf.gz'
                with open(annotation_file, 'wb') as arab_fp:
                    ebi_ftp.retrbinary(f'RETR {arab_uri}', arab_fp.write)
            print("Download complete.")

        geneid2gos = parse_gaf_for_agi_codes(annotation_file)
        if not geneid2gos:
            print("Error: Could not parse any AGI-to-GO mappings from the annotation file. Check the file format or parsing function.")
            return pd.DataFrame()

    except Exception as e:
        print(f"An error occurred during file parsing: {e}")
        return pd.DataFrame()

    population = set(geneid2gos.keys())
    study_genes = {gene.upper() for gene in gene_list}

    genes_in_pop = study_genes.intersection(population)
    print(f"\nOriginal gene list size: {len(study_genes)}")
    if study_genes:
        print(f"Genes found in annotation population: {len(genes_in_pop)} ({len(genes_in_pop)/len(study_genes):.1%})")

    if not genes_in_pop:
        print("Warning: None of the genes in your list were found in the annotation database.")
        return pd.DataFrame()

    print("\nPerforming GO enrichment study...")
    goeaobj = GOEnrichmentStudy(
        population,
        geneid2gos,
        obodag,
        propagate_counts=True,
        alpha=0.05,
        methods=['fdr_bh']
    )
    # use Gene Set Enerichment Analysis gene enrichment 
    #! order maters DO IT WITH THIS

    results = goeaobj.run_study(genes_in_pop)
    
    if not results:
        return pd.DataFrame()

    results_df = pd.DataFrame([
        {
            "GO ID": r.GO, "GO Name": r.name, "Namespace": r.NS,
            "p-value": r.p_uncorrected, "Adjusted p-value": r.p_fdr_bh,
            "Ratio in Study": f"{r.study_count}/{len(genes_in_pop)}",
            "Ratio in Population": f"{r.pop_count}/{len(population)}",
            "Study Genes": ", ".join(sorted(list(r.study_items)))
        } for r in results
    ])

    return results_df[results_df['Adjusted p-value'] < 0.5].sort_values(by='Adjusted p-value')
    # return results_df.sort_values(by='Adjusted p-value')

def run_test():
    """Runs a test with a known list of cold-related genes."""
    print("--- Running Built-in Test ---")
    test_genes = ['AT4G25470', 'AT4G25480', 'AT4G25490', 'AT1G20440', 'AT3G50970']
    print(f"Test gene list: {test_genes}")
    
    enrichment_results = perform_arabidopsis_enrichment(test_genes, 'data/go-basic.obo', 'data/goa_arabidopsis.gaf.gz')

    if not enrichment_results.empty:
        print("\n✅ Test Passed! Found significant enrichment results.")
        cold_related = enrichment_results[enrichment_results['GO Name'].str.contains("cold|stress|temperature", case=False)]
        print(cold_related[['GO ID', 'GO Name', 'Adjusted p-value']].head())
    else:
        print("\n❌ Test Failed. No significant enrichment results were found for the test list.")
    print("--- Test Complete ---\n")

def plot_gene_enrichment(final_results,path,treatment,tissue):

    x = list(final_results['GO Name'])
    y = -np.log(np.array(list(final_results['Adjusted p-value'])))

    plt.bar(x, y)
    
    plt.xticks(range(len(x)), x, rotation='vertical')

    plt.xlabel('Enriched terms')
    plt.ylabel('-log(adjusted P-value)')
    plt.title(f'{tissue} exposef to {treatment}.svg')
    plt.tight_layout()
    plt.savefig(path+f'/{treatment}_{tissue}')
    # x = 0
    plt.close()

def plot_gene_enrichment_with_ratios(final_results, path, treatment, tissue):
    """
    Generates and saves a bar plot of gene enrichment results, with ratios
    annotated on top of each bar.

    Args:
        final_results (pd.DataFrame): DataFrame with enrichment results.
                                      Must contain 'GO Name', 'Adjusted p-value',
                                      'Ratio in Study', and 'Ratio in Population'.
        path (str): The directory path to save the plot.
        treatment (str): The treatment name for the plot title and filename.
        tissue (str): The tissue name for the plot title and filename.
    """
    # Ensure the output directory exists
    if not os.path.exists(path):
        os.makedirs(path)

    # Prepare data for plotting
    # Sorting by p-value makes the plot more intuitive
    df_sorted = final_results.sort_values('Adjusted p-value').head(15) # Plot top 15 for readability
    
    go_terms = df_sorted['GO Name']
    y_values = -np.log10(df_sorted['Adjusted p-value']) # Use log10 for standard visualization
    ratio_study = df_sorted['Ratio in Study']
    ratio_pop = df_sorted['Ratio in Population']

    # Create figure and axes for more control
    fig, ax = plt.subplots(figsize=(12, 8))
    bars = ax.bar(go_terms, y_values, color='skyblue')

    # --- Code to add labels on top of each bar ---
    for i, bar in enumerate(bars):
        height = bar.get_height()
        # Format the text with the two ratios
        label_text = f"{ratio_study.iloc[i]}\n{ratio_pop.iloc[i]}"
        
        ax.text(
            bar.get_x() + bar.get_width() / 2.0,  # X position (center of the bar)
            height,                               # Y position (top of the bar)
            label_text,                           # The text to display
            ha='center',                          # Horizontal alignment
            va='bottom',                          # Vertical alignment
            fontsize=8,                           # Adjust font size for readability
            rotation=90,                          # Rotate text to fit better
            color='black'
        )

    # --- Formatting the plot for readability ---
    ax.set_xticklabels(go_terms, rotation=45, ha='right') # Rotate labels for better fit
    ax.set_xlabel('Enriched GO Terms', fontsize=12)
    ax.set_ylabel('-log10(Adjusted P-value)', fontsize=12)
    ax.set_title(f'GO Enrichment for {tissue} exposed to {treatment}', fontsize=14)
    
    # Adjust y-axis limit to make space for the text
    ax.set_ylim(top=ax.get_ylim()[1] * 1.2)
    
    fig.tight_layout() # Adjust layout to prevent labels from being cut off

    # Save the figure
    file_path = os.path.join(path, f'{treatment}_{tissue}_improved.svg')
    plt.savefig(file_path, format='svg')
    plt.close(fig) # Close the figure to free up memory

def evaluate_diff_exp(constant_data_path,diff_exp_path,output_path,plot_path,keys:list=None, ids:list=None):
    for file in os.listdir(diff_exp_path):
        if file.endswith('.csv'):
            gene_list_path = os.path.join(diff_exp_path, file)

            treatment_name:str = gene_list_path.split('/')[-1].split('_')[1]
            tissue_name:str = gene_list_path.split('/')[-1].split('_')[0]
            if tissue_name !='':
                continue
            try:
                print("--- Running Analysis on Your Gene List ---")
                user_genes_df = pd.read_csv(gene_list_path)
                user_genes = list(user_genes_df['ID'])
                
                final_results = perform_arabidopsis_enrichment(user_genes, f'{constant_data_path}/go-basic.obo', f'{constant_data_path}/goa_arabidopsis.gaf.gz')
                # TODO: Keep only the enrichment results that are relevant to the treatment set we are using
                keywords:list =keys# stress_keyword_map.get(treatment_name)
                # for n in stress_keyword_map.keys():
                #     keywords.extend(stress_keyword_map.get(n))
                
                if keywords:
                    # print(f"Filtering for keywords: {keywords}")
                    pattern = '|'.join(keywords)
                    final_results = final_results[final_results['GO Name'].str.contains(pattern, case=False, na=False)]
                elif ids:
                    pattern = '|'.join(ids)
                    final_results = final_results[final_results['GO ID'].str.contains(pattern, case=False, na=False)]
                else:
                    print(f"Warning: No keywords defined for treatment '{treatment_name}'. Saving all significant results.")
                    final_results = final_results

                if not final_results.empty:
                    print("\nSignificant Enrichment Results from Your List:")
                    print(final_results)
                    final_results.to_csv(f"{output_path}/{tissue_name}_{treatment_name}_enrichment_results.csv", index=False)
                    # plot_gene_enrichment(final_results,plot_path,treatment_name,tissue_name)
                    plot_gene_enrichment_with_ratios(final_results,plot_path,treatment_name,tissue_name)
                    print("\nResults saved to enrichment_results.csv")
                else:
                    print("\nNo significant enrichment results were found for your gene list.")
                    
            except FileNotFoundError:
                print(f"\nError: {gene_list_path} not found. Please place it in the same directory.")

if __name__ == '__main__':
    # run_test()

    constant_data_path:str = 'data'
    data_type = '2_way_norm'
    diff_exp_path = f'diff_exp/1.0_0.05_{data_type}'
    output_path = f'enrichment_results/new_{data_type}'
    # d = pd.read_csv(f'enrichment_results/corrected/_Chemical Stress_enrichment_results.csv')
    plot_path = f'{output_path}/plots'
    os.makedirs(plot_path,exist_ok=True)
    go_terms = pd.read_excel(f'{constant_data_path}/GO_terms_to_enrich_2.xlsx')
    keys = None#go_terms['GO Term Name']
    ids =  go_terms['GO Term ID']
    evaluate_diff_exp(constant_data_path,diff_exp_path,output_path,plot_path,keys,list(ids))
    # pd.read_csv('enrichment_results_full/Drought Stress_enrichment_results.csv')