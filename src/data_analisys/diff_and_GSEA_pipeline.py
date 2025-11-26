import pandas as pd
import numpy as np
import sys
import os # <-- Added this import
module_dir = './'
sys.path.append(module_dir)
from src.constants import *
from src.data_analisys.diff_exp_and_enrichment.pr_rank_gene_enrich import get_go_data,perform_gsea_enrichment
from src.data_analisys.diff_exp_and_enrichment.plot_enrichment_new import plot_enrichment_scatter_interactive, create_gsea_spider_plot
from src.data_analisys.diff_exp_and_enrichment.diff_expr import diff_exp_combine_tissues


iterations = 1000

treatments = [
    "Drought Stress",
    "Salinity Stress",
    "Heat Stress",
    "Cold Stress",
    "High Light Stress",
    ]
STRESS_GO_ROOTS = {
    "GO:0009414", # response to drought
    "GO:0009651", # response to salt stress
    "GO:0009408", # response to heat
    "GO:0009409", # response to cold
    "GO:0009644", # response to high light intensity
    "GO:0009611", # response to wounding
    "GO:0001666", # response to hypoxia (flooding)
    "GO:0009411", # response to UV
    "GO:0009636", # response to toxic substance
}

STRESS_GO_ROOTS_small = {
    "GO:0009414": "Drought Stress", # response to drought
    "GO:0009651": "Salinity Stress", # response to salt stress
    "GO:0009408": "Heat Stress", # response to heat
    "GO:0009409": "Cold Stress", # response to cold
    "GO:0009644": "High Light Stress", # response to high light intensity
}
def get_spider_plots(path, results_path, data_types, Fulls, tissues, pures,filter,full=False):
    """
    Generates spider plots comparing a single GO term across all experiment combinations.
    """
    if not os.path.exists(path):
        os.makedirs(path)
        print(f"Created directory: {path}")
    print("\n--- Generating Spider Plots ---")
    for term in STRESS_GO_ROOTS_small.keys():
        print(f"Gathering data for GO term: {term}")
        all_term_rows = [] # List to store rows for the current term from all experiments

        # Iterate through all experiment combinations to find the GSEA results
        # for fil in filter_low_combination:
        for data_type in data_types:
                # for pure in pures:
            for Full in Fulls:
                for tissue in tissues:
                    # if (not Full) and (tissue is None):
                    #     break
                    # Reconstruct the experiment name and path from the main loop
                    exp_name = f'{EXPERIMENT_NAME}_{data_type}_{tissue if tissue else 'All_tissues'}_{'full' if Full else 'sanity'}_{'pure' if pures else 'mixed'}_min_group_{filter}'
                    out_path = f'{results_path}GSEA_enrichment_{data_type}_{exp_name}/'

                    stress = STRESS_GO_ROOTS_small[term]
                    csv_file = f'{out_path}{stress}_gsea_go_enrichment_results_{iterations}.csv'
                    
                    if os.path.isfile(csv_file):
                        try:
                            gsea_results_df = pd.read_csv(csv_file)
                            # Find the row for the specific GO term we're plotting
                            term_row = gsea_results_df[gsea_results_df['go_id'] == term].copy()
                            
                            if not term_row.empty:
                                # This is the key step: create a unique, descriptive name 
                                # for this specific experiment to use in the plot legend.
                                legend_name = f'{data_type} {'Single Tissue' if tissue else "Combined Tissues"} {"Combined Studies" if Full else "Single study"}'
                                
                                # Set the 'Name' column, which create_gsea_spider_plot uses for labels
                                term_row.loc[:, 'Name'] = legend_name
                                all_term_rows.append(term_row)
                        except Exception as e:
                            print(f"Error reading or processing file {csv_file}: {e}")
                    else:
                        # This is normal if a combination was skipped or failed
                        pass 

        if not all_term_rows:
            print(f"No data found for GO term {term}. Skipping plot.")
            continue # Move to the next GO term

        # Combine all the found rows (one per experiment) into a single DataFrame
        df = pd.concat(all_term_rows, ignore_index=True)
        
        if df.empty:
            print(f"DataFrame is empty for term {term} after concatenation.")
            continue
        
        print(f"Generating spider plot for {term} with {len(df)} entries.")
        
        # Create the plot. We assume create_gsea_spider_plot saves the plot
        # as an HTML file (like plot_enrichment_scatter_interactive).
        try:
            clean_term = STRESS_GO_ROOTS_small[term].replace(':', '_')
            plot_filename = f"{path}{clean_term}_spider_plot{'_full' if full else ''}.svg"
            #TODO: check why this is needed.
            if full:
                plot = df[df['Name'].isin(['study_corrected Combined Tissues Combined Studies',
                                        'study_corrected Single Tissue Combined Studies',
                                        'study_corrected Single Tissue Single study',

                                        'imputed Combined Tissues Combined Studies',
                                        'imputed Single Tissue Combined Studies',
                                        'imputed Single Tissue Single study',

                                        '2_way_norm Combined Tissues Combined Studies',
                                        '2_way_norm Single Tissue Combined Studies',
                                        '2_way_norm Single Tissue Single study',

                                        'standardized Combined Tissues Combined Studies',
                                        'standardized Single Tissue Combined Studies',
                                        'standardized Single Tissue Single study',

                                        'robust Combined Tissues Combined Studies',
                                        'robust Single Tissue Combined Studies',
                                        'robust Single Tissue Single study'

                                        'tissue_normalized Combined Tissues Combined Studies',
                                        'tissue_normalized Single Tissue Combined Studies',
                                        'tissue_normalized Single Tissue Single study'])]

            else:
                plot = df[df['Name'].isin([
                                        'study_corrected Combined Tissues Combined Studies',
                                        'study_corrected Single Tissue Combined Studies',
                                        'study_corrected Single Tissue Single study',

                                        'tissue_normalized Combined Tissues Combined Studies',
                                        'tissue_normalized Single Tissue Combined Studies',
                                        'tissue_normalized Single Tissue Single study'

                                        'tissue_normalized_2 Combined Tissues Combined Studies',
                                        'tissue_normalized_2 Single Tissue Combined Studies',
                                        'tissue_normalized_2 Single Tissue Single study'

                                        'imputed Combined Tissues Combined Studies',
                                        'imputed Single Tissue Combined Studies',
                                        'imputed Single Tissue Single study'])]
            create_gsea_spider_plot(plot, save_path=plot_filename, term=STRESS_GO_ROOTS_small[term])
            print(f"Saved plot to {plot_filename}")
        except Exception as e:
            print(f"Error generating spider plot for {term}: {e}")
            print("DataFrame columns:", df.columns)
            print("DataFrame head:")
            print(df.head())


def run_diff_exp_and_enrichment(save_dir:str=PROCESSED_DATA_FOLDER,
                                data_types = ['tissue_normalized','study_corrected','imputed','2_way_norm','tissue_normalized_2'],#,],#,'standardized+'],#,'standardized','standardized', 'robust'],#'robust+','2_way_norm_og',,
                                pures = [True,False],
                                Fulls = [True,False],
                                filter_low_combination = [0,15],
                                tissues = [None,'leaf'],
                                just_plot=False):    
    for fil in filter_low_combination:
        for pure in pures:
            for data_type in data_types:
                for Full in Fulls:
                    for tissue in tissues:
                        exp_name = f'{EXPERIMENT_NAME}_{data_type}_{tissue if tissue else 'All_tissues'}_{'full' if Full else 'sanity'}_{'pure' if pure else 'mixed'}_min_group_{fil}'
                        
                        diff_exp_out_dir = f'{FIGURES_DIR}dif_expression_results/{exp_name}/'
                        if not os.path.isfile(f'{diff_exp_out_dir}/done.txt'):
                            if not just_plot:
                                if not Full:
                                    samples_GSE41935 = ['GSM1027688', 'GSM1027720', 'GSM1027701', 'GSM1027875', 'GSM1027745', 'GSM1027710', 'GSM1027860', 'GSM1027768', 'GSM1027781', 'GSM1027855', 'GSM1027682', 'GSM1027801', 'GSM1027843', 'GSM1027746', 'GSM1027854', 'GSM1027695', 'GSM1027690', 'GSM1027804', 'GSM1027730', 'GSM1027719', 'GSM1027773', 'GSM1027871', 'GSM1027691', 'GSM1027856', 'GSM1027738', 'GSM1027742', 'GSM1027800', 'GSM1027806', 'GSM1027815', 'GSM1027863', 'GSM1027776', 'GSM1027760', 'GSM1027882', 'GSM1027817', 'GSM1027876', 'GSM1027739', 'GSM1027721', 'GSM1027861', 'GSM1027704', 'GSM1027826', 'GSM1027799', 'GSM1027808', 'GSM1027759', 'GSM1027818', 'GSM1027830', 'GSM1027733', 'GSM1027751', 'GSM1027715', 'GSM1027697', 'GSM1027728', 'GSM1027785', 'GSM1027824', 'GSM1027835', 'GSM1027740', 'GSM1027767', 'GSM1027884', 'GSM1027716', 'GSM1027868', 'GSM1027880', 'GSM1027849', 'GSM1027698', 'GSM1027807', 'GSM1027689', 'GSM1027812', 'GSM1027680', 'GSM1027709', 'GSM1027862', 'GSM1027877', 'GSM1027832', 'GSM1027726', 'GSM1027839', 'GSM1027753', 'GSM1027718', 'GSM1027810', 'GSM1027797', 'GSM1027834', 'GSM1027761', 'GSM1027802', 'GSM1027732', 'GSM1027821', 'GSM1027793', 'GSM1027684', 'GSM1027853', 'GSM1027712', 'GSM1027703', 'GSM1027789', 'GSM1027692', 'GSM1027722', 'GSM1027833', 'GSM1027827', 'GSM1027699', 'GSM1027873', 'GSM1027735', 'GSM1027774', 'GSM1027841', 'GSM1027820', 'GSM1027743', 'GSM1027696', 'GSM1027837', 'GSM1027846', 'GSM1027744', 'GSM1027694', 'GSM1027850', 'GSM1027828', 'GSM1027805', 'GSM1027847', 'GSM1027794', 'GSM1027829', 'GSM1027811', 'GSM1027845', 'GSM1027842', 'GSM1027754', 'GSM1027678', 'GSM1027859', 'GSM1027787', 'GSM1027708', 'GSM1027852', 'GSM1027792', 'GSM1027874', 'GSM1027878', 'GSM1027883', 'GSM1027700', 'GSM1027848', 'GSM1027679', 'GSM1027879', 'GSM1027851', 'GSM1027756', 'GSM1027724', 'GSM1027764', 'GSM1027723', 'GSM1027750', 'GSM1027798', 'GSM1027775', 'GSM1027791', 'GSM1027795', 'GSM1027831', 'GSM1027731', 'GSM1027685', 'GSM1027870', 'GSM1027736', 'GSM1027836', 'GSM1027872', 'GSM1027734', 'GSM1027866', 'GSM1027725', 'GSM1027749', 'GSM1027757', 'GSM1027687', 'GSM1027788', 'GSM1027681', 'GSM1027858', 'GSM1027881', 'GSM1027840', 'GSM1027786', 'GSM1027763', 'GSM1027705', 'GSM1027778', 'GSM1027765', 'GSM1027867', 'GSM1027825', 'GSM1027752', 'GSM1027803', 'GSM1027790', 'GSM1027823', 'GSM1027809', 'GSM1027822', 'GSM1027783', 'GSM1027771', 'GSM1027755', 'GSM1027747', 'GSM1027814', 'GSM1027686', 'GSM1027844', 'GSM1027741', 'GSM1027706', 'GSM1027762', 'GSM1027780', 'GSM1027702', 'GSM1027777', 'GSM1027714', 'GSM1027782', 'GSM1027717', 'GSM1027865', 'GSM1027713', 'GSM1027758', 'GSM1027838', 'GSM1027796', 'GSM1027769', 'GSM1027857', 'GSM1027770', 'GSM1027707', 'GSM1027748', 'GSM1027813', 'GSM1027784', 'GSM1027766', 'GSM1027869', 'GSM1027772', 'GSM1027779', 'GSM1027816', 'GSM1027727', 'GSM1027693', 'GSM1027737', 'GSM1027864', 'GSM1027683', 'GSM1027819', 'GSM1027729', 'GSM1027711']
                                else:
                                    samples_GSE41935 = None                            
                                diff_exp_combine_tissues(treatments,save_dir,data_type,out_dir=diff_exp_out_dir,samples=samples_GSE41935,pure=pure,tissue=tissue,filter_low_combination=fil)
                            f = open(os.path.join(diff_exp_out_dir, 'done.txt'), 'w')
                            f.write('this diff exp has been completed.')
                            f.close()
                            print(f"DONE WITH {diff_exp_out_dir}")
                        else:
                            print(f"ALREADY DID {diff_exp_out_dir}")
                        # GSEA
                            # 1. Define file paths
                        constant_data_path = CORE_DATA_DIR
                        diff_exp_out_path = f'{FIGURES_DIR}GSEA_enrichment_results/GSEA_enrichment_{data_type}_{exp_name}/'
                        GO_OBO_FILE = f'{constant_data_path}/go-basic.obo'
                        ANNOTATION_FILE = f'{constant_data_path}/tair.gaf.gz'
                        
                        # 2. Load the necessary GO data first
                        # This only needs to be done once

                        # go_terms = pd.read_excel(f'{constant_data_path}GO_terms_to_enrich_2.xlsx')
                        for stress in treatments:
                            try:
                                if (not just_plot) and (not os.path.isfile(f'{diff_exp_out_path}{stress}_gsea_go_enrichment_results_{iterations}.csv')):
                                    STRESS_IDS = {key for key, val in STRESS_GO_ROOTS_small.items() if True}
                                    # STRESS_IDS =None
                                    obodag, geneid2gos = get_go_data(GO_OBO_FILE, ANNOTATION_FILE,stress_root_go_ids=STRESS_IDS)
                                    # filter
                                    obodag = {key: value for key, value in obodag.items() if key in STRESS_IDS}  
                                    diff_exp_results = pd.read_csv(f'{diff_exp_out_dir}{tissue if tissue is not None else 'All-Tissues'}_{stress}_genes.csv') #TODO: this is hardcoded from {output_filename} in {diff_exp_combine_tissues and diff_exp}

                                    ids = None
                                    # diff_exp_results['rank'] = abs(diff_exp_results['logFC']) *(-np.log(diff_exp_results['adj.P.Val']))
                                    diff_exp_results['rank'] = diff_exp_results['logFC'] *(-np.log10(diff_exp_results['adj.P.Val']))
                                    
                                    gsea_results_df = perform_gsea_enrichment(
                                        ranked_gene_df = diff_exp_results,
                                        gene_col = 'ID',         # The column with AGI codes in your data
                                        rank_col = 'rank',       # The column to rank by (logFC is perfect for this)
                                        obodag = obodag,
                                        geneid2gos = geneid2gos,
                                        keys = ids,
                                        stress = stress,
                                        out_path = diff_exp_out_path,
                                        permutations = iterations
                                    )
                                    gsea_results_df.sort_values(by=['FDR q-val'])
                                    print("\n--- GSEA Enrichment Results ---")
                                    print(gsea_results_df.head())
                                    
                                    # Save results to a file
                                    gsea_results_df.to_csv(f'{diff_exp_out_path}{stress}_gsea_go_enrichment_results_{iterations}.csv', index=False)
                                    # del gsea_results_df
                                    print(f"\nResults saved to {f'{diff_exp_out_path}{stress}_gsea_go_enrichment_results_{iterations}.csv'}")
                                else:
                                    print(f"\nLoading pre existing results from: {f'{diff_exp_out_path}{stress}_gsea_go_enrichment_results_{iterations}.csv'}")
                                    gsea_results_df = pd.read_csv(f'{diff_exp_out_path}{stress}_gsea_go_enrichment_results_{iterations}.csv')
                                plot_enrichment_scatter_interactive(gsea_results_df,save_path=f'{FIGURES_DIR}plots_enrichment/{exp_name.split('_')[0]}/{'full' if Full else 'sanity'}/{tissue if tissue is not None else 'All-Tissues'}/{data_type}/{fil}/{'pure' if pure else 'mixed'}/{stress}.html',
                                                                    title=f'GSEA for {stress} on {tissue if tissue is not None else 'All-Tissues'} on {'full' if Full else 'sanity'} with {'pure' if pure else 'mixed'} treatments with a filter of >{fil}',treatments = treatments, normalizations=data_types)


                            except Exception as e:
                                print(f"An error occurred during the analysis: {e}")
    
            # Updated call to get_spider_plots, passing all the loop parameters
            get_spider_plots(
                path=f'{FIGURES_DIR}GSEA_radar_comparison/{EXPERIMENT_NAME}_{fil}_{pure}/',
                results_path=f'{FIGURES_DIR}GSEA_enrichment_results/',
                data_types=data_types,
                Fulls=Fulls,
                tissues=tissues,
                pures=pure,
                filter=fil,
                full=False
            )

    print("DONE with enrichment pipeline")



if __name__ == "__main__":
    run_diff_exp_and_enrichment(just_plot=False)