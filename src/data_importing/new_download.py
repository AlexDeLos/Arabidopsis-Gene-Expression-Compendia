from Bio import Entrez
import sys
import argparse
import pandas as pd
# Ensure we can import the local modules
module_dir = './'
sys.path.append(module_dir)

# Import the processing function from your other file
from src.data_importing.microarray_data_processing import Microarray_tracker,download_experiments_microarray
from src.data_importing.RNA_seq_processing import RNASeq_tracker, download_experiments_RNA_seq
from src.data_importing.helpers import plot_tracker_results,plot_tracker_results_RNA, combine_files_microarray,plot_study_distributions_incremental,plot_study_distributions_seaborn
from src.data_importing.download_helper import search_geo_accessions
from src.constants import *

# --- CONFIGURATION ---
Entrez.email = "your.email@example.com" 

MICROARRAY_QUERY = '"Arabidopsis thaliana"[Organism] AND "Expression profiling by array"[DataSet Type] AND "GSE"[Entry Type]'
RNASEQ_QUERY = '"Arabidopsis thaliana"[Organism] AND "Expression profiling by high throughput sequencing"[DataSet Type] AND "GSE"[Entry Type]'

# --- MAIN EXECUTION ---
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--out_dir", help="output_dir", default='new_storage/')
    args = parser.parse_args()

    root_storage_dir = args.out_dir
    # 1. Microarray Data
    print("--- STARTING MICROARRAY SEARCH ---")
    
    # Initialize Tracker

    # ma_tracker = Microarray_tracker.load_from_json("tracker_stats.json")
    ma_tracker = Microarray_tracker()
    ma = 20000
    scan_folder = f'.{root_storage_dir}microarray_scan_{ma}/'
    processed_folder = f'{root_storage_dir}processed_microarray_data/'
    downloads_folder = f"{root_storage_dir}microarray_data/"
    # microarray_ids = search_geo_accessions(MICROARRAY_QUERY, max_results=ma)
    
    # # Pass tracker to the download function
    # ma_tracker.sync_with_filesystem(downloads_folder,processed_folder)
    # saved_tracker = Microarray_tracker.load_from_json('new_storage/microarray_data/tracker_stats.json')
    # test = ma_tracker.compare_states(saved_tracker)
    # valid_microarray_ids = download_experiments_microarray(microarray_ids, downloads_folder, saved_tracker, download_raw=True, scan=False,output_folder=processed_folder)
    # ma_tracker.print_summary()
    # ma_tracker.save_to_json(f"{root_storage_dir}{scan_folder}tracker_stats.json")
    # plot_tracker_results(f"{root_storage_dir}{scan_folder}tracker_stats.json", output_dir= scan_folder)

    # combined,map = combine_files_microarray(processed_folder, "RMA_Microarray_Combined.csv", f"{root_storage_dir}final_data",combination_method='max',combine_genes=True)
    # combined = pd.read_csv(f'{root_storage_dir}final_data/RMA_Microarray_Combined.csv')

    # plot_study_distributions_seaborn(processed_folder, "new_storage/new_plots/intensity/intensity_plot_matplot_test.svg")
    # plot_study_distributions_incremental(processed_folder, "new_storage/new_plots/intensity/intensity_plot_incremental")
    # raise ValueError("DONE")



    # print("\n--- STARTING RNA-SEQ SEARCH ---")
    RNA_tracker = RNASeq_tracker()
    rnaseq_ids = search_geo_accessions(RNASEQ_QUERY, max_results=10, filter_organism="Arabidopsis thaliana")#= ['GSE299572']# 
    RNA_tracker = RNASeq_tracker.load_from_json('new_storage/rnaseq_data/rnaseq_tracker_stats.json')
    download_experiments_RNA_seq(rnaseq_ids,root_storage_dir, f"{root_storage_dir}/rnaseq_data",RNA_tracker, download_raw=True, scan=False,run_and_delete=False)
    # RNA_tracker.print_summary()
    #
    RNA_tracker.save_to_json(f"{root_storage_dir}/RNA_seq_scan/RNA_tracker_stats.json")
    plot_tracker_results_RNA(f"{root_storage_dir}/RNA_seq_scan/RNA_tracker_stats.json",output_dir=f'{root_storage_dir}/RNA_seq_scan')
    print("\nDone!")