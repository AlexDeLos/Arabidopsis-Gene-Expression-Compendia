from Bio import Entrez
import sys

# Ensure we can import the local modules
module_dir = './'
sys.path.append(module_dir)

# Import the processing function from your other file
from src.data_importing.microarray_data_processing import Data_processing, Microarray_tracker,RNASeq_tracker
from src.data_importing.helpers import plot_tracker_results,plot_tracker_results_RNA, combine_files_microarray,plot_intensity_distributions
from src.data_importing.download_helper import search_geo_accessions,download_experiments_microarray,download_experiments_RNA_seq

# --- CONFIGURATION ---
Entrez.email = "your.email@example.com" 

MICROARRAY_QUERY = '"Arabidopsis thaliana"[Organism] AND "Expression profiling by array"[DataSet Type] AND "GSE"[Entry Type]'
RNASEQ_QUERY = '"Arabidopsis thaliana"[Organism] AND "Expression profiling by high throughput sequencing"[DataSet Type] AND "GSE"[Entry Type]'

# --- MAIN EXECUTION ---
if __name__ == "__main__":
    root_storage_dir = 'new_storage/'
    # 1. Microarray Data
    print("--- STARTING MICROARRAY SEARCH ---")
    
    # Initialize Tracker

    # ma_tracker = Microarray_tracker.load_from_json("tracker_stats.json")
    ma_tracker = Microarray_tracker()
    ma = 100
    scan_folder = f'.{root_storage_dir}microarray_scan_{ma}/'
    processed_folder = f'{root_storage_dir}processed_microarray_data_{ma}/'
    microarray_ids = search_geo_accessions(MICROARRAY_QUERY, max_results=ma)
    
    # # Pass tracker to the download function
    valid_microarray_ids = download_experiments_microarray(microarray_ids, f"{root_storage_dir}microarray_data/", ma_tracker, download_raw=True, scan=False,output_folder=processed_folder)
    ma_tracker.print_summary()
    ma_tracker.save_to_json(f"{root_storage_dir}{scan_folder}tracker_stats.json")
    plot_tracker_results(f"{root_storage_dir}{scan_folder}tracker_stats.json", output_dir= scan_folder)

    combined = combine_files_microarray(processed_folder, "RMA_Microarray_Combined.csv", f"{root_storage_dir}final_data")
    # combined = pd.read_csv(f'{root_storage}final_data/RMA_Microarray_Combined.csv')

    plot_intensity_distributions(combined, "new_storage/intensity_plots")



    # print("\n--- STARTING RNA-SEQ SEARCH ---")
    # RNA_tracker = RNASeq_tracker()
    # rnaseq_ids = search_geo_accessions(RNASEQ_QUERY, max_results=10000)
    # download_experiments_RNA_seq(rnaseq_ids, f"{root_storage}/rnaseq_data",RNA_tracker, download_raw=True, scan=True)
    # RNA_tracker.print_summary()
    # RNA_tracker = RNASeq_tracker.load_from_json(f"{root_storage}/RNA_seq_scan/RNA_tracker_stats.json")
    # RNA_tracker.save_to_json(f"{root_storage}/RNA_seq_scan/RNA_tracker_stats.json")
    # plot_tracker_results_RNA(f"{root_storage}/RNA_seq_scan/RNA_tracker_stats.json",output_dir=f'{root_storage}/RNA_seq_scan')
    print("\nDone!")