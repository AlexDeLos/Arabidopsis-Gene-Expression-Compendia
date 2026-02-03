from Bio import Entrez
import sys
import argparse
# Ensure we can import the local modules
module_dir = './'
sys.path.append(module_dir)

# Import the processing function from your other file
from src.data_importing.microarray_data_processing import Microarray_tracker,download_experiments_microarray
from src.data_importing.RNA_seq_processing import download_experiments_RNA_seq_nf_core
from src.data_importing.helpers.helpers import plot_tracker_results,plot_tracker_results_RNA, combine_files_microarray,plot_study_distributions_incremental,plot_study_distributions_seaborn
from src.data_importing.helpers.download_helper import search_geo_accessions
from src.constants import *

from src.data_importing.helpers.file_tracker import FileTracker
from src.data_importing.helpers.scan_tracker import RNASeq_tracker

# --- CONFIGURATION ---
Entrez.email = "your.email@example.com" 

MICROARRAY_QUERY = '"Arabidopsis thaliana"[Organism] AND "Expression profiling by array"[DataSet Type] AND "GSE"[Entry Type]'
RNASEQ_QUERY = '"Arabidopsis thaliana"[Organism] AND "Expression profiling by high throughput sequencing"[DataSet Type] AND "GSE"[Entry Type]'
import GEOparse
import os
import requests
import re

def download_processed_counts(gse_id, output_dir):
    print(f"Checking metadata for {gse_id}...")
    
    # 1. Parse the GEO metadata (downloads a small XML/Soft file)
    try:
        gse = GEOparse.get_GEO(geo=gse_id, destdir=output_dir, silent=True)
    except Exception as e:
        print(f"Error connecting to GEO: {e}")
        return False

    # 2. Look for Supplementary Files
    # The metadata contains links to files uploaded by authors
    if 'supplementary_files' not in gse.metadata:
        print(f"No supplementary files found for {gse_id}.")
        return False

    downloaded = False
    
    for url in gse.metadata['supplementary_files']:
        filename = url.split('/')[-1]
        
        # 3. Filter: We only want Count Matrices, not raw tars or READMEs
        # Common patterns for count matrices: txt, csv, tsv, xls, tab, count
        if re.search(r'(count|fpkm|tpm|expression|matrix|table)', filename, re.IGNORECASE):
            if re.search(r'(tar|xml|json)', filename, re.IGNORECASE):
                continue # Skip archives or metadata

            print(f"  > Found candidate: {filename}")
            
            # 4. Download
            save_path = os.path.join(output_dir, filename)
            try:
                response = requests.get(url, stream=True)
                if response.status_code == 200:
                    with open(save_path, 'wb') as f:
                        for chunk in response.iter_content(chunk_size=1024):
                            f.write(chunk)
                    print(f"    Downloaded: {filename}")
                    downloaded = True
                else:
                    print(f"    Failed to download link.")
            except Exception as e:
                print(f"    Download Error: {e}")

    return downloaded

# Usage
# ids = ['GSE77815', 'GSE44053'] # Your list
# for i in ids:
#     download_processed_counts(i, "./count_data_folder")
# --- MAIN EXECUTION ---
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--out_dir", help="output_dir", default='./new_storage/')
    parser.add_argument("--array_index", type=int, default=None, help="SLURM Array Task ID")
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
    # tracker_loc = f"{root_storage_dir}rnaseq_data/RNA_tracker_stats_temp.json"
    # RNA_tracker = RNASeq_tracker.load_from_json(tracker_loc)


    file_tracker_loc = f"{root_storage_dir}rnaseq_data/file_tracker/"
    rnaseq_ids = search_geo_accessions(RNASEQ_QUERY, max_results=50, filter_organism="Arabidopsis thaliana")#= ['GSE299572']# 
    RNA_tracker = FileTracker(file_tracker_loc)
    if args.array_index is not None:
        # --- PARALLEL MODE ---
        # Python checks if the index is valid
        if 0 <= args.array_index < len(rnaseq_ids):
            target_id = rnaseq_ids[args.array_index]
            print(f"--- ARRAY JOB: Processing ID #{args.array_index}: {target_id} ---")
            # Process ONLY this one ID
            download_experiments_RNA_seq_nf_core(target_id,root_storage_dir, f"{root_storage_dir}rnaseq_data",RNA_tracker, download_raw=True, scan=False,run_and_delete=False)
        else:
            print(f"Index {args.array_index} is out of bounds for {len(rnaseq_ids)} studies.")
    else:
        # --- SERIAL MODE (Original behavior) ---
        download_experiments_RNA_seq_nf_core(rnaseq_ids,root_storage_dir, f"{root_storage_dir}rnaseq_data",RNA_tracker, download_raw=True, scan=False,run_and_delete=False)
    
    print("\nDone!")