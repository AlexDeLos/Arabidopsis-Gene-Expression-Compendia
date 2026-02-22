from Bio import Entrez
import sys
import argparse
# Ensure we can import the local modules
module_dir = './'
sys.path.append(module_dir)

from src.constants import *

# Import the processing function from your other file
from src.data_importing.microarray_data_processing import Microarray_data_processing, Microarray_tracker,download_experiments_microarray
from src.data_importing.RNA_seq_processing_batch import download_experiments_RNA_seq_nf_core
from src.data_importing.helpers.helpers import plot_tracker_results,plot_tracker_results_RNA, combine_files_microarray,plot_study_distributions_incremental,plot_study_distributions_seaborn
from src.data_importing.helpers.download_helper import search_geo_accessions
from src.data_importing.helpers.file_tracker import FileTracker
# from src.data_importing.helpers.scan_tracker import RNASeq_tracker

# --- CONFIGURATION ---
Entrez.email = "your.email@example.com" 

def save_list_to_txt(data_list, filename):
    """Saves a list to a text file, one item per line."""
    with open(filename, 'w') as f:
        for item in data_list:
            f.write(f"{item}\n")
    print(f"Saved {len(data_list)} IDs to {filename}")

def load_list_from_txt(filename):
    """Reads a text file into a list, stripping newlines."""
    if not os.path.exists(filename):
        return []
    
    with open(filename, 'r') as f:
        # .strip() removes the \n character from the end
        return [line.strip() for line in f]

MICROARRAY_QUERY = '"Arabidopsis thaliana"[Organism] AND "Expression profiling by array"[DataSet Type] AND "GSE"[Entry Type]'#cel"[Supplementary Files]'
RNASEQ_QUERY = '"Arabidopsis thaliana"[Organism] AND "Expression profiling by high throughput sequencing"[DataSet Type] AND "GSE"[Entry Type]'
STRESS_QUERY = ' AND ("stress"[Title] OR "response"[Title] OR "abiotic"[Title] OR "biotic"[Title])'
FULL_QUERY_RNA = RNASEQ_QUERY + STRESS_QUERY
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
    parser.add_argument("-b", "--batch_size", help="output_dir", default=2,type=int)
    parser.add_argument("--array_index", type=int, default=0, help="SLURM Array Task ID")
    parser.add_argument("--ma", action="store_true", default=False)
    parser.add_argument("--rna", action="store_true", default=False)
    args = parser.parse_args()

    root_storage_dir = args.out_dir
    # 1. Microarray Data
    run_microarray:bool = args.ma
    run_rna_seq:bool = args.rna
    
    # Initialize Tracker
    if run_microarray:
        print("--- STARTING MICROARRAY SEARCH ---")
        ma = 20000
        scan_folder = f'.{root_storage_dir}microarray_scan/'
        processed_folder = f'{root_storage_dir}processed_microarray_data/'
        downloads_folder = f"{root_storage_dir}microarray_data/"
        # microarray_ids = search_geo_accessions(MICROARRAY_QUERY, max_results=ma)
        filename = './microarray_ids.txt'
        # save_list_to_txt(microarray_ids,filename)
        microarray_ids = load_list_from_txt(filename)
        
        # Pass tracker to the download function
        # ma_tracker.sync_with_filesystem(downloads_folder,processed_folder)
        saved_tracker = Microarray_tracker.load_from_json(root_storage_dir+'microarray_data/tracker_stats.json')
        data_processor = Microarray_data_processing()
        print(f'processing {len(microarray_ids)} studies')
        valid_microarray_ids = download_experiments_microarray(data_processor,microarray_ids, downloads_folder, saved_tracker, download_raw=True, scan=False,output_folder=processed_folder)

        plot_tracker_results(f"{root_storage_dir}{scan_folder}tracker_stats.json", output_dir= scan_folder)

        combined,map = combine_files_microarray(processed_folder, "RMA_Microarray_Combined.csv", f"{root_storage_dir}final_data",combination_method='max',combine_genes=True)
    # raise ValueError('DONE')
    # combined = pd.read_csv(f'{root_storage_dir}final_data/RMA_Microarray_Combined.csv')

    # plot_study_distributions_seaborn(processed_folder, "new_storage/new_plots/intensity/intensity_plot_matplot_test.svg")
    # plot_study_distributions_incremental(processed_folder, "new_storage/new_plots/intensity/intensity_plot_incremental")
    # raise ValueError("DONE")



    # print("\n--- STARTING RNA-SEQ SEARCH ---")
    # tracker_loc = f"{root_storage_dir}rnaseq_data/RNA_tracker_stats_temp.json"
    # RNA_tracker = RNASeq_tracker.load_from_json(tracker_loc)
    def read_id(path):
        with open(path, 'r') as f:
                return (f.read().strip())
    if run_rna_seq:
        print('Running RNA SEQ')
        file_tracker_loc = f"{root_storage_dir}rnaseq_data/file_tracker/"
        
        # Load your IDs
        rnaseq_ids: list[str] = eval(read_id('RNA_seq_ids.txt'))
        # query_ids = search_geo_accessions(RNASEQ_QUERY, max_results=200000, filter_organism="Arabidopsis thaliana")
        
        # --- BATCH CONFIGURATION ---
        BATCH_SIZE:int = int(args.batch_size)
        RNA_tracker = FileTracker(file_tracker_loc)

        if args.array_index is not None:
            # --- PARALLEL BATCH MODE ---
            # Logic: Array Index 0 processes IDs 0-4, Index 1 processes 5-9, etc.
            
            start_idx:int = int(args.array_index) * BATCH_SIZE
            end_idx:int = start_idx + BATCH_SIZE
            
            # Use rnaseq_ids (from file) or query_ids (from search) depending on your goal
            target_list = rnaseq_ids 

            if start_idx < len(target_list):
                # Slice the list to get the batch
                current_batch = target_list[start_idx : end_idx]
                
                print(f"--- ARRAY JOB #{args.array_index} ---")
                print(f"Processing batch of {len(current_batch)} studies")
                print(f"Range: {start_idx} to {end_idx}")
                print(f"IDs: {current_batch}")
                print(f'tracker dir:{RNA_tracker.tracker_dir}')

                # Pass the BATCH list to the function
                download_experiments_RNA_seq_nf_core(
                    gse_list=current_batch,
                    root_storage_dir=root_storage_dir, 
                    output_dir=f"{root_storage_dir}rnaseq_data",
                    tracker=RNA_tracker, 
                    download_raw=True, 
                    scan=False,
                    run_and_delete=True, # Enable deletion to save space after batch
                    batch_size=BATCH_SIZE, # Should match slice size
                    debug=True
                )
            else:
                print(f"Index {args.array_index} (Start ID {start_idx}) is out of bounds for {len(target_list)} studies.")
        
        else:
            raise ValueError('NO array index')
            # --- SERIAL MODE ---
            # Pass the FULL list, the function handles the batching loop internally
            # print(f"--- SERIAL JOB: Processing {len(query_ids)} studies in batches of {BATCH_SIZE} ---")
            
            # download_experiments_RNA_seq_nf_core(
            #     gse_list=query_ids,
            #     root_storage_dir=root_storage_dir, 
            #     output_dir=f"{root_storage_dir}rnaseq_data",
            #     tracker=RNA_tracker, 
            #     download_raw=False, # Usually False for serial scanning? Set True if you want to download.
            #     scan=True,
            #     run_and_delete=True,
            #     batch_size=BATCH_SIZE
            # )

        # 2. Make the plots (Optional: might want to skip this in Array jobs to avoid write conflicts)
        if args.array_index is None:
            RNA_tracker.get_pie_charts()
            RNA_tracker.produce_study_dis()
            RNA_tracker.produce_platform_dis()
            
        print("\nDone!")