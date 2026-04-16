import argparse
import os
import sys

import dotenv
from Bio import Entrez

# Ensure we can import the local modules
module_dir = "./"
sys.path.append(module_dir)

from src.data_importing.get_count_matrices import download_processed_counts  # noqa: E402
from src.data_importing.helpers.file_tracker import FileTracker  # noqa: E402
from src.data_importing.helpers.helpers import combine_files_microarray, plot_tracker_results  # noqa: E402

# Import the processing function from your other file
from src.data_importing.microarray_data_processing import Microarray_data_processing, Microarray_tracker, download_experiments_microarray  # noqa: E402
from src.data_importing.RNA_seq_processing_batch import download_experiments_RNA_seq_nf_core  # noqa: E402

# --- CONFIGURATION ---
# dotenv.load_dotenv()
dotenv.load_dotenv(os.path.join(os.path.dirname(__file__), "../../.env"))
Entrez.email = os.getenv("EMAIL")
Entrez.api_key = os.getenv("NCBI_APY_KEY")


def save_list_to_txt(data_list, filename):
    """Saves a list to a text file, one item per line."""
    with open(filename, "w") as f:
        for item in data_list:
            f.write(f"{item}\n")
    print(f"Saved {len(data_list)} IDs to {filename}")


def load_list_from_txt(filename):
    """Reads a text file into a list, stripping newlines."""
    if not os.path.exists(filename):
        return []

    with open(filename) as f:
        # .strip() removes the \n character from the end
        return [line.strip() for line in f]


MICROARRAY_QUERY = '"Arabidopsis thaliana"[Organism] AND "Expression profiling by array"[DataSet Type] AND "GSE"[Entry Type]'  # cel"[Supplementary Files]'
RNASEQ_QUERY = '"Arabidopsis thaliana"[Organism] AND "Expression profiling by high throughput sequencing"[DataSet Type] AND "GSE"[Entry Type]'
RNASEQ_QUERY = (
    '"Arabidopsis thaliana"[Organism] AND '
    '"Expression profiling by high throughput sequencing"[DataSet Type] AND '
    '"GSE"[Entry Type] '
    "NOT ("
    '"Non-coding RNA profiling by high throughput sequencing"[DataSet Type] OR '
    '"Methylation profiling by high throughput sequencing"[DataSet Type] OR '
    '"Genome binding/occupancy profiling by high throughput sequencing"[DataSet Type] OR '
    '"Genome variation profiling by high throughput sequencing"[DataSet Type] OR '
    '"Expression profiling by array"[DataSet Type] OR '
    '"Other"[DataSet Type]'
    ")"
)

# Usage
# ids = ['GSE77815', 'GSE44053'] # Your list
# for i in ids:
#     download_processed_counts(i, "./count_data_folder")
# --- MAIN EXECUTION ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--out_dir", help="output_dir", default="./new_storage/")
    parser.add_argument("-b", "--batch_size", help="output_dir", default=100, type=int)
    parser.add_argument("--array_index", type=int, default=10, help="SLURM Array Task ID")
    parser.add_argument("--ma", action="store_true", default=False)
    parser.add_argument("--rna", action="store_true", default=False)
    parser.add_argument("--container", action="store_true", default=False)
    args = parser.parse_args()

    root_storage_dir = args.out_dir
    # 1. Microarray Data
    run_microarray: bool = args.ma
    run_rna_seq: bool = args.rna

    # Initialize Tracker
    if run_microarray:
        print("--- STARTING MICROARRAY SEARCH ---")
        ma = 20000
        scan_folder = f".{root_storage_dir}microarray_scan/"
        processed_folder = f"{root_storage_dir}processed_microarray_data/"
        downloads_folder = f"{root_storage_dir}microarray_data/"
        # microarray_ids = search_geo_accessions(MICROARRAY_QUERY, max_results=ma)
        filename = "./study_ids/microarray_ids.txt"
        # save_list_to_txt(microarray_ids,filename)
        microarray_ids = load_list_from_txt(filename)

        # Pass tracker to the download function
        # ma_tracker.sync_with_filesystem(downloads_folder,processed_folder)
        saved_tracker = Microarray_tracker.load_from_json(root_storage_dir + "microarray_data/tracker_stats.json")
        data_processor = Microarray_data_processing()
        print(f"processing {len(microarray_ids)} studies")
        valid_microarray_ids = download_experiments_microarray(data_processor, microarray_ids, downloads_folder, saved_tracker, download_raw=True, scan=False, output_folder=processed_folder)

        plot_tracker_results(f"{root_storage_dir}{scan_folder}tracker_stats.json", output_dir=scan_folder)

        combined, map = combine_files_microarray(processed_folder, "RMA_Microarray_Combined.csv", f"{root_storage_dir}final_data", combination_method="max", combine_genes=True)
    # raise ValueError('DONE')
    # combined = pd.read_csv(f'{root_storage_dir}final_data/RMA_Microarray_Combined.csv')

    # plot_study_distributions_seaborn(processed_folder, "new_storage/new_plots/intensity/intensity_plot_matplot_test.svg")
    # plot_study_distributions_incremental(processed_folder, "new_storage/new_plots/intensity/intensity_plot_incremental")
    # raise ValueError("DONE")

    # print("\n--- STARTING RNA-SEQ SEARCH ---")
    # tracker_loc = f"{root_storage_dir}rnaseq_data/RNA_tracker_stats_temp.json"
    # RNA_tracker = RNASeq_tracker.load_from_json(tracker_loc)
    def read_id(path):
        with open(path) as f:
            return f.read().strip()

    def save_id_list(data_list, path):
        with open(path, "w") as f:
            # This writes "GSE123\nGSE456\nGSE789" to the file
            f.write(",".join(data_list))

    if run_rna_seq:
        print("Running RNA SEQ")
        file_tracker_loc = f"{root_storage_dir}rnaseq_data/file_tracker/"

        # Load your IDs
        rnaseq_ids: list[str] = eval("['" + read_id("./study_ids/RNA_seq_ids.txt").replace(",", "','") + "']")
        # query_ids = search_geo_accessions(RNASEQ_QUERY, max_results=200000, filter_organism="Arabidopsis thaliana")

        # TODO: change this to that we filter out all files that are being tried or tried???

        # --- BATCH CONFIGURATION ---
        BATCH_SIZE: int = int(args.batch_size)
        RNA_tracker = FileTracker(file_tracker_loc)

        if args.array_index is not None:
            # --- PARALLEL BATCH MODE ---
            # Logic: Array Index 0 processes IDs 0-4, Index 1 processes 5-9, etc.

            start_idx: int = int(args.array_index) * BATCH_SIZE
            end_idx: int = start_idx + BATCH_SIZE

            # Use rnaseq_ids (from file) or query_ids (from search) depending on your goal
            target_list = rnaseq_ids
            if False:
                ret = []
                for id in target_list:
                    ret.append(download_processed_counts(id, "./temp/"))
            if start_idx < len(target_list):
                # Slice the list to get the batch
                current_batch = target_list[start_idx:end_idx]
                # for id in current_batch:
                #     RNA_tracker.set_status(id,0)
                print(f"--- ARRAY JOB #{args.array_index} ---")
                print(f"Processing batch of {len(current_batch)} studies")
                print(f"Range: {start_idx} to {end_idx}")
                print(f"IDs: {current_batch}")
                print(f"tracker dir:{RNA_tracker.tracker_dir}")

                # Pass the BATCH list to the function
                download_experiments_RNA_seq_nf_core(
                    gse_list=current_batch,
                    root_storage_dir=root_storage_dir,
                    output_dir=f"{root_storage_dir}rnaseq_data",
                    tracker=RNA_tracker,
                    download_raw=True,
                    metadata_only=False,
                    run_and_delete=True,  # Enable deletion to save space after batch
                    batch_size=BATCH_SIZE,  # Should match slice size
                    debug=False,
                    container=args.container,
                )
            else:
                print(f"Index {args.array_index} (Start ID {start_idx}) is out of bounds for {len(target_list)} studies.")

        else:
            msg = "NO array index"
            raise ValueError(msg)

        # 2. Make the plots (Optional: might want to skip this in Array jobs to avoid write conflicts)
        if args.array_index is None:
            RNA_tracker.get_pie_charts()
            RNA_tracker.produce_study_dis()
            RNA_tracker.produce_platform_dis()

        print("\nDone!")
