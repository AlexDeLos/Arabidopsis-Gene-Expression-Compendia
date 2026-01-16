import os
import json
import pandas as pd
import numpy as np
import sys
import shutil
from tqdm import tqdm
import json
import GEOparse

module_dir = './'
sys.path.append(module_dir)

# Import the helper to process metadata (Same as Microarray)
from src.data_importing.helpers import process_metadata
from src.data_importing.download_helper import clean_files,check_metadata_for_counts

class RNASeq_tracker:
    def __init__(self) -> None:
        self.platform_counts: dict = {}
        self.totals: dict = {
            'total_studies_seen': 0,
            'total_sample_seen': 0,
            'total_samples_used': 0, 
            'total_studies_used': 0
        }
        # State tracking sets (Mutually Exclusive)
        self.states: dict = {
            'ignore': set(),
            'downloaded': set(),
            'processed': set()
        }
        self.study_tracker: dict = {}

    def update_platform(self, platform: str, samples: int, has_counts: bool):
        if platform not in self.platform_counts:
            self.platform_counts[platform] = {
                'studies_seen': 0, 'samples_seen': 0,
                'studies_with_counts': 0, 'samples_with_counts': 0
            }
        
        self.totals['total_studies_seen'] += 1
        self.totals['total_sample_seen'] += samples
        self.platform_counts[platform]['studies_seen'] += 1
        self.platform_counts[platform]['samples_seen'] += samples

        if has_counts:
            self.totals['total_studies_used'] += 1
            self.totals['total_samples_used'] += samples
            self.platform_counts[platform]['studies_with_counts'] += 1
            self.platform_counts[platform]['samples_with_counts'] += samples

    # --- STATE MANAGEMENT ---

    def mark_ignore(self, gse_id):
        self.states['ignore'].add(gse_id)
        self.states['downloaded'].discard(gse_id)
        self.states['processed'].discard(gse_id)

    def mark_downloaded(self, gse_id):
        self.states['downloaded'].add(gse_id)
        self.states['ignore'].discard(gse_id)
        self.states['processed'].discard(gse_id)

    def mark_processed(self, gse_id):
        self.states['processed'].add(gse_id)
        self.states['downloaded'].discard(gse_id)
        self.states['ignore'].discard(gse_id)

    def is_ignored(self, gse_id):
        return gse_id in self.states['ignore']

    def is_processed(self, gse_id):
        return gse_id in self.states['processed']
    
    def is_downloaded(self, gse_id):
        return gse_id in self.states['downloaded']

    def save_to_json(self, filename="rnaseq_tracker_results.json"):
        serializable_states = {k: list(v) for k, v in self.states.items()}
        data = {
            "totals": self.totals,
            "platform_counts": self.platform_counts,
            "states": serializable_states
        }
        with open(filename, 'w') as f:
            json.dump(data, f, indent=4)
    
    @classmethod
    def load_from_json(cls, filename="rnaseq_tracker_results.json"):
        if not os.path.exists(filename):
            print(f"File {filename} not found. Returning empty tracker.")
            return cls()
        with open(filename, 'r') as f:
            data = json.load(f)
        tracker = cls()
        tracker.totals = data.get("totals", tracker.totals)
        tracker.platform_counts = data.get("platform_counts", {})
        loaded_states = data.get("states", {})
        tracker.states = {
            'ignore': set(loaded_states.get('ignore', [])),
            'downloaded': set(loaded_states.get('downloaded', [])),
            'processed': set(loaded_states.get('processed', []))
        }
        return tracker

    def print_summary(self):
        print("\n" + "="*40)
        print("       RNA-SEQ PROCESSING SUMMARY")
        print("="*40)
        print(f"Total Studies: {self.totals['total_studies_seen']} Seen | {self.totals['total_studies_used']} Has Counts")
        print(f"Ignored: {len(self.states['ignore'])} | Downloaded: {len(self.states['downloaded'])} | Processed: {len(self.states['processed'])}")
        print("="*40 + "\n")

class RNASeq_processor:
    def __init__(self):
        pass

    def extract_tar_files(self, folder):
        """
        Finds any .tar or .tar.gz files in the folder and extracts them.
        """
        for f in os.listdir(folder):
            if f.endswith(('.tar', '.tar.gz')):
                full_path = os.path.join(folder, f)
                print(f"    > Extracting archive: {f}...")
                try:
                    with tarfile.open(full_path, "r:*") as tar:
                        tar.extractall(path=folder)
                except Exception as e:
                    print(f"    ! Failed to extract {f}: {e}")

    def aggregate_individual_files(self, folder, gse):
        """
        If no master matrix is found, try to find 1 file per sample (GSM)
        and merge them into a matrix.
        """
        print("    > Attempting to aggregate individual sample files...")
        
        sample_dfs = []
        found_samples = []
        
        # Get all files in the folder
        all_files = [f for f in os.listdir(folder) if os.path.isfile(os.path.join(folder, f))]
        
        # Iterate over expected GSMs
        for gsm_id in gse.gsms.keys():
            # Look for a file that contains this GSM ID
            # Priority: Exact match > Contains ID
            match = None
            for f in all_files:
                if gsm_id in f and f.lower().endswith(('.txt', '.tsv', '.csv', '.tab', '.txt.gz', '.tsv.gz')):
                    match = f
                    break
            
            if match:
                file_path = os.path.join(folder, match)
                try:
                    # Read individual file. Usually no header, or header with "Gene, Count"
                    # We assume Col 0 is Gene ID, Col 1 is Count.
                    # 'header=None' might be needed if there are no headers, but 'infer' usually works
                    df = pd.read_csv(file_path, sep=None, engine='python', index_col=0)
                    
                    # If the dataframe has multiple columns, try to guess which is the count
                    if df.shape[1] > 1:
                        # Heuristic: Take the first numeric column
                        df = df.select_dtypes(include=[np.number])
                        if df.shape[1] > 0:
                            df = df.iloc[:, [0]] # Take first numeric col
                    
                    # Rename the column to the GSM ID
                    df.columns = [gsm_id]
                    
                    # Deduplicate indices (genes) if necessary
                    if df.index.duplicated().any():
                         df = df.groupby(df.index).mean()

                    sample_dfs.append(df)
                    found_samples.append(gsm_id)
                except Exception:
                    continue
        
        if not sample_dfs:
            return None

        print(f"    > Found individual files for {len(sample_dfs)} / {len(gse.gsms)} samples.")
        
        # Merge all into one dataframe (Outer join to keep all genes)
        master_df = pd.concat(sample_dfs, axis=1, join='outer')
        return master_df

    def process_rnaseq_metadata_and_counts(self, experiment_id, raw_data_folder, output_folder, gse) -> int:
        """
        1. Extracts Metadata.
        2. Extracts TAR files (if any).
        3. Tries to find a MASTER matrix.
        4. If failing, tries to AGGREGATE individual files.
        5. Saves result.
        """
        print(f"\n[RNA-Seq Pipeline] Processing {experiment_id}...")
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        # 1. Process Metadata
        print("  - Processing metadata...")
        try:
            for gsm_name, gsm in gse.gsms.items():
                process_metadata(experiment_id, gse, gsm, save_path=output_folder)
        except Exception as e:
            print(f"    > Metadata processing failed: {e}")

        # 2. Extract any TAR/GZ archives downloaded
        self.extract_tar_files(raw_data_folder)

        # 3. Strategy A: Search for a Single Master Matrix
        print("  - Strategy A: Searching for a single master count matrix...")
        candidates = []
        for root, dirs, files in os.walk(raw_data_folder):
            for f in files:
                if f.lower().endswith(('.csv', '.tsv', '.txt', '.tab', '.xlsx', '.csv.gz')):
                    # Exclude small files (likely READMEs)
                    if os.path.getsize(os.path.join(root, f)) > 1000: 
                        candidates.append(os.path.join(root, f))
        
        valid_matrix = None
        expected_samples = set(gse.gsms.keys())

        for file_path in candidates:
            try:
                if file_path.endswith('.xlsx'):
                    df = pd.read_excel(file_path, index_col=0)
                else:
                    df = pd.read_csv(file_path, sep=None, engine='python', index_col=0)
                
                columns_clean = [str(c).strip() for c in df.columns]
                matches = set(columns_clean).intersection(expected_samples)
                
                # If we match > 2 samples or 50% of expected
                if len(matches) > 2 or len(matches) >= len(expected_samples) * 0.5:
                    valid_matrix = df
                    print(f"    > Found master matrix: {os.path.basename(file_path)}")
                    break
            except Exception:
                continue

        # 4. Strategy B: Aggregate Individual Files
        if valid_matrix is None:
            print("    > No master matrix found. Trying Strategy B (Aggregate individual files)...")
            valid_matrix = self.aggregate_individual_files(raw_data_folder, gse)

        if valid_matrix is None:
            raise ValueError("Failed to construct a count matrix (Neither master file nor individual sample files found).")

        # 5. Save Cleaned Matrix
        valid_matrix.index.name = "GeneID"
        output_file = os.path.join(output_folder, f"{experiment_id}_counts.csv")
        valid_matrix.to_csv(output_file)
        print(f"  - SUCCESS: Saved count matrix to {output_file}")
        
        return 0
    
def download_experiments_RNA_seq(gse_list, download_dir, tracker, download_raw=True, scan=True, output_folder='./rnaseq_processed_data/'):
    """
    Downloads RNA-Seq studies that provide Count Matrices (not FASTQ).
    """
    processor = RNASeq_processor()
    if not os.path.exists(download_dir):
        os.makedirs(download_dir)

    valid_gse_ids = []
    
    # Path where the tracker will be auto-saved
    tracker_save_path = os.path.join(download_dir, "../rnaseq_tracker_stats.json")

    for gse_id in tqdm(gse_list, desc='Processing RNA-Seq (Counts)', unit='study'):
        
        # 1. CHECK TRACKER
        if tracker.is_ignored(gse_id):
            continue
        if tracker.is_processed(gse_id):
            valid_gse_ids.append(gse_id)
            continue
            
        soft_file_path = os.path.join(download_dir, f'{gse_id}_family.soft.gz')
        exp_folder = os.path.join(download_dir, gse_id)
        
        try:
            gse = None
            
            # --- PHASE 1: DOWNLOAD & CHECK METADATA ---
            if not tracker.is_downloaded(gse_id):
                try:
                    gse = GEOparse.get_GEO(geo=gse_id, destdir=download_dir, silent=True)
                except Exception as e:
                    tqdm.write(f"  [{gse_id}] Metadata download failed: {e}")
                    tracker.mark_ignore(gse_id)
                    tracker.save_to_json(tracker_save_path)
                    clean_files([soft_file_path])
                    continue
                
                # Check Platform / Samples
                try:
                    platform = gse.metadata.get('platform_id', ['Unknown'])[0]
                    num_samples = len(gse.gsms)
                except:
                    platform = 'Unknown'
                    num_samples = 0
                
                # Check for Count Files
                has_counts = check_metadata_for_counts(gse) if download_raw else False

                # SCAN MODE
                if scan:
                    tracker.update_platform(platform, num_samples, has_counts=has_counts)
                    if has_counts:
                        valid_gse_ids.append(gse_id)
                    clean_files([soft_file_path])
                    continue

                # NORMAL MODE
                if download_raw:
                    if not has_counts:
                        tqdm.write(f'  [{gse_id}] SKIPPING: No Count Matrix found in metadata.')
                        tracker.update_platform(platform, num_samples, has_counts=False)
                        tracker.mark_ignore(gse_id)
                        tracker.save_to_json(tracker_save_path)
                        clean_files([soft_file_path])
                        continue
                    
                    # Download Supplementary Files
                    if not os.path.exists(exp_folder):
                        os.makedirs(exp_folder)
                    
                    # We download supplementary files. 
                    # Note: We disable SRA download, we only want attached txt/csv
                    gse.download_supplementary_files(directory=exp_folder, download_sra=False)

                    # Check if anything downloaded
                    if len(os.listdir(exp_folder)) > 0:
                        valid_gse_ids.append(gse_id)
                        tracker.update_platform(platform, num_samples, has_counts=True)
                        
                        tracker.mark_downloaded(gse_id)
                        tracker.save_to_json(tracker_save_path)
                    else:
                        tqdm.write(f'  [{gse_id}] WARNING: Metadata claimed counts, but no files downloaded.')
                        tracker.update_platform(platform, num_samples, has_counts=False)
                        tracker.mark_ignore(gse_id)
                        tracker.save_to_json(tracker_save_path)
                        clean_files([exp_folder, soft_file_path])
                        continue
            
            # --- PHASE 2: PROCESSING ---
            if download_raw and not scan:
                if gse is None:
                    try:
                        gse = GEOparse.get_GEO(filepath=soft_file_path, destdir=download_dir, silent=True)
                    except:
                        gse = GEOparse.get_GEO(geo=gse_id, destdir=download_dir, silent=True)

                tqdm.write(f'  [{gse_id}] Processing Count Matrix...')
                try:
                    # Call the new processor
                    processor.process_rnaseq_metadata_and_counts(gse_id, exp_folder, f'{output_folder}/{gse_id}', gse)
                    
                    tracker.mark_processed(gse_id)
                    tracker.save_to_json(tracker_save_path)
                    
                    # Optional: Clean up the downloaded raw files after processing to save space?
                    # clean_files([exp_folder]) 
                    
                except Exception as e_proc:
                    tqdm.write(f'  [{gse_id}] PROCESSING FAILED: {e_proc}')
                    tracker.mark_ignore(gse_id)
                    tracker.save_to_json(tracker_save_path)
                    clean_files([soft_file_path, exp_folder])
                    continue

        except Exception as e:
            tqdm.write(f'  [{gse_id}] FAILED: {e}')
            tracker.mark_ignore(gse_id)
            tracker.save_to_json(tracker_save_path)
            clean_files([soft_file_path, exp_folder])

    return valid_gse_ids