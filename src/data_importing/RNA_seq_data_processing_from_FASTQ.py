import os
import json
import pandas as pd
import subprocess
import shutil
from Bio import Entrez
from tqdm import tqdm
import GEOparse
import sys
module_dir = './'
sys.path.append(module_dir)
from src.data_importing.download_helper import check_metadata_for_sra

# Standard Nextflow command template
# Adjust 'profile' to 'singularity' or 'conda' if needed
NF_CORE_COMMAND = """
nextflow run nf-core/rnaseq \
    -profile docker \
    --input {samplesheet} \
    --outdir {outdir} \
    --skip_alignment \
    --pseudo_aligner salmon \
    --remove_ribo_rna \
    -resume
"""
# Note: I selected --skip_alignment and --pseudo_aligner salmon for speed (TPM generation),
# as suggested by modern protocols. Remove these flags if you want full STAR alignment.

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

    def update_platform(self, platform: str, samples: int, has_raw: bool):
        if platform not in self.platform_counts:
            self.platform_counts[platform] = {
                'studies_seen': 0, 'samples_seen': 0,
                'studies_with_raw': 0, 'samples_with_raw': 0
            }
        
        self.totals['total_studies_seen'] += 1
        self.totals['total_sample_seen'] += samples
        self.platform_counts[platform]['studies_seen'] += 1
        self.platform_counts[platform]['samples_seen'] += samples

        if has_raw:
            self.totals['total_studies_used'] += 1
            self.totals['total_samples_used'] += samples
            self.platform_counts[platform]['studies_with_raw'] += 1
            self.platform_counts[platform]['samples_with_raw'] += samples

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

    def sync_with_filesystem(self, download_dir: str, processed_dir: str):
        """
        Syncs state based on files on disk.
        """
        print(f"Syncing RNA-Seq tracker with files...")
        
        # 1. Identify Processed (Look for Salmon/Star output folders or summary CSVs)
        # Assuming nf-core output structure: processed_dir/GSE.../salmon/
        actual_processed_ids = set()
        if os.path.exists(processed_dir):
            for gse_id in os.listdir(processed_dir):
                # Check if Salmon quantification exists
                salmon_path = os.path.join(processed_dir, gse_id, "salmon")
                # Or check for the final combined matrix if you generated one
                if os.path.exists(salmon_path):
                    actual_processed_ids.add(gse_id)

        # 2. Identify Downloaded (Look for FASTQ files)
        actual_downloaded_ids = set()
        if os.path.exists(download_dir):
            for d in os.listdir(download_dir):
                full_path = os.path.join(download_dir, d)
                if d.startswith("GSE") and os.path.isdir(full_path):
                    # Check if folder is not empty
                    if len(os.listdir(full_path)) > 0:
                        actual_downloaded_ids.add(d)

        # 3. Update Tracker
        changes = 0
        for gse in actual_processed_ids:
            if not self.is_processed(gse):
                self.mark_processed(gse)
                changes += 1
        
        for gse in actual_downloaded_ids:
            if gse not in actual_processed_ids:
                if not self.is_downloaded(gse):
                    self.mark_downloaded(gse)
                    changes += 1
        
        print(f"  - RNA-Seq Sync complete. Updated {changes} records.")

    def print_summary(self):
        print("\n" + "="*40)
        print("       RNA-SEQ PROCESSING SUMMARY")
        print("="*40)
        print(f"Total Studies: {self.totals['total_studies_seen']} Seen | {self.totals['total_studies_used']} Has Raw (SRA)")
        print(f"Ignored: {len(self.states['ignore'])} | Downloaded: {len(self.states['downloaded'])} | Processed: {len(self.states['processed'])}")
        print("="*40 + "\n")


class RNASeq_processor:
    def __init__(self):
        # check for fasterq-dump
        if not shutil.which("fasterq-dump"):
            print("WARNING: 'fasterq-dump' (SRA Toolkit) not found in PATH. Downloading will fail.")
        # check for nextflow
        if not shutil.which("nextflow"):
            print("WARNING: 'nextflow' not found in PATH. Processing will fail.")

    def get_srr_ids(self, gsm_id):
        """Maps a GEO GSM ID to SRA Run IDs (SRR)."""
        try:
            # 1. Search SRA for the GSM ID
            handle = Entrez.esearch(db="sra", term=gsm_id)
            record = Entrez.read(handle)
            handle.close()
            
            if not record['IdList']:
                return []
            
            # 2. Fetch the Run Info
            sra_ids = record['IdList']
            handle = Entrez.esummary(db="sra", id=",".join(sra_ids))
            summaries = Entrez.read(handle)
            handle.close()
            
            run_ids = []
            for summary in summaries:
                # Parse the weird XML-like string in 'Runs' field
                # Format usually: '<Run acc="SRR12345" ...'
                runs = summary.get('Runs', '')
                import re
                matches = re.findall(r'acc="([A-Z0-9]+)"', runs)
                run_ids.extend(matches)
            
            return list(set(run_ids))
        except Exception as e:
            print(f"    ! Error fetching SRA IDs for {gsm_id}: {e}")
            return []

    def download_fastq_from_sra(self, gse, output_folder):
        """
        Downloads all FASTQ files for a GSE study using fasterq-dump.
        """
        print(f"  - Resolving SRA IDs for {len(gse.gsms)} samples...")
        
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        sra_map = {} # GSM -> [SRR1, SRR2]
        
        for gsm_name, gsm in tqdm(gse.gsms.items(), desc="Fetching Metadata", leave=False):
            srrs = self.get_srr_ids(gsm_name)
            if srrs:
                sra_map[gsm_name] = srrs
        
        if not sra_map:
            raise ValueError("No linked SRA entries found for this study.")

        print(f"  - Downloading FASTQ files for {len(sra_map)} samples...")
        
        success_count = 0
        for gsm, srrs in tqdm(sra_map.items(), desc="Downloading FASTQ"):
            for srr in srrs:
                # Check if file already exists
                # fasterq-dump output is usually SRRxxxx.fastq or SRRxxxx_1.fastq
                expected_1 = os.path.join(output_folder, f"{srr}.fastq")
                expected_2 = os.path.join(output_folder, f"{srr}_1.fastq")
                
                if os.path.exists(expected_1) or os.path.exists(expected_2):
                    continue

                # Run fasterq-dump
                # --split-files ensures paired-end data is split into _1 and _2
                cmd = ["fasterq-dump", "--split-files", "--outdir", output_folder, "--threads", "4", srr]
                try:
                    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                    success_count += 1
                except subprocess.CalledProcessError:
                    print(f"    ! Failed to download {srr}")
        
        return success_count

    def generate_nfcore_samplesheet(self, fastq_folder, output_csv):
        """
        Scans a folder for FASTQ files and generates a samplesheet.csv 
        compatible with nf-core/rnaseq.
        Columns: sample,fastq_1,fastq_2,strandedness
        """
        samples = {}
        
        for f in os.listdir(fastq_folder):
            if f.endswith(".fastq") or f.endswith(".fq") or f.endswith(".gz"):
                # Guess sample name (SRRxxxx)
                parts = f.split('_')
                sample = parts[0].split('.')[0] # Get SRR ID
                
                if sample not in samples:
                    samples[sample] = {'1': '', '2': ''}
                
                # Assign read 1 or 2
                if "_1" in f:
                    samples[sample]['1'] = os.path.join(fastq_folder, f)
                elif "_2" in f:
                    samples[sample]['2'] = os.path.join(fastq_folder, f)
                else:
                    # Single end
                    samples[sample]['1'] = os.path.join(fastq_folder, f)

        rows = []
        for sample, files in samples.items():
            fq1 = files['1']
            fq2 = files['2']
            if fq1:
                # Auto-detect strandedness is 'auto' in nf-core
                rows.append({'sample': sample, 'fastq_1': fq1, 'fastq_2': fq2, 'strandedness': 'auto'})

        if not rows:
            return False

        df = pd.DataFrame(rows)
        # Ensure fastq_2 is empty string if NaN
        df['fastq_2'] = df['fastq_2'].fillna('')
        df.to_csv(output_csv, index=False)
        return True

    def run_nfcore_rnaseq(self, gse_id, fastq_folder, output_folder):
        """
        Runs the nf-core/rnaseq pipeline.
        """
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
            
        samplesheet_path = os.path.join(output_folder, "samplesheet.csv")
        
        # 1. Generate Samplesheet
        if not self.generate_nfcore_samplesheet(fastq_folder, samplesheet_path):
            raise ValueError("Could not generate samplesheet (no FASTQ found?)")
            
        # 2. Build Command
        # Using formatted string defined at top of file
        cmd_str = NF_CORE_COMMAND.format(
            samplesheet=samplesheet_path,
            outdir=output_folder
        )
        
        # 3. Execute
        print(f"  - Launching nf-core/rnaseq for {gse_id}...")
        try:
            # We use shell=True to handle the multiline string and pipes if needed
            subprocess.run(cmd_str, shell=True, check=True)
            print("  - Pipeline finished successfully.")
        except subprocess.CalledProcessError as e:
            print(f"  - Pipeline FAILED: {e}")
            raise e
        

import os
import shutil # Required for folder deletion
from tqdm import tqdm

def download_experiments_RNA_seq(gse_list, output_dir, tracker, download_raw=True, scan=False):
    """
    Downloads and Processes RNA-Seq experiments.
    Implements 'Stream-Process-Delete' to minimize storage footprint.
    """
    # from src.data_importing.RNA_seq_data_processing import RNASeq_processor
    
    # Initialize Processor
    processor = RNASeq_processor()
    
    # Tracker File Path
    tracker_save_path = os.path.join(output_dir, "../rnaseq_tracker_stats.json")
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    valid_gse_ids = []

    for gse_id in tqdm(gse_list, desc="Processing RNA-Seq", unit="study"):
        
        # 1. CHECK TRACKER
        if tracker.is_ignored(gse_id):
            continue
        if tracker.is_processed(gse_id):
            valid_gse_ids.append(gse_id)
            continue
            
        soft_file_path = os.path.join(output_dir, f"{gse_id}_family.soft.gz")
        exp_folder = os.path.join(output_dir, gse_id) 
        fastq_folder = os.path.join(output_dir, "fastq_storage", gse_id) 
        results_folder = os.path.join(output_dir, "processed_rnaseq", gse_id) 

        try:
            # --- METADATA ---
            try:
                gse = GEOparse.get_GEO(geo=gse_id, destdir=output_dir, silent=True)
            except Exception as e:
                tqdm.write(f"  [{gse_id}] Metadata failed: {e}")
                tracker.mark_ignore(gse_id)
                if os.path.exists(soft_file_path): os.remove(soft_file_path)
                continue

            platform = gse.metadata.get('platform_id', ['Unknown'])[0]
            num_samples = len(gse.gsms)
            
            # Check for SRA
            has_raw_sra = check_metadata_for_sra(gse)

            # --- SCAN MODE ---
            if scan:
                tracker.update_platform(platform, num_samples, has_raw=has_raw_sra)
                if has_raw_sra:
                    valid_gse_ids.append(gse_id)
                
                # Cleanup Metadata
                if os.path.exists(soft_file_path): os.remove(soft_file_path)
                continue

            # --- NORMAL MODE ---
            if download_raw:
                if not has_raw_sra:
                    tqdm.write(f"  [{gse_id}] SKIPPING: No SRA link.")
                    tracker.mark_ignore(gse_id)
                    tracker.save_to_json(tracker_save_path)
                    if os.path.exists(soft_file_path): os.remove(soft_file_path)
                    continue
                
                tracker.update_platform(platform, num_samples, has_raw=True)
                
                # A. DOWNLOAD FASTQ (If not already done)
                # We download strictly for immediate processing
                if not tracker.is_downloaded(gse_id) and not os.path.exists(fastq_folder):
                    tqdm.write(f"  [{gse_id}] Downloading FASTQ from SRA...")
                    try:
                        processor.download_fastq_from_sra(gse, fastq_folder)
                        tracker.mark_downloaded(gse_id)
                        tracker.save_to_json(tracker_save_path)
                    except Exception as e:
                        tqdm.write(f"  [{gse_id}] DOWNLOAD FAILED: {e}")
                        tracker.mark_ignore(gse_id) 
                        tracker.save_to_json(tracker_save_path)
                        # Clean partial download
                        if os.path.exists(fastq_folder): shutil.rmtree(fastq_folder)
                        if os.path.exists(soft_file_path): os.remove(soft_file_path)
                        continue

                # B. PROCESS PIPELINE
                if not tracker.is_processed(gse_id):
                    tqdm.write(f"  [{gse_id}] Running nf-core pipeline...")
                    try:
                        processor.run_nfcore_rnaseq(gse_id, fastq_folder, results_folder)
                        
                        # --- C. CRITICAL: VERIFY & DELETE ---
                        # Check if pipeline generated the expected output folder (e.g. 'salmon' or 'star_salmon')
                        # and that it is not empty.
                        has_results = False
                        if os.path.exists(results_folder) and len(os.listdir(results_folder)) > 0:
                            has_results = True

                        if has_results:
                            tqdm.write(f"  [{gse_id}] SUCCESS. Deleting raw FASTQ to save space.")
                            
                            # 1. Delete Huge FASTQ Files
                            if os.path.exists(fastq_folder):
                                shutil.rmtree(fastq_folder)
                            
                            # 2. Update Tracker
                            tracker.mark_processed(gse_id)
                            tracker.save_to_json(tracker_save_path)
                            valid_gse_ids.append(gse_id)
                        else:
                            raise Exception("Pipeline finished but output folder is empty.")
                        
                    except Exception as e:
                        tqdm.write(f"  [{gse_id}] PIPELINE FAILED: {e}")
                        # If pipeline failed, we usually KEEP the FASTQ so we can debug/retry later
                        # without re-downloading. But if you are strict on space, uncomment below:
                        # if os.path.exists(fastq_folder): shutil.rmtree(fastq_folder)
                        continue
            
            # Cleanup Metadata (Normal Mode) to keep folder clean
            if os.path.exists(soft_file_path):
                os.remove(soft_file_path)

        except Exception as e:
            tqdm.write(f"  [{gse_id}] FAILED: {e}")
            tracker.mark_ignore(gse_id)
            if os.path.exists(soft_file_path): os.remove(soft_file_path)

    return valid_gse_ids