import os
import subprocess
import shutil
import csv
import pandas as pd
from Bio import Entrez
from tqdm import tqdm
import GEOparse
from GEOparse.GEOTypes import GSE, GSM, GPL, GDS
import sys
import time
import glob

module_dir = './'
sys.path.append(module_dir)
from src.data_importing.helpers.download_helper import check_metadata_for_sra_boolean
from src.constants import *
from src.data_importing.helpers.file_tracker import FileTracker


class RNASeq_processor:
    def __init__(self, threads=4, genome_index=None, gtf_annotation=None, profile='docker'):
        self.threads = str(threads)
        self.profile = profile
        
        # Verify Tools
        required = ['fasterq-dump', 'trimmomatic', 'hisat2', 'samtools', 'featureCounts']
        for tool in required:
            if not shutil.which(tool):
                print(f"WARNING: {tool} not found in PATH. Pipeline may fail.")

    def get_srr_ids(self, gsm_id):
        """Fetches SRR IDs for a GSM from NCBI."""
        try:
            handle = Entrez.esearch(db="sra", term=gsm_id)
            record = Entrez.read(handle)
            handle.close()
            if not record['IdList']: return []
            
            handle = Entrez.esummary(db="sra", id=",".join(record['IdList']))
            summaries = Entrez.read(handle)
            handle.close()
            
            run_ids = []
            import re
            for summary in summaries:
                run_ids.extend(re.findall(r'acc="([A-Z0-9]+)"', summary.get('Runs', '')))
            return list(set(run_ids))
        except:
            return []

    def download_fastq_faster(self, gse, output_folder, temp_files):
        """
        Downloads using fasterq-dump (fast) and immediately compresses using pigz (fast).
        """
        if not os.path.exists(output_folder): os.makedirs(output_folder)
        if not os.path.exists(temp_files): os.makedirs(temp_files)

        sra_map = {gsm: self.get_srr_ids(gsm) for gsm in gse.gsms.keys()}
        sra_map = {k:v for k,v in sra_map.items() if v}
        
        fasterq_path = os.path.expanduser("~/fasterq-dump")
        zip_cmd = "pigz" if shutil.which("pigz") else "gzip"

        for gsm, srrs in tqdm(sra_map.items(), desc="Downloading SRRs", leave=False):
            for srr in srrs:
                # 1. Check if ZIPPED file exists
                existing_gz = [f for f in os.listdir(output_folder) if f.startswith(srr) and f.endswith('.gz')]
                if existing_gz: 
                    continue

                # 2. Run fasterq-dump
                cmd = [fasterq_path if CLUSTER_RUN else "fasterq-dump", 
                       "--split-files", "--outdir", output_folder, 
                       "--temp", temp_files , "--threads", self.threads, srr]
                
                try:
                    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
                    
                    # 3. Compress
                    raw_fastqs = [f for f in os.listdir(output_folder) if f.startswith(srr) and f.endswith('.fastq')]
                    if not raw_fastqs:
                        continue

                    for fq in raw_fastqs:
                        full_path = os.path.join(output_folder, fq)
                        if zip_cmd == "pigz":
                            subprocess.run(["pigz", "-f", "-p", self.threads, full_path], check=True)
                        else:
                            subprocess.run(["gzip", "-f", full_path], check=True)
                                
                except subprocess.CalledProcessError as e:
                    print(f"\nCRITICAL ERROR downloading {srr}: {e}")
                    raise e
                
    def download_fastq(self, gse, output_folder, temp_files):
        """Downloads using fastq-dump with built-in compression."""
        if not os.path.exists(output_folder): os.makedirs(output_folder)
        if not os.path.exists(temp_files): os.makedirs(temp_files)

        sra_map = {gsm: self.get_srr_ids(gsm) for gsm in gse.gsms.keys()}
        sra_map = {k:v for k,v in sra_map.items() if v}
        
        fastq_dump_cmd = "fastq-dump"

        for gsm, srrs in tqdm(sra_map.items(), desc="Downloading SRRs", leave=False):
            for srr in srrs:
                existing_gz = [f for f in os.listdir(output_folder) if f.startswith(srr) and f.endswith('.gz')]
                if existing_gz: continue

                cmd = [fastq_dump_cmd, "--gzip", "--split-files", "--outdir", output_folder, srr]
                max_retries = 3
                success = False

                for attempt in range(max_retries):
                    try:
                        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
                        if any(f.startswith(srr) and f.endswith('.gz') for f in os.listdir(output_folder)):
                            success = True
                            break
                    except subprocess.CalledProcessError:
                        print(f"Retrying {srr} ({attempt+1}/{max_retries})...")
                        # Cleanup partials
                        for pf in [f for f in os.listdir(output_folder) if f.startswith(srr)]:
                            try: os.remove(os.path.join(output_folder, pf))
                            except: pass
                        time.sleep(10)

                if not success:
                    print(f"\nFailed to download {srr} after retries.")

    def get_samplesheet_rows(self, gse_id, fastq_folder):
        """
        Generates samplesheet rows for a single study.
        Returns: List of rows [sample, fq1, fq2, strandedness]
        """
        rows = []
        samples = {}
        # Identify samples and pair files
        files = [f for f in os.listdir(fastq_folder) if f.endswith('.fastq') or f.endswith('.fq') or f.endswith('.gz')]
        
        for f in files:
            path = os.path.abspath(os.path.join(fastq_folder, f)) # Use absolute paths for batching
            if '_1' in f:
                srr = f.split('_1')[0]
                samples.setdefault(srr, {'1': None, '2': None})['1'] = path
            elif '_2' in f:
                srr = f.split('_2')[0]
                samples.setdefault(srr, {'1': None, '2': None})['2'] = path
            else:
                srr = f.split('.')[0]
                samples.setdefault(srr, {'1': None, '2': None})['1'] = path

        for srr, paths in samples.items():
            # Create a unique sample name combining GSE and SRR to avoid collisions in batches
            unique_sample_name = f"{gse_id}_{srr}"
            fq1 = paths['1']    
            fq2 = paths['2']
            if not fq1: continue 
            
            if fq2:
                rows.append([unique_sample_name, fq1, fq2, 'auto'])
            else:
                rows.append([unique_sample_name, fq1, '', 'auto'])
        
        return rows

    def run_pipeline_batch(self, samplesheet_path, batch_out_dir):
        """
        Runs nf-core/rnaseq on a combined samplesheet.
        """
        os.makedirs(batch_out_dir, exist_ok=True)
        project_root = os.getcwd() 
        config_path = os.path.join(project_root, ".nextflow.config")
        print(f"Running nf-core/rnaseq (Batch Mode) in {batch_out_dir}...")

        cmd = [
            "nextflow", "run", "nf-core/rnaseq",
            "-profile", self.profile,
            "-c", config_path,
            "-revision", "3.14.0",
            "-ansi-log", "false",
            "--slurm_account", "ewi-insy-prb",
            "--slurm_partition", "ewi-insy-prb,prb,ewi-insy,insy,general",
            "--input", samplesheet_path,
            "--outdir", batch_out_dir,
            "--pseudo_aligner", "salmon",
            "--skip_alignment",
            "--skip_biotype_qc",
            "--skip_stringtie",
            "--skip_bigwig",
            "--genome", "TAIR10",
            "--skip_fastqc",
            "--skip_multiqc",
            "--skip_dupradar",
            "--skip_qualimap",
            "--skip_rseqc",
            "-resume"
        ]
        
        try:
            subprocess.run(cmd, check=True)
            return True
        except subprocess.CalledProcessError as e:
            print(f"Nextflow Batch Error: {e}")
            return False

# --- BATCH HELPER FUNCTIONS ---

def split_merged_counts(batch_results_dir, study_map, output_root):
    """
    Reads the huge salmon.merged.gene_counts.tsv from the batch run
    and saves individual copies for each GSE study.
    
    study_map: dict {gse_id: [list of sample_ids (e.g. GSE123_SRR456)]}
    """
    merged_file = os.path.join(batch_results_dir, "star_salmon", "salmon.merged.gene_counts.tsv")
    if not os.path.exists(merged_file):
        merged_file = os.path.join(batch_results_dir, "salmon", "salmon.merged.gene_counts.tsv")
        
    if not os.path.exists(merged_file):
        print("Error: Merged count file not found in batch output.")
        return False

    print("Demultiplexing batch results...")
    df = pd.read_csv(merged_file, sep='\t')
    
    # Gene info columns
    meta_cols = ['gene_id', 'gene_name']
    
    for gse_id, samples in study_map.items():
        # Determine which columns belong to this study
        # The samplesheet named them f"{gse_id}_{srr}", so they should match df columns
        study_cols = [c for c in df.columns if c in samples]
        
        if not study_cols:
            print(f"Warning: No samples found in results for {gse_id}")
            continue
            
        # Create study output folder
        study_out = os.path.join(output_root, "processed_rnaseq", gse_id)
        os.makedirs(os.path.join(study_out, "star_salmon"), exist_ok=True)
        
        # Subset and Save
        study_df = df[meta_cols + study_cols]
        target_file = os.path.join(study_out, "star_salmon", "salmon.merged.gene_counts.tsv")
        study_df.to_csv(target_file, sep='\t', index=False)
        print(f"  Saved {gse_id} counts to {target_file}")
        
    return True

def chunk_list(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def download_experiments_RNA_seq_nf_core(gse_list:list[str], root_storage_dir:str, output_dir:str, tracker:FileTracker, download_raw:bool=True, scan:bool=True, run_and_delete:bool=True, batch_size:int=5):
    """
    Orchestrates the download and processing of RNA-Seq studies in BATCHES.
    """
    PATH_TO_INDEX = f"{root_storage_dir}genome_index/tair10"
    PATH_TO_GTF = f"{root_storage_dir}genome_index/Arabidopsis_thaliana.TAIR10.56.gtf"
    
    processor = RNASeq_processor(threads=4, genome_index=PATH_TO_INDEX, gtf_annotation=PATH_TO_GTF, profile='singularity')
    tracker_save_path = os.path.join(output_dir, "rnaseq_tracker_stats.json")
    valid_gse_ids = []

    # Filter list for things already processed
    todos = [g for g in gse_list if not tracker.is_processed(g) and not tracker.is_ignored(g)]
    
    print(f"Found {len(todos)} studies to process. Running in batches of {batch_size}.")

    for batch in chunk_list(todos, batch_size):
        
        batch_samplesheet_rows = []
        batch_study_map = {} # {gse_id: [sample_names]}
        batch_fastq_dirs = []
        
        print(f"\n=== Processing Batch: {batch} ===")
        
        # --- PHASE 1: DOWNLOAD & PREPARE ---
        for gse_id in batch:
            try:
                # 1. Setup Folders
                fastq_folder = os.path.join(output_dir, "fastq_storage", gse_id)
                cluster_temp = os.environ.get('TMPDIR', '/tmp')
                
                # 2. Metadata Check
                try:
                    gse = GEOparse.get_GEO(geo=gse_id, destdir=output_dir, silent=True)
                except:
                    print(f"Metadata error for {gse_id}")
                    tracker.mark_ignore(gse_id); continue

                if not check_metadata_for_sra_boolean(gse):
                    print(f"No SRA data for {gse_id}")
                    tracker.mark_ignore(gse_id); continue
                
                # 3. Download
                if download_raw:
                    if not tracker.is_downloaded(gse_id):
                        try:
                            processor.download_fastq(gse, fastq_folder, cluster_temp)
                            tracker.mark_downloaded(gse_id)
                        except Exception as e:
                            print(f"Download failed for {gse_id}: {e}")
                            tracker.mark_ignore(gse_id)
                            shutil.rmtree(fastq_folder, ignore_errors=True)
                            continue

                # 4. Generate Samplesheet Rows
                if os.path.exists(fastq_folder) and os.listdir(fastq_folder):
                    rows = processor.get_samplesheet_rows(gse_id, fastq_folder)
                    if rows:
                        batch_samplesheet_rows.extend(rows)
                        batch_study_map[gse_id] = [r[0] for r in rows] # Store sample names (col 0)
                        batch_fastq_dirs.append(fastq_folder)
                    else:
                        print(f"No valid FASTQ pairs found for {gse_id}")
                        tracker.mark_ignore(gse_id)
                else:
                    tracker.mark_ignore(gse_id)
            
            except Exception as e:
                print(f"Error preparing {gse_id}: {e}")
                tracker.mark_ignore(gse_id)

        # --- PHASE 2: EXECUTE BATCH ---
        if not batch_samplesheet_rows:
            print("Skipping batch execution (no valid samples).")
            continue

        # Write Combined Samplesheet
        batch_id = f"batch_{batch[0]}_{len(batch)}"
        batch_dir = os.path.join(output_dir, "batch_processing", batch_id)
        os.makedirs(batch_dir, exist_ok=True)
        
        samplesheet_path = os.path.join(batch_dir, "samplesheet.csv")
        with open(samplesheet_path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['sample', 'fastq_1', 'fastq_2', 'strandedness'])
            writer.writerows(batch_samplesheet_rows)

        # Run Nextflow
        success = processor.run_pipeline_batch(samplesheet_path, batch_dir)

        # --- PHASE 3: DISTRIBUTE RESULTS & CLEANUP ---
        if success:
            # Split the big results file back into study folders
            split_success = split_merged_counts(batch_dir, batch_study_map, output_dir)
            
            if split_success:
                for gse_id in batch_study_map.keys():
                    tracker.mark_processed(gse_id)
                    valid_gse_ids.append(gse_id)
                    
                    if run_and_delete:
                        fq_dir = os.path.join(output_dir, "fastq_storage", gse_id)
                        if os.path.exists(fq_dir):
                            print(f"Cleaning FASTQs for {gse_id}")
                            shutil.rmtree(fq_dir)
            
            tracker.save_to_json(tracker_save_path)
            
            # Optional: Cleanup batch folder to save space
            # shutil.rmtree(batch_dir) 
        else:
            print("Batch execution failed. Marking studies for review.")
            # Logic: You might want to retry individually or mark all as ignored.
            # Currently leaving them as downloaded but not processed.

    return valid_gse_ids