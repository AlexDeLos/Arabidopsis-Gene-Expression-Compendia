import os
import subprocess
import shutil
import csv
import re
import pandas as pd
from Bio import Entrez
from tqdm import tqdm
import GEOparse
from urllib.error import HTTPError
import sys
import time
from Bio import Entrez
Entrez.email = "A.DeLosSantosSubirats@tudelft.nl"
import ssl
# Force Python to ignore SSL certificate verification globally
ssl._create_default_https_context = ssl._create_unverified_context
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
        required = ['fastq-dump','trimmomatic', 'hisat2', 'samtools', 'featureCounts']
        for tool in required:
            if not shutil.which(tool):
                print(f"WARNING: {tool} not found in PATH. Pipeline may fail.")

    def get_srr_ids(self, gsm_id, max_retries=5):
        for attempt in range(max_retries):
            try:
                # 1. Proactive Rate Limiting (ensures < 3 requests per second)
                time.sleep(0.4) 
                
                handle = Entrez.esearch(db="sra", term=gsm_id)
                record = Entrez.read(handle)
                handle.close()
                
                if not record['IdList']:
                    return []
        
                handle = Entrez.esummary(db="sra", id=",".join(record['IdList']))
                summaries = Entrez.read(handle)
                handle.close()
                
                run_ids = []
                import re
                for summary in summaries:
                    run_ids.extend(re.findall(r'acc="([A-Z0-9]+)"', summary.get('Runs', '')))
                return list(set(run_ids))
            
            except HTTPError as e:
                # 2. Reactive Backoff: If we still hit the limit, wait and retry
                if e.code == 429:
                    wait_time = 2 ** attempt  # Waits 1s, 2s, 4s, 8s, 16s...
                    print(f"    [!] HTTP 429 (Too Many Requests) for {gsm_id}. Waiting {wait_time}s before retry {attempt+1}/{max_retries}...")
                    time.sleep(wait_time)
                else:
                    print(f"    [!] HTTP Error {e.code} for {gsm_id}. Skipping.")
                    return []
            except Exception as e:
                print(f"    [!] Connection error for {gsm_id}: {e}. Retrying in 5s...")
                time.sleep(5)
                
        print(f"    [!] Failed to retrieve SRR for {gsm_id} after {max_retries} retries.")
        return []

                
    def download_fastq(self, gse, output_folder, temp_files):
        """Downloads using fastq-dump via parallel SLURM array jobs."""
        if not os.path.exists(output_folder): os.makedirs(output_folder)
        if not os.path.exists(temp_files): os.makedirs(temp_files)
        
        # Directory for the SLURM out/err logs
        logs_folder = os.path.join('./logs_slurm/', "download_logs")
        if not os.path.exists(logs_folder): os.makedirs(logs_folder)

        print(f"Fetching SRR IDs for {len(gse.gsms)} samples...")
        sra_map = {gsm: self.get_srr_ids(gsm) for gsm in gse.gsms.keys()}
        
        # DEBUG 2: Print the raw map before filtering
        print(f"RAW SRA MAP: {sra_map}")

        sra_map = {k:v for k,v in sra_map.items() if v}
        
        # DEBUG 3: Print the map after filtering
        if not sra_map:
            print("CRITICAL WARNING: sra_map is empty. get_srr_ids failed to find runs.")
            return

        # 1. Collect all SRRs that need downloading
        srrs_to_download = []
        for gsm, srrs in sra_map.items():
            for srr in srrs:
                existing_gz = [f for f in os.listdir(output_folder) if f.startswith(srr) and f.endswith('.gz')]
                if not existing_gz:
                    srrs_to_download.append(srr)

        if not srrs_to_download:
            print(f"All SRRs for {gse} are already downloaded.")
            return

        # 2. Write missing SRRs to a text file for the SLURM array to read
        srr_list_path = os.path.join(output_folder, "srr_list.txt")
        with open(srr_list_path, 'w') as f:
            for srr in srrs_to_download:
                f.write(f"{srr}\n")

        print(f"Submitting SLURM array job for {len(srrs_to_download)} SRRs...")
        sbatch_script = os.path.abspath(os.path.join(module_dir, "slurm_jobs/download_srr.sbatch"))
        
        if not os.path.exists(sbatch_script):
            print(f"CRITICAL ERROR: {sbatch_script} not found! Cannot execute download.")
            return

        # 3. Call sbatch and WAIT for all array jobs to finish
        cmd = [
            "sbatch", 
            "--wait", # Blocks the python script until the downloads finish
            f"--array=1-{len(srrs_to_download)}",
            f"--output={logs_folder}/fastq_dump_%A_%a.out",
            f"--error={logs_folder}/fastq_dump_%A_%a.err",
            sbatch_script,
            srr_list_path,
            output_folder
        ]
        
        try:
            subprocess.run(cmd, check=True)
            print("All SLURM download array jobs completed.")
        except subprocess.CalledProcessError as e:
            print(f"Error executing sbatch job array: {e}")

        # 4. Verify post-download and check data integrity
        for srr in srrs_to_download:
            existing_gz = [f for f in os.listdir(output_folder) if f.startswith(srr) and f.endswith('.gz')]
            
            valid_gz = []
            for gz_file in existing_gz:
                gz_path = os.path.join(output_folder, gz_file)
                try:
                    # Run a quiet integrity check on the gzip archive
                    subprocess.run(['gzip', '-t', '-q', gz_path], check=True, stderr=subprocess.PIPE)
                    valid_gz.append(gz_file)
                except subprocess.CalledProcessError:
                    print(f"    [!] CORRUPTION DETECTED: {gz_file} is incomplete/corrupted. Deleting to prevent pipeline crash.")
                    os.remove(gz_path)
            
            if not valid_gz:
                print(f"    [!] Failed to download valid files for {srr} after sbatch completion.")
                
        print(f'Done downloading for {gse.name}')
                
        print(f'Done downloading for {gse}')
            
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

    def run_pipeline_batch(self, samplesheet_path, batch_out_dir,refs:dict):
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
            # "-preview",
            "-with-dag", f'{batch_out_dir}/flow_diagram.svg',
            "-revision", "3.14.0",
            "-ansi-log", "false",
            "--slurm_account", "ewi-insy-prb",
            "--slurm_partition", "ewi-insy-prb,prb,ewi-insy,insy",
            "--input", samplesheet_path,
            "--outdir", batch_out_dir,
            
            # --- ALIGNMENT STRATEGY (MAX SPEED) ---
            "--pseudo_aligner", "salmon",
            "--skip_alignment",            
            
            # --- REFERENCE GENOME & ANNOTATIONS ---
            "--fasta",        refs['fasta'],
            "--gtf",          refs['gtf'],
            "--salmon_index", refs['salmon_index'],
            "--gtf_group_features_type", "mRNA",
            # "--fasta", "/tudelft.net/staff-umbrella/GeneExpressionStorage/files_for_rna_seq/col-0.fasta",
            # "--gtf", "/tudelft.net/staff-umbrella/GeneExpressionStorage/files_for_rna_seq/col-0_liftoff_polished_sorted_tbtools_clean.gtf",
            
            # Note: --salmon_index and --transcript_fasta have been intentionally removed.
            # Nextflow will automatically build a new index from the fasta and gtf above.
            
            # --- SKIP UNNECESSARY HEAVY QC & STEPS ---
            "--skip_biotype_qc",
            "--skip_stringtie",
            "--skip_bigwig",
            "--skip_fastqc",
            "--skip_multiqc",
            "--skip_dupradar",
            "--skip_qualimap",
            "--skip_rseqc",
            
            # "-resume"  <-- Keep this commented out for the very first run to ensure a clean start
        ]
        
        try:
            clean_env = os.environ.copy()

            # 2. Delete any variable that mentions CONDA so it doesn't leak
            for key in list(clean_env.keys()):
                if 'CONDA' in key:
                    del clean_env[key]
            subprocess.run(cmd, check=True, env=clean_env, cwd=batch_out_dir)
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

def get_ecotype_from_gse(gse) -> str:
    """
    Extracts the reference genome/ecotype from GEO metadata.
    Returns a string key like 'col-0', 'ler', 'ws', or 'unknown'.
    """
    # Check study-level characteristics
    sources_to_check = []
    
    # 1. Check overall study metadata
    for field in ['overall_design', 'summary', 'title']:
        sources_to_check.extend(gse.metadata.get(field, []))
    
    # 2. Check per-sample characteristics (more reliable)
    for gsm in list(gse.gsms.values())[:3]:  # Check first 3 samples
        for key, val_list in gsm.metadata.items():
            if 'characteristics' in key or 'source' in key:
                sources_to_check.extend(val_list)
    
    text = ' '.join(sources_to_check).lower()
    
    # Order matters: check specific ecotypes before generic 'col'
    ECOTYPE_PATTERNS = {
        'col-0':  [r'col-0', r'columbia-0', r'columbia\b', r'col\b'],
        'ler':    [r'landsberg', r'\bler\b', r'ler-0'],
        'ws':     [r'\bws\b', r'wassilewskija'],
        'c24':    [r'\bc24\b'],
        'cvi':    [r'\bcvi\b', r'cape verde'],
    }
    
    for ecotype, patterns in ECOTYPE_PATTERNS.items():
        if any(re.search(p, text) for p in patterns):
            return ecotype
    
    return 'col-0'  # Default to Col-0 — safe assumption for most studies

def save_rnaseq_sample_metadata(gse_id: str, gse, output_dir: str) -> str | None:
    """
    Extracts and saves sample-level metadata for an RNA-seq study to a CSV file.
    Mirrors the structure produced by process_metadata() in the microarray pipeline.

    Saved to: {output_dir}/metadata_rnaseq/{gse_id}/{gse_id}_sample_metadata.csv

    Columns extracted per GSM (when available):
        geo_accession, title, source_name, organism, platform,
        library_strategy, library_layout, instrument_model,
        + all 'characteristics_ch1' fields parsed as individual key:value columns

    Returns the path to the saved CSV, or None on failure.
    """
    try:
        save_dir = os.path.join(output_dir, "metadata_rnaseq", gse_id)
        os.makedirs(save_dir, exist_ok=True)
        save_path = os.path.join(save_dir, f"{gse_id}_sample_metadata.csv")

        # Skip if already saved (idempotent)
        if os.path.exists(save_path):
            print(f"  Metadata already saved for {gse_id}, skipping.")
            return save_path

        study_title    = gse.metadata.get('title', [''])[0]
        study_platform = gse.metadata.get('platform_id', [''])[0]
        study_summary  = gse.metadata.get('summary', [''])[0]

        rows = []
        for gsm_id, gsm in gse.gsms.items():
            m = gsm.metadata

            # --- Fixed fields (always present in GEO) ---
            row = {
                'study_id':        gse_id,
                'study_title':     study_title,
                'study_summary':   study_summary,
                'geo_accession':   gsm_id,
                'title':           m.get('title',            [''])[0],
                'source_name':     m.get('source_name_ch1',  [''])[0],
                'organism':        m.get('organism_ch1',     [''])[0],
                'platform':        study_platform,
                'library_strategy':   m.get('library_strategy',   [''])[0],
                'library_layout':     m.get('library_layout',     [''])[0],
                'instrument_model':   m.get('instrument_model',   [''])[0],
            }

            # --- Parse characteristics_ch1 as individual columns ---
            # GEOparse stores these as a list of "key: value" strings, e.g.:
            # ["genotype: wild type", "treatment: drought stress", "tissue: leaf"]
            for char in m.get('characteristics_ch1', []):
                if ':' in char:
                    key, _, value = char.partition(':')
                    col_name = key.strip().lower().replace(' ', '_')
                    row[col_name] = value.strip()
                else:
                    # Unparseable characteristic — store as-is under a generic key
                    row[f'characteristic_{len(row)}'] = char.strip()

            rows.append(row)

        if not rows:
            print(f"  No samples found for {gse_id}, skipping metadata save.")
            return None

        df = pd.DataFrame(rows)
        df.to_csv(save_path, index=False)
        print(f"  Sample metadata saved: {save_path} ({len(df)} samples)")
        return save_path

    except Exception as e:
        print(f"  WARNING: Could not save metadata for {gse_id}: {e}")
        return None


def download_experiments_RNA_seq_nf_core(gse_list:list[str], root_storage_dir:str, output_dir:str, tracker:FileTracker, download_raw:bool=True, scan:bool=True, run_and_delete:bool=True, batch_size:int=5,debug:bool=False):
    """
    Orchestrates the download and processing of RNA-Seq studies in BATCHES.
    """
    PATH_TO_INDEX = f"{root_storage_dir}genome_index/tair10"
    PATH_TO_GTF = f"{root_storage_dir}genome_index/Arabidopsis_thaliana.TAIR10.56.gtf"
    REFERENCE_MAP = {
        'col-0': {
            'fasta':         f"{root_storage_dir}files_for_rna_seq/col-0/col-0.fasta",
            'gtf':           f"{root_storage_dir}files_for_rna_seq/col-0/col-0_nfcore.gtf",
            'salmon_index':  f"{root_storage_dir}files_for_rna_seq/col-0/salmon_index",
        },
        #TODO: for now it is pointing to col-0 ecotype
        'ler': {
            'fasta':        f"{root_storage_dir}files_for_rna_seq/col-0/col-0.fasta",
            'gtf':          f"{root_storage_dir}files_for_rna_seq/col-0/col-0_nfcore.gtf",
            'salmon_index': f"{root_storage_dir}files_for_rna_seq/col-0/salmon_index",
        },
        # Default fallback
        'unknown': {
            'fasta':         f"{root_storage_dir}files_for_rna_seq/col-0/col-0.fasta",
            'gtf':           f"{root_storage_dir}files_for_rna_seq/col-0/col-0_nfcore.gtf",
            'salmon_index':  f"{root_storage_dir}files_for_rna_seq/col-0/salmon_index",
        },
    }
    processor = RNASeq_processor(threads=4, genome_index=PATH_TO_INDEX, gtf_annotation=PATH_TO_GTF, profile='singularity,slurm')
    tracker_save_path = os.path.join(output_dir, "rnaseq_tracker_stats.json")
    valid_gse_ids = []

    # Filter list for things already processed
    todos = [g for g in gse_list if not tracker.is_processed(g) and not tracker.is_ignored(g)]
        # --- PHASE 0: DETECT AND GROUP BY ECOTYPE ---
    # Group studies by their ecotype so each batch shares one reference genome
    from collections import defaultdict
    ecotype_groups: dict[str, list[str]] = defaultdict(list)
    
    print("Detecting ecotypes for all studies...")
    for gse_id in todos:
        cached = tracker.get_ecotype(gse_id)
        if cached is not None:
            ecotype = cached
            print(f"  {gse_id} -> {ecotype} (cached)")
        else:
            try:
                gse = GEOparse.get_GEO(geo=gse_id, destdir=output_dir, silent=True)
                ecotype = get_ecotype_from_gse(gse)
            except Exception:
                ecotype = 'col-0'  # Safe default
            tracker.mark_ecotype(gse_id, ecotype)
            print(f"  {gse_id} -> {ecotype} (detected)")

        ecotype_groups[ecotype].append(gse_id)

    # Log ecotype summary before processing
    for ecotype, ids in ecotype_groups.items():
        print(f"  Ecotype '{ecotype}': {len(ids)} studies -> {ids}")

    # Iterate over each ecotype group, then chunk within that group.
    # This guarantees every batch is homogeneous and uses the correct reference.
    for ecotype, gse_ids_for_ecotype in ecotype_groups.items():
        refs = REFERENCE_MAP.get(ecotype, REFERENCE_MAP['col-0'])
        print(f"\n{'='*60}")
        print(f"  ECOTYPE GROUP: {ecotype} ({len(gse_ids_for_ecotype)} studies)")
        print(f"  Reference FASTA:  {refs['fasta']}")
        print(f"  Salmon Index:     {refs['salmon_index']}")
        print(f"{'='*60}")

        for batch in chunk_list(gse_ids_for_ecotype, batch_size):
            
            batch_samplesheet_rows = []
            batch_study_map = {} # {gse_id: [sample_names]}
            batch_fastq_dirs = []
            
            print(f"\n=== Processing Batch [{ecotype}]: {batch} ===")
            
            # --- PHASE 1: DOWNLOAD & PREPARE ---
            for gse_id in batch:
                print(f"\n=== Processing study: {gse_id} ===")
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
                    if len(gse.gsms) < 5:
                        tracker.mark_ignore(gse_id); continue

                    # 2b. Save sample-level metadata (mirrors microarray process_metadata)
                    save_rnaseq_sample_metadata(gse_id, gse, output_dir)

                    if debug:
                        # --- DEBUG: KEEP ONLY 1 SAMPLE ---
                        if len(gse.gsms) > 1:
                            first_sample_id = list(gse.gsms.keys())[0]
                            gse.gsms = {first_sample_id: gse.gsms[first_sample_id]}
                            print(f"DEBUG MODE: Reduced {gse_id} to single sample: {first_sample_id}")
                        # ---------------------------------

                    # 3. Download
                    if download_raw:
                        if not tracker.is_downloaded(gse_id):
                            try:
                                processor.download_fastq(gse, fastq_folder, cluster_temp)
                                tracker.mark_downloaded(gse_id)
                                print(f"Download completed for {gse_id}")
                            except Exception as e:
                                print(f"Download failed for {gse_id}: {e}")
                                tracker.mark_ignore(gse_id)
                                shutil.rmtree(fastq_folder, ignore_errors=True)
                                continue

                    # 4. Generate Samplesheet Rows
                    if os.path.exists(fastq_folder) and os.listdir(fastq_folder):
                        print(f'Generating sample sheet rows for: {gse_id}')
                        rows = processor.get_samplesheet_rows(gse_id, fastq_folder)
                        if rows:
                            batch_samplesheet_rows.extend(rows)
                            batch_study_map[gse_id] = [r[0] for r in rows]
                            batch_fastq_dirs.append(fastq_folder)
                            print(f'DONE generating sample sheet for: {gse_id}')
                        else:
                            print(f"No valid FASTQ pairs found for {gse_id}")
                            tracker.mark_ignore(gse_id)
                    else:
                        print(f'fastq_folder={fastq_folder} | exists={os.path.exists(fastq_folder)}')
                        tracker.mark_ignore(gse_id)
                
                except Exception as e:
                    print(f"Error preparing {gse_id}: {e}")
                    tracker.mark_ignore(gse_id)

            # --- PHASE 2: EXECUTE BATCH ---
            if not batch_samplesheet_rows:
                print("Skipping batch execution (no valid samples).")
                continue

            # Write Combined Samplesheet
            # Batch ID encodes ecotype so dirs are never shared across ecotype runs
            batch_id = f"batch_{ecotype}_{batch[0]}_{len(batch)}"
            batch_dir = os.path.join(output_dir, "batch_processing", batch_id)
            os.makedirs(batch_dir, exist_ok=True)
            
            samplesheet_path = os.path.join(batch_dir, "samplesheet.csv")
            with open(samplesheet_path, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(['sample', 'fastq_1', 'fastq_2', 'strandedness'])
                writer.writerows(batch_samplesheet_rows)
                f.flush()
                os.fsync(f.fileno())

            print("Waiting 10 seconds for umbrella drive to sync samplesheet...")
            time.sleep(10)

            # BUG FIX: pass `refs` (single ecotype dict), not the entire REFERENCE_MAP
            success = processor.run_pipeline_batch(samplesheet_path, batch_dir, refs)

            # --- PHASE 3: DISTRIBUTE RESULTS & CLEANUP ---
            if success:
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
                                
                            soft_file = os.path.join(output_dir, f"{gse_id}_family.soft.gz")
                            if os.path.exists(soft_file):
                                print(f"Cleaning SOFT file for {gse_id}")
                                os.remove(soft_file)
                
                tracker.save_to_json(tracker_save_path)
                
                if split_success and run_and_delete:
                    print(f"Trimming batch directory {batch_dir} to save space (Keeping only QC logs)...")
                    for root, dirs, files in os.walk(batch_dir, topdown=False):
                        for name in files:
                            filepath = os.path.join(root, name)
                            keep = False
                            if name == "deseq2.plots.pdf" and "deseq2_qc" in root:
                                keep = True
                            elif name == "meta_info.json" and "aux_info" in root:
                                keep = True
                            if not keep:
                                try:
                                    os.remove(filepath)
                                except OSError:
                                    pass
                        for name in dirs:
                            dirpath = os.path.join(root, name)
                            try:
                                os.rmdir(dirpath)
                            except OSError:
                                pass
                else:
                    print('Count matrices for studies were not generated.')
                    print('Nothing was deleted and study tracker states were not changed.')
            else:
                print("Batch execution failed. Marking studies for review.")
                for gse_id in batch_study_map.keys():
                    tracker.mark_error(gse_id)

    return valid_gse_ids