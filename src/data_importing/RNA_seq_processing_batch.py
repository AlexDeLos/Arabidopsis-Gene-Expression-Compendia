import os
import subprocess
import shutil
import csv
import re
import pandas as pd
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
from src.data_importing.helpers.helpers import process_metadata
from src.constants import *
from src.data_importing.helpers.file_tracker import FileTracker


# ---------------------------------------------------------------------------
# FIX 1: Instrument/technology blocklist
# Studies sequenced on these platforms cannot be processed by this pipeline.
# SOLiD uses colorspace encoding; CAGE/RAMPAGE are not standard RNA-seq.
# ---------------------------------------------------------------------------
INCOMPATIBLE_INSTRUMENT_PATTERNS = [
    r'solid',       # AB SOLiD System, AB SOLiD 4 System, etc.
    r'ab solid',
    r'colorspace',
]

INCOMPATIBLE_LIBRARY_STRATEGIES = {
    'CAGE', 'RAMPAGE', 'ChIP-Seq', 'ATAC-seq', 'AMPLICON',
}

def is_study_compatible(gse) -> tuple[bool, str]:
    """
    Check whether a study uses a sequencing technology compatible with
    this pipeline (standard Illumina short-read RNA-seq).

    Returns (is_compatible: bool, reason: str).
    The reason is empty when compatible.
    """
    for gsm_id, gsm in list(gse.gsms.items())[:5]:   # check first 5 samples
        m = gsm.metadata

        instrument = m.get('instrument_model', [''])[0].lower()
        for pattern in INCOMPATIBLE_INSTRUMENT_PATTERNS:
            if re.search(pattern, instrument):
                return False, f"Incompatible instrument '{instrument}' detected in {gsm_id}"

        strategy = m.get('library_strategy', [''])[0].strip()
        if strategy in INCOMPATIBLE_LIBRARY_STRATEGIES:
            return False, f"Incompatible library strategy '{strategy}' detected in {gsm_id}"

    return True, ""


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
            raise  # propagate so caller marks study as error (retryable), not downloaded

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

    # ------------------------------------------------------------------
    # Private helpers for download_fastq
    # ------------------------------------------------------------------

    def _run_sbatch_download(self, srrs: list, output_folder: str, logs_folder: str):
        """Write an SRR list file and block until the SLURM array finishes."""
        srr_list_path = os.path.join(output_folder, "srr_list.txt")
        with open(srr_list_path, 'w') as f:
            for srr in srrs:
                f.write(f"{srr}\n")

        print(f"Submitting SLURM array job for {len(srrs)} SRRs...")
        sbatch_script = os.path.abspath(os.path.join(module_dir, "slurm_jobs/download_srr.sbatch"))

        if not os.path.exists(sbatch_script):
            print(f"CRITICAL ERROR: {sbatch_script} not found! Cannot execute download.")
            return
        
        cmd = [
            "sbatch",
            "--wait",
            f"--array=1-{len(srrs)}",
            f"--output={logs_folder}/fastq_dump_%A_%a.out",
            f"--error={logs_folder}/fastq_dump_%A_%a.err",
            sbatch_script,
            srr_list_path,
            output_folder
        ]

        # --- NEW HANDOFF LOGIC ---
        # Get the current array ID to prevent parallel jobs from overwriting each other
        task_id = os.environ.get('SLURM_ARRAY_TASK_ID', 'dev')
        submit_script_path = f"submit_slurm_array_{task_id}.sh"

        try:
            with open(submit_script_path, 'w') as f:
                f.write("#!/bin/bash\n")
                # Join the command list into a single string
                f.write(" ".join(cmd) + "\n")
            print(f"Apptainer finished. Generated host submission script: {submit_script_path}")
        except Exception as e:
            print(f"Error generating submission script: {e}")

    def _get_missing_or_corrupt_srrs(self, srrs: list, output_folder: str) -> list:
        """
        For each SRR, verify that at least one valid (non-corrupt) gz file
        exists.  Corrupt files are deleted so they can be re-downloaded.
        Returns the list of SRRs that still need downloading.
        """
        still_needed = []
        for srr in srrs:
            existing_gz = [
                f for f in os.listdir(output_folder)
                if f.startswith(srr) and f.endswith('.gz')
            ]

            valid_gz = []
            for gz_file in existing_gz:
                gz_path = os.path.join(output_folder, gz_file)
                try:
                    subprocess.run(
                        ['gzip', '-t', '-q', gz_path],
                        check=True, stderr=subprocess.PIPE
                    )
                    valid_gz.append(gz_file)
                except subprocess.CalledProcessError:
                    print(f"    [!] CORRUPTION DETECTED: {gz_file} is incomplete/corrupted. Deleting.")
                    os.remove(gz_path)

            if not valid_gz:
                still_needed.append(srr)

        return still_needed

            
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

    def run_pipeline_batch(self, samplesheet_path, batch_out_dir, refs: dict):
        os.makedirs(batch_out_dir, exist_ok=True)
        project_root = os.getcwd()
        config_path = os.path.join(project_root, ".nextflow.config")

        # Unique log file inside batch_out_dir — safe across concurrent array jobs
        batch_name = os.path.basename(batch_out_dir)
        log_path = os.path.join(batch_out_dir, f"nextflow_{batch_name}.log")
        print(f"Running nf-core/rnaseq (Batch Mode) in {batch_out_dir}...")
        print(f"Nextflow log: {log_path}")

        all_bad_samples = set()
        current_samplesheet = samplesheet_path
        max_isolation_rounds = 10

        for round_num in range(max_isolation_rounds):
            success = self._run_nextflow(current_samplesheet, batch_out_dir, refs, config_path, log_path)

            if success:
                return True, all_bad_samples

            bad_samples = self._extract_failed_samples(log_path)   # <-- pass path directly

            if not bad_samples:
                print("  [!] No bad samples identified — unrecoverable failure.")
                return False, all_bad_samples

            new_bad = bad_samples - all_bad_samples
            if not new_bad:
                print("  [!] Same samples failing again — unrecoverable failure.")
                return False, all_bad_samples

            all_bad_samples.update(new_bad)
            print(f"  [!] Round {round_num+1}: isolated {len(new_bad)} new failing sample(s): {new_bad}")
            print(f"  [!] Total bad samples so far: {all_bad_samples}")
            print(f"  [!] Removing them from the samplesheet and retrying...")

            retry_path = samplesheet_path.replace('.csv', f'_retry{round_num+1}.csv')
            removed = self._write_samplesheet_without(samplesheet_path, retry_path, all_bad_samples)

            if removed == 0:
                print("  [!] Could not match failing samples to samplesheet rows — cannot retry.")
                return False, all_bad_samples

            current_samplesheet = retry_path

        print(f"  [!] Reached max isolation rounds ({max_isolation_rounds}).")
        return False, all_bad_samples

    # ------------------------------------------------------------------
    # Private helpers for run_pipeline_batch
    # ------------------------------------------------------------------

    def _run_nextflow(self, samplesheet_path: str, batch_out_dir: str, refs: dict, config_path: str, log_path: str) -> bool:
        """Execute the nextflow command and return True on success."""
        cmd = [
            "nextflow",
            "-log", log_path,          # <-- explicit log file, isolated per batch
            "run", "nf-core/rnaseq",
            "-profile", self.profile,
            "-c", config_path,
            "-with-dag", f'{batch_out_dir}/flow_diagram.svg',
            "-revision", "3.14.0",
            "-ansi-log", "false",
            "--slurm_account", "ewi-insy-prb",
            "--slurm_partition", "ewi-insy-prb,prb,ewi-insy,insy",
            "--input", samplesheet_path,
            "--outdir", batch_out_dir,
            "--pseudo_aligner", "salmon",
            "--skip_alignment",
            "--fasta",        refs['fasta'],
            "--gtf",          refs['gtf'],
            "--salmon_index", refs['salmon_index'],
            "--gtf_group_features_type", "mRNA",
            "--skip_biotype_qc",
            "--skip_stringtie",
            "--skip_bigwig",
            "--skip_fastqc",
            "--skip_multiqc",
            "--skip_dupradar",
            "--skip_qualimap",
            "--skip_rseqc",
        ]

        try:
            clean_env = os.environ.copy()
            for key in list(clean_env.keys()):
                if 'CONDA' in key:
                    del clean_env[key]
            subprocess.run(cmd, check=True, env=clean_env, cwd=batch_out_dir)
            return True
        except subprocess.CalledProcessError as e:
            print(f"Nextflow Batch Error: {e}")
            return False

    def _extract_failed_samples(self, log_path: str) -> set:
        """
        Parse a Nextflow log file to extract sample names from ERROR lines like:
        ERROR ~ Error executing process > 'TRIMGALORE (GSE280648_SRR31170500)'
        Returns a set of sample name strings.
        """
        bad = set()
        if not os.path.exists(log_path):
            print(f"  [!] Nextflow log not found at {log_path}; cannot identify failing samples.")
            return bad

        pattern = re.compile(r"Error executing process\s*>\s*'[^']*\(([^)]+)\)'")
        with open(log_path, 'r', errors='replace') as fh:
            for line in fh:
                m = pattern.search(line)
                if m:
                    bad.add(m.group(1).strip())
        return bad

    def _write_samplesheet_without(self, src_path: str, dst_path: str, bad_samples: set) -> int:
        """
        Copy the samplesheet from src_path to dst_path, omitting rows
        whose sample name (column 0) is in bad_samples.
        Returns the number of rows removed.
        """
        removed = 0
        with open(src_path, newline='') as fin, open(dst_path, 'w', newline='') as fout:
            reader = csv.reader(fin)
            writer = csv.writer(fout)
            for i, row in enumerate(reader):
                if i == 0:   # header
                    writer.writerow(row)
                    continue
                if row and row[0] in bad_samples:
                    print(f"    Removing bad sample from samplesheet: {row[0]}")
                    removed += 1
                else:
                    writer.writerow(row)
        return removed


# --- BATCH HELPER FUNCTIONS ---

def split_merged_counts(batch_results_dir, study_map, output_root):
    merged_file = os.path.join(batch_results_dir, "star_salmon", "salmon.merged.gene_counts.tsv")
    if not os.path.exists(merged_file):
        merged_file = os.path.join(batch_results_dir, "salmon", "salmon.merged.gene_counts.tsv")
    if not os.path.exists(merged_file):
        print("Error: Merged count file not found in batch output.")
        return False

    print("Demultiplexing batch results...")
    df = pd.read_csv(merged_file, sep='\t')

    # Only keep meta columns that actually exist in the dataframe
    # pseudo-alignment produces 'gene_id' only; full alignment adds 'gene_name'
    meta_cols = [c for c in ['gene_id', 'gene_name'] if c in df.columns]

    saved = []
    for gse_id, samples in study_map.items():
        study_cols = [c for c in df.columns if c in samples]

        if not study_cols:
            print(f"  Warning: No samples found in results for {gse_id}")
            continue

        study_out = os.path.join(output_root, "processed_rnaseq", gse_id)
        os.makedirs(os.path.join(study_out, "star_salmon"), exist_ok=True)

        study_df = df[meta_cols + study_cols]
        target_file = os.path.join(study_out, "star_salmon", "salmon.merged.gene_counts.tsv")
        study_df.to_csv(target_file, sep='\t', index=False)
        print(f"  Saved {gse_id} counts to {target_file}")
        saved.append(gse_id)

    return saved  # return list of saved GSEs instead of bare True

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
            for char in m.get('characteristics_ch1', []):
                if ':' in char:
                    key, _, value = char.partition(':')
                    col_name = key.strip().lower().replace(' ', '_')
                    row[col_name] = value.strip()
                else:
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


def download_experiments_RNA_seq_nf_core(gse_list:list[str], root_storage_dir:str, output_dir:str, tracker:FileTracker, download_raw:bool=True, metadata_only:bool=True, run_and_delete:bool=True, batch_size:int=5,debug:bool=False):
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
    todos = [g for g in gse_list if (not tracker.is_processed(g) or metadata_only) and not tracker.is_ignored(g) and not tracker.is_error(g)]
    print(f'we are going to process these studies {todos}')
    from collections import defaultdict
    ecotype_groups: dict[str, list[str]] = defaultdict(list)
    
    print("Detecting ecotypes for all studies...")
    for gse_id in todos:
        ecotype = 'col-0'  # safe default always set first
        try:
            cached = tracker.get_ecotype(gse_id)
            if cached is not None:
                ecotype = cached
                print(f"  {gse_id} -> {ecotype} (cached)")
            else:
                try:
                    gse = GEOparse.get_GEO(geo=gse_id, destdir=output_dir, silent=True)
                    ecotype = get_ecotype_from_gse(gse)
                except Exception as e:
                    print(f"  WARNING: ecotype detection failed for {gse_id} ({e}), defaulting to col-0")
                    ecotype = 'col-0'
                tracker.mark_ecotype(gse_id, ecotype)
                print(f"  {gse_id} -> {ecotype} (detected)")
        except Exception as e:
            print(f"  WARNING: unexpected error for {gse_id} ({e}), defaulting to col-0")
            ecotype = 'col-0'
        finally:
            ecotype_groups[ecotype].append(gse_id)  # ALWAYS appended no matter what

    for ecotype, ids in ecotype_groups.items():
        print(f"  Ecotype '{ecotype}': {len(ids)} studies -> {ids}")

    for ecotype, gse_ids_for_ecotype in ecotype_groups.items():
        refs = REFERENCE_MAP.get(ecotype, REFERENCE_MAP['col-0'])
        print(f"\n{'='*60}")
        print(f"  ECOTYPE GROUP: {ecotype} ({len(gse_ids_for_ecotype)} studies)")
        print(f"  Reference FASTA:  {refs['fasta']}")
        print(f"  Salmon Index:     {refs['salmon_index']}")
        print(f"{'='*60}")

        for batch in chunk_list(gse_ids_for_ecotype, batch_size):
            
            batch_samplesheet_rows = []
            batch_study_map = {}
            batch_fastq_dirs = []
            
            print(f"\n=== Processing Batch [{ecotype}]: {batch} ===")
            
            # --- PHASE 1: DOWNLOAD & PREPARE ---
            for gse_id in batch:
                print(f"\n=== Processing study: {gse_id} ===")
                if tracker.is_ignored(gse_id):
                    print(f"Ignore made it this far, why? {gse_id}")
                try:
                    fastq_folder = os.path.join(output_dir, "fastq_storage", gse_id)
                    cluster_temp = os.environ.get('TMPDIR', '/tmp')
                    
                    try:
                        gse = GEOparse.get_GEO(geo=gse_id, destdir=output_dir, silent=True)
                    except:
                        print(f"Metadata error for {gse_id}")
                        tracker.mark_ignore(gse_id); continue

                    if not check_metadata_for_sra_boolean(gse):
                        print(f"No SRA data for {gse_id}")
                        tracker.mark_ignore(gse_id); continue
                    if len(gse.gsms) < 5: # type: ignore
                        tracker.mark_ignore(gse_id); continue

                    # --- FIX 1: technology compatibility check ---
                    compatible, reason = is_study_compatible(gse)
                    if not compatible:
                        print(f"  [!] Skipping {gse_id}: {reason}")
                        tracker.mark_ignore(gse_id)
                        continue
                    
                    print("  - Processing metadata for all samples in study...")
                    try:
                        for _, gsm in gse.gsms.items(): # type: ignore
                            process_metadata(gse_id, gse, gsm, save_path=os.path.join(output_dir, "metadata", gse_id))
                    except Exception as e:
                        print(f"    > Metadata processing failed: {e}")
                    
                    if metadata_only:
                        continue #TODO: test if this works, should only download the metadata, nothing more
                    # save_rnaseq_sample_metadata(gse_id, gse, output_dir)
                    # 3. Download (FIX 2: retry on corruption baked into download_fastq)
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
            if metadata_only:
                continue
            # --- PHASE 2: EXECUTE BATCH ---
            if not batch_samplesheet_rows:
                print("Skipping batch execution (no valid samples).")
                continue

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

            # FIX 3: run_pipeline_batch now returns (success, bad_samples)
            success, bad_samples = processor.run_pipeline_batch(samplesheet_path, batch_dir, refs)

            # --- PHASE 3: DISTRIBUTE RESULTS & CLEANUP ---
            if success:
                # Exclude bad samples from demultiplexing
                effective_study_map = {}
                for gse_id, samples in batch_study_map.items():
                    clean_samples = [s for s in samples if s not in bad_samples]
                    if clean_samples:
                        effective_study_map[gse_id] = clean_samples
                    else:
                        print(f"  [!] All samples for {gse_id} were removed — marking as error.")
                        tracker.mark_error(gse_id)

                split_success = split_merged_counts(batch_dir, effective_study_map, output_dir)
                # split_success is a list of saved GSE IDs, or False if merged file not found

                if split_success is not False:
                    for gse_id in split_success:
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

                    # Studies with no matching columns stay at STATUS_DOWNLOADED for retry
                    not_saved = set(effective_study_map.keys()) - set(split_success)
                    for gse_id in not_saved:
                        print(f"  [!] No columns found for {gse_id} in merged counts — leaving as downloaded.")

                    tracker.save_to_json(tracker_save_path)

                    if run_and_delete and split_success:
                        print(f"Trimming batch directory {batch_dir} to save space (Keeping only QC logs)...")
                        for root, dirs, files in os.walk(batch_dir, topdown=False):
                            for name in files:
                                filepath = os.path.join(root, name)
                                keep = (name == "deseq2.plots.pdf" and "deseq2_qc" in root) or \
                                       (name == "meta_info.json" and "aux_info" in root)
                                if not keep:
                                    try:
                                        os.remove(filepath)
                                    except OSError:
                                        pass
                            for name in dirs:
                                try:
                                    os.rmdir(os.path.join(root, name))
                                except OSError:
                                    pass
                else:
                    # Merged count file not found at all — leave everything as downloaded for retry
                    print("Error: Merged count file not found — leaving all studies as downloaded.")
                    tracker.save_to_json(tracker_save_path)
            else:
                # Pipeline failed even after bad-sample isolation.
                # bad_samples contains the GSE_SRR names that caused failures.
                # Only mark the GSEs that owned those bad samples as error.
                # All other GSEs in the batch stay at STATUS_DOWNLOADED so
                # they will be retried in the next run.

                # Identify which GSEs are responsible for the bad samples
                guilty_gses = set()
                for bad_sample in bad_samples:
                    # Sample names are formatted as GSE{id}_SRR{id}
                    for gse_id in batch_study_map.keys():
                        if bad_sample.startswith(gse_id):
                            guilty_gses.add(gse_id)
                            break

                # If we couldn't identify any guilty GSE (e.g. no bad_samples
                # were extracted from the log), mark everything as error to
                # avoid an infinite retry loop
                if not guilty_gses:
                    print("  [!] Could not identify guilty GSEs — marking all as error.")
                    guilty_gses = set(batch_study_map.keys())

                for gse_id in batch_study_map.keys():
                    if gse_id in guilty_gses:
                        print(f"  Marking {gse_id} as error (caused pipeline failure).")
                        tracker.mark_error(gse_id)
                    else:
                        # Leave at STATUS_DOWNLOADED — will be retried next run
                        print(f"  Leaving {gse_id} as downloaded (collateral, will retry).")

    return valid_gse_ids
