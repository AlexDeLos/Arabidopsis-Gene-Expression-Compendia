import contextlib
import csv
import glob
import json
import os
import re
import shutil
import ssl
import subprocess
import sys
import time
from collections import defaultdict
from collections.abc import Callable
from urllib.error import HTTPError

import GEOparse
import pandas as pd
from Bio import Entrez

Entrez.email = "A.DeLosSantosSubirats@tudelft.nl"


# Force Python to ignore SSL certificate verification globally
ssl._create_default_https_context = ssl._create_unverified_context
module_dir = "./"
sys.path.append(module_dir)

from src.data_importing.helpers.download_helper import check_metadata_for_sra_boolean  # noqa: E402
from src.data_importing.helpers.file_tracker import FileTracker  # noqa: E402
from src.data_importing.helpers.helpers import process_metadata  # noqa: E402

# ---------------------------------------------------------------------------
# FIX 1: Instrument/technology blocklist
# Studies sequenced on these platforms cannot be processed by this pipeline.
# SOLiD uses colorspace encoding; CAGE/RAMPAGE are not standard RNA-seq.
# ---------------------------------------------------------------------------
INCOMPATIBLE_INSTRUMENT_PATTERNS = [
    r"solid",  # AB SOLiD System, AB SOLiD 4 System, etc.
    r"ab solid",
    r"colorspace",
]

INCOMPATIBLE_LIBRARY_STRATEGIES = {
    "CAGE",
    "RAMPAGE",
    "ChIP-Seq",
    "ATAC-seq",
    "AMPLICON",
}


def is_study_compatible(gse) -> tuple[bool, str]:
    """
    Check whether a study uses a sequencing technology compatible with
    this pipeline (standard Illumina short-read RNA-seq).

    Returns (is_compatible: bool, reason: str).
    The reason is empty when compatible.
    """
    for gsm_id, gsm in list(gse.gsms.items()):
        m = gsm.metadata

        instrument = m.get("instrument_model", [""])[0].lower()
        for pattern in INCOMPATIBLE_INSTRUMENT_PATTERNS:
            if re.search(pattern, instrument):
                return False, f"Incompatible instrument '{instrument}' detected in {gsm_id}"

        strategy = m.get("library_strategy", [""])[0].strip()
        if strategy in INCOMPATIBLE_LIBRARY_STRATEGIES:
            return False, f"Incompatible library strategy '{strategy}' detected in {gsm_id}"

    return True, ""


class RNASeq_processor:
    def __init__(self, threads=4, genome_index=None, gtf_annotation=None, profile="docker"):
        self.threads = str(threads)
        self.profile = profile
        self.genome_index = genome_index
        self.gtf_annotation = gtf_annotation

        # Verify Tools
        required = ["fastq-dump", "trimmomatic", "hisat2", "samtools", "featureCounts"]
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

                if not record["IdList"]:
                    return []

                handle = Entrez.esummary(db="sra", id=",".join(record["IdList"]))
                summaries = Entrez.read(handle)
                handle.close()

                run_ids = []
                for summary in summaries:
                    run_ids.extend(re.findall(r'acc="([A-Z0-9]+)"', summary.get("Runs", "")))
                return list(set(run_ids))

            except HTTPError as e:
                # 2. Reactive Backoff: If we still hit the limit, wait and retry
                if e.code == 429:
                    wait_time = 2**attempt  # Waits 1s, 2s, 4s, 8s, 16s...
                    print(f"    [!] HTTP 429 (Too Many Requests) for {gsm_id}. Waiting {wait_time}s before retry {attempt + 1}/{max_retries}...")
                    time.sleep(wait_time)
                else:
                    print(f"    [!] HTTP Error {e.code} for {gsm_id}. Skipping.")
                    return []
            except Exception as e:
                print(f"    [!] Connection error for {gsm_id}: {e}. Retrying in 5s...")
                time.sleep(5)

        print(f"    [!] Failed to retrieve SRR for {gsm_id} after {max_retries} retries.")
        return []

    def download_fastq(self, gse, output_folder, temp_files, container,old):
        """Downloads using fastq-dump via parallel SLURM array jobs."""
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        if not os.path.exists(temp_files):
            os.makedirs(temp_files)

        # Directory for the SLURM out/err logs
        big_storage = os.environ.get("BIG_STORAGE", "/tudelft.net/staff-umbrella/GeneExpressionStorage")
        logs_folder = os.path.join(big_storage, "logs_slurm", "download_logs")
        if not os.path.exists(logs_folder):
            os.makedirs(logs_folder, exist_ok=True)

        print(f"Fetching SRR IDs for {len(gse.gsms)} samples...")
        sra_map = {gsm: self.get_srr_ids(gsm) for gsm in gse.gsms}

        # DEBUG 2: Print the raw map before filtering
        print(f"RAW SRA MAP: {sra_map}")

        sra_map = {k: v for k, v in sra_map.items() if v}

        # DEBUG 3: Print the map after filtering
        if not sra_map:
            print("CRITICAL WARNING: sra_map is empty. get_srr_ids failed to find runs.")
            return
        # Pre-validate any existing files and delete corrupt ones
        print(f"Pre-validating existing files in {output_folder}...")
        for fname in list(os.listdir(output_folder)):
            if not fname.endswith(".gz"):
                continue
            gz_path = os.path.join(output_folder, fname)
            try:
                subprocess.run(["gzip", "-t", "-q", gz_path], check=True, stderr=subprocess.PIPE)
                result = subprocess.run(  # noqa: PLW1510
                    ["bash", "-c", f'zcat "{gz_path}" | head -4 | wc -l'], capture_output=True, text=True, timeout=60
                )
                if result.returncode != 0 or int(result.stdout.strip()) < 4:
                    raise subprocess.CalledProcessError(1, "zcat")
            except (subprocess.CalledProcessError, subprocess.TimeoutExpired, ValueError):
                print(f"    [!] PRE-SCAN: {fname} is corrupt. Deleting.")
                os.remove(gz_path)

        # 1. Collect all SRRs that need downloading
        srrs_to_download = []
        for _gsm, srrs in sra_map.items():
            for srr in srrs:
                existing_gz = [f for f in os.listdir(output_folder) if f.startswith(srr) and f.endswith(".gz")]
                if not existing_gz:
                    srrs_to_download.append(srr)

        if not srrs_to_download:
            print(f"All SRRs for {gse} are already downloaded.")
            return

        # 2. Write missing SRRs to a text file for the SLURM array to read
        srr_list_path = os.path.join(output_folder, "srr_list.txt")
        with open(srr_list_path, "w") as f:
            for srr in srrs_to_download:
                f.write(f"{srr}\n")

        print(f"Submitting SLURM array job for {len(srrs_to_download)} SRRs...")
        if container:
            if old:
                pa = "slurm_jobs_2/download_srr_container.sbatch"
            else:
                pa = "slurm_jobs/download_srr_container.sbatch"
            print(f'using: {pa}')
            sbatch_script = os.path.abspath(os.path.join(module_dir, pa))
        else:
            sbatch_script = os.path.abspath(os.path.join(module_dir, "slurm_jobs/download_srr.sbatch"))

        if not os.path.exists(sbatch_script):
            print(f"CRITICAL ERROR: {sbatch_script} not found! Cannot execute download.")
            return

        # 3. Call sbatch and WAIT for all array jobs to finish
        cmd = [
            "sbatch",
            "--wait",  # Blocks the python script until the downloads finish
            f"--array=1-{len(srrs_to_download)}",
            f"--output={logs_folder}/fastq_dump_%A_%a.out",
            f"--error={logs_folder}/fastq_dump_%A_%a.err",
            sbatch_script,
            srr_list_path,
            output_folder,
        ]

        try:
            subprocess.run(cmd, check=True)
            print("All SLURM download array jobs completed.")
        except subprocess.CalledProcessError as e:
            print(f"Error executing sbatch job array: {e}")
            raise  # propagate so caller marks study as error (retryable), not downloaded

        # 4. Verify post-download and check data integrity
        successful_downloads = 0  # <-- Track successful SRRs

        for srr in srrs_to_download:
            existing_gz = [f for f in os.listdir(output_folder) if f.startswith(srr) and f.endswith(".gz")]

            valid_gz = []
            for gz_file in existing_gz:
                gz_path = os.path.join(output_folder, gz_file)
                try:
                    # Two-stage check: header test + partial decompression
                    subprocess.run(["gzip", "-t", "-q", gz_path], check=True, stderr=subprocess.PIPE)
                    # Also verify actual content is readable (catches truncated-but-valid-header files)
                    result = subprocess.run(  # noqa: PLW1510
                        ["bash", "-c", f'zcat "{gz_path}" | head -4 | wc -l'], capture_output=True, text=True, timeout=60
                    )
                    if result.returncode != 0 or int(result.stdout.strip()) < 4:
                        raise subprocess.CalledProcessError(1, "zcat")
                    valid_gz.append(gz_file)
                except (subprocess.CalledProcessError, subprocess.TimeoutExpired, ValueError):
                    print(f"    [!] CORRUPTION DETECTED: {gz_file} fails decompression check. Deleting.")
                    os.remove(gz_path)

            if not valid_gz:
                print(f"    [!] Failed to download valid files for {srr} after sbatch completion.")
            else:
                successful_downloads += 1  # <-- Increment if at least 1 valid file exists for this SRR

        # --- NEW BLOCK ---
        # If we attempted to download SRRs, but absolutely none survived the validation:
        if len(srrs_to_download) > 0 and successful_downloads == 0:
            error_msg = f"CRITICAL: All {len(srrs_to_download)} SRR downloads for {gse.name} failed or were completely corrupted."
            print(error_msg)
            raise RuntimeError(error_msg)
        # -----------------

        print(f"Done downloading for {gse.name}")
        print(f"Done downloading for {gse}")

    # ------------------------------------------------------------------
    # Private helpers for download_fastq
    # ------------------------------------------------------------------


    def _get_missing_or_corrupt_srrs(self, srrs: list, output_folder: str) -> list:
        """
        For each SRR, verify that at least one valid (non-corrupt) gz file
        exists.  Corrupt files are deleted so they can be re-downloaded.
        Returns the list of SRRs that still need downloading.
        """
        still_needed = []
        for srr in srrs:
            existing_gz = [f for f in os.listdir(output_folder) if f.startswith(srr) and f.endswith(".gz")]

            valid_gz = []
            for gz_file in existing_gz:
                gz_path = os.path.join(output_folder, gz_file)
                try:
                    subprocess.run(["gzip", "-t", "-q", gz_path], check=True, stderr=subprocess.PIPE)
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
        files = [f for f in os.listdir(fastq_folder) if f.endswith((".fastq", ".fq", ".gz"))]

        for f in files:
            path = os.path.abspath(os.path.join(fastq_folder, f))  # Use absolute paths for batching
            if "_1" in f:
                srr = f.split("_1")[0]
                samples.setdefault(srr, {"1": None, "2": None})["1"] = path
            elif "_2" in f:
                srr = f.split("_2")[0]
                samples.setdefault(srr, {"1": None, "2": None})["2"] = path
            else:
                srr = f.split(".")[0]
                samples.setdefault(srr, {"1": None, "2": None})["1"] = path

        for srr, paths in samples.items():
            # Create a unique sample name combining GSE and SRR to avoid collisions in batches
            unique_sample_name = f"{gse_id}_{srr}"
            fq1 = paths["1"]
            fq2 = paths["2"]
            if not fq1:
                continue

            if fq2:
                rows.append([unique_sample_name, fq1, fq2, "auto"])
            else:
                rows.append([unique_sample_name, fq1, "", "auto"])
        return rows

    def run_pipeline_batch(self, samplesheet_path, batch_out_dir, refs: dict,old:bool):
        os.makedirs(batch_out_dir, exist_ok=True)
        project_root = os.getcwd()
        if old:
            config_path = os.path.join(project_root, ".new_nextflow.config")
        else:
            config_path = os.path.join(project_root, ".nextflow.config")
        print(f"using config: {config_path}")
        # Unique log file inside batch_out_dir — safe across concurrent array jobs
        batch_name = os.path.basename(batch_out_dir)
        log_path = os.path.join(batch_out_dir, f"nextflow_{batch_name}.log")
        print(f"Running nf-core/rnaseq (Batch Mode) in {batch_out_dir}...")
        print(f"Nextflow log: {log_path}")

        all_bad_samples = set()
        current_samplesheet = samplesheet_path
        max_isolation_rounds = 10

        for round_num in range(max_isolation_rounds):
            success = self._run_nextflow(current_samplesheet, batch_out_dir, refs, config_path, log_path,old)

            if success:
                return True, all_bad_samples

            bad_samples = self._extract_failed_samples(log_path)
            if not bad_samples:
                print("  [!] No bad samples identified — unrecoverable failure.")
                return False, all_bad_samples

            new_bad = bad_samples - all_bad_samples
            if not new_bad:
                print("  [!] Same samples failing again — unrecoverable failure.")
                return False, all_bad_samples

            all_bad_samples.update(new_bad)
            print(f"  [!] Round {round_num + 1}: isolated {len(new_bad)} new failing sample(s): {new_bad}")
            print(f"  [!] Total bad samples so far: {all_bad_samples}")
            print("  [!] Removing them from the samplesheet and retrying...")

            retry_path = samplesheet_path.replace(".csv", f"_retry{round_num + 1}.csv")
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

    def _run_nextflow(self, samplesheet_path: str, batch_out_dir: str, refs: dict, config_path: str, log_path: str,old:bool) -> bool:
        """Execute the nextflow command and return True on success."""
        cmd = [
            "nextflow",
            "-log",
            log_path,  # <-- explicit log file, isolated per batch
            "run",
            "nf-core/rnaseq",
            "-resume",
            "-profile",
            self.profile,
            "-c",
            config_path,
            "-with-dag",
            f"{batch_out_dir}/flow_diagram.svg",
            "-revision",
            "3.14.0",  # New version might be worth use
            "-ansi-log",
            "false",
            "--slurm_account",
            "ewi-insy-prb" if old else 'testusers',
            "--slurm_partition",
            "general" if old else 'all',
            "--input",
            samplesheet_path,
            "--outdir",
            batch_out_dir,
            "--pseudo_aligner",
            "salmon",
            "--skip_alignment",
            "--fasta",
            refs["fasta"],
            "--gff",
            refs["gtf"],
            "--salmon_index",
            refs["salmon_index"],
            "--gtf_group_features_type",
            "mRNA",
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
                if "CONDA" in key:
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
        with open(log_path, errors="replace") as fh:
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
        with open(src_path, newline="") as fin, open(dst_path, "w", newline="") as fout:
            reader = csv.reader(fin)
            writer = csv.writer(fout)
            for i, row in enumerate(reader):
                if i == 0:  # header
                    writer.writerow(row)
                    continue
                if row and row[0] in bad_samples:
                    print(f"    Removing bad sample from samplesheet: {row[0]}")
                    removed += 1
                else:
                    writer.writerow(row)
        return removed


# --- BATCH HELPER FUNCTIONS ---


def split_merged_counts(batch_results_dir, study_map, output_root, batch_id=None):
    # nf-core rnaseq with --skip_alignment outputs to star_salmon/;
    # salmon-only runs output to salmon/. Check both, prefer star_salmon.
    merged_counts_file = os.path.join(batch_results_dir, "star_salmon", "salmon.merged.gene_counts.tsv")
    merged_tpm_file = os.path.join(batch_results_dir, "star_salmon", "salmon.merged.gene_tpm.tsv")
    salmon_subdir = "star_salmon"
    if not os.path.exists(merged_counts_file):
        merged_counts_file = os.path.join(batch_results_dir, "salmon", "salmon.merged.gene_counts.tsv")
        merged_tpm_file = os.path.join(batch_results_dir, "salmon", "salmon.merged.gene_tpm.tsv")
        salmon_subdir = "salmon"
    if not os.path.exists(merged_counts_file):
        print("Error: Merged count file not found in batch output.")
        return False

    # salmon metadata JSONs live in the same subdir as the merged counts
    salmon_base_dir = os.path.join(batch_results_dir, salmon_subdir)

    print(f"Demultiplexing batch results from {salmon_subdir}/...")
    df_counts = pd.read_csv(merged_counts_file, sep="\t")

    # Load TPM data if it exists
    df_tpm = None
    if os.path.exists(merged_tpm_file):
        df_tpm = pd.read_csv(merged_tpm_file, sep="\t")
    else:
        print("  [!] Warning: Merged TPM file not found in batch output.")

    # Only keep meta columns that actually exist in the dataframe
    meta_cols = [c for c in ["gene_id", "gene_name"] if c in df_counts.columns]

    saved = []
    for gse_id, samples in study_map.items():
        study_cols = [c for c in df_counts.columns if c in samples]

        if not study_cols:
            print(f"  Warning: No samples found in results for {gse_id}")
            continue

        study_out = os.path.join(output_root, "processed_rnaseq", gse_id)
        os.makedirs(os.path.join(study_out, "star_salmon"), exist_ok=True)

        meta_out_dir = os.path.join(study_out, "salmon_metadata")
        os.makedirs(meta_out_dir, exist_ok=True)

        # 1. Write a richer batch_info.txt that survives batch dir deletion
        batch_record_file = os.path.join(study_out, "batch_info.txt")
        with open(batch_record_file, "w") as f:
            f.write(f"Batch_ID:                {batch_id or 'unknown'}\n")
            f.write(f"Origin_Batch_Directory:  {batch_results_dir}\n")
            f.write(f"Salmon_Subdir:           {salmon_subdir}\n")
            f.write(f"Studies_In_Batch:        {', '.join(study_map.keys())}\n")
            f.write(f"Samples_This_Study:      {', '.join(study_cols)}\n")
            f.write(f"Timestamp_UTC:           {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}\n")

        # 2. Extract metrics AND copy the raw JSON files
        coverage_data = []
        for sample_col in study_cols:
            n_processed = "unspecified"
            n_mapped = "unspecified"
            pct_mapped = "unspecified"

            # Search in the correct subdir (star_salmon OR salmon)
            sample_folders = glob.glob(os.path.join(salmon_base_dir, f"*{sample_col}*"))

            if sample_folders:
                meta_path = os.path.join(sample_folders[0], "aux_info", "meta_info.json")
                if os.path.exists(meta_path):
                    # A) Copy the entire JSON file into the study directory
                    target_json = os.path.join(meta_out_dir, f"{sample_col}_meta_info.json")
                    shutil.copy2(meta_path, target_json)

                    # B) Extract key metrics for the summary CSV
                    try:
                        with open(meta_path) as fh:
                            meta_json = json.load(fh)
                            n_processed = meta_json.get("num_processed", "unspecified")
                            n_mapped = meta_json.get("num_mapped", "unspecified")
                            pct_mapped = meta_json.get("percent_mapped", "unspecified")
                    except Exception as e:
                        print(f"  Warning: Could not read meta_info for {sample_col}: {e}")
            else:
                print(f"  Warning: No salmon output folder found for sample {sample_col}")

            coverage_data.append(
                {
                    "Sample": sample_col,
                    "Total_Reads_Processed": n_processed,
                    "Reads_Mapped": n_mapped,
                    "Percent_Mapped": pct_mapped,
                    "Batch_ID": batch_id or "unknown",
                }
            )

        # Save coverage CSV
        coverage_df = pd.DataFrame(coverage_data)
        coverage_file = os.path.join(study_out, "sample_coverage.csv")
        coverage_df.to_csv(coverage_file, index=False)

        # 3. Save the Count Matrix
        study_df_counts = df_counts[meta_cols + study_cols]
        target_counts_file = os.path.join(study_out, "star_salmon", "salmon.merged.gene_counts.tsv")
        study_df_counts.to_csv(target_counts_file, sep="\t", index=False)

        # 4. Save the TPM Matrix
        # Re-derive study_cols from df_tpm in case column names differ slightly
        if df_tpm is not None:
            tpm_cols = [c for c in df_tpm.columns if c in samples]
            if tpm_cols:
                study_df_tpm = df_tpm[meta_cols + tpm_cols]
                target_tpm_file = os.path.join(study_out, "star_salmon", "salmon.merged.gene_tpm.tsv")
                study_df_tpm.to_csv(target_tpm_file, sep="\t", index=False)
                print(f"  Saved {gse_id} counts AND TPM to {study_out}/star_salmon/")
            else:
                print(f"  Saved {gse_id} counts to {target_counts_file} (TPM columns not found in df_tpm)")
        else:
            print(f"  Saved {gse_id} counts to {target_counts_file} (TPM file was missing)")

        print(f"  Saved coverage + Salmon metadata to {study_out}/")
        saved.append(gse_id)

    return saved


def chunk_list(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i : i + n]


def get_ecotype_from_gse(gse) -> str:
    """
    Extracts the reference genome/ecotype from GEO metadata.
    Returns a string key like 'col-0', 'ler', 'ws', or 'unknown'.
    """
    # Check study-level characteristics
    sources_to_check = []

    # 1. Check overall study metadata
    for field in ["overall_design", "summary", "title"]:
        sources_to_check.extend(gse.metadata.get(field, []))

    # 2. Check per-sample characteristics (more reliable)
    for gsm in list(gse.gsms.values())[:3]:  # Check first 3 samples
        for key, val_list in gsm.metadata.items():
            if "characteristics" in key or "source" in key:
                sources_to_check.extend(val_list)

    text = " ".join(sources_to_check).lower()

    # Order matters: check specific ecotypes before generic 'col'
    ECOTYPE_PATTERNS = {
        "col-0": [r"col-0", r"columbia-0", r"columbia\b", r"col\b"],
        "ler": [r"landsberg", r"\bler\b", r"ler-0"],
        "ws": [r"\bws\b", r"wassilewskija"],
        "c24": [r"\bc24\b"],
        "cvi": [r"\bcvi\b", r"cape verde"],
    }

    for ecotype, patterns in ECOTYPE_PATTERNS.items():
        if any(re.search(p, text) for p in patterns):
            return ecotype

    return "col-0"  # Default to Col-0 — safe assumption for most studies


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

        study_title = gse.metadata.get("title", [""])[0]
        study_platform = gse.metadata.get("platform_id", [""])[0]
        study_summary = gse.metadata.get("summary", [""])[0]

        rows = []
        for gsm_id, gsm in gse.gsms.items():
            m = gsm.metadata

            # --- Fixed fields (always present in GEO) ---
            row = {
                "study_id": gse_id,
                "study_title": study_title,
                "study_summary": study_summary,
                "geo_accession": gsm_id,
                "title": m.get("title", [""])[0],
                "source_name": m.get("source_name_ch1", [""])[0],
                "organism": m.get("organism_ch1", [""])[0],
                "platform": study_platform,
                "library_strategy": m.get("library_strategy", [""])[0],
                "library_layout": m.get("library_layout", [""])[0],
                "instrument_model": m.get("instrument_model", [""])[0],
            }

            # --- Parse characteristics_ch1 as individual columns ---
            for char in m.get("characteristics_ch1", []):
                if ":" in char:
                    key, _, value = char.partition(":")
                    col_name = key.strip().lower().replace(" ", "_")
                    row[col_name] = value.strip()
                else:
                    row[f"characteristic_{len(row)}"] = char.strip()

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


def download_experiments_RNA_seq_nf_core(
    gse_list: list[str],
    root_storage_dir: str,
    output_dir: str,
    tracker: FileTracker,
    download_raw: bool = True,
    metadata_only: bool = True,
    run_and_delete: bool = True,
    batch_size: int = 5,
    debug: bool = False,
    container: bool = False,
    old:bool = False
):
    """
    Orchestrates the download and processing of RNA-Seq studies in BATCHES.
    """
    _debug = debug
    PATH_TO_INDEX = f"{root_storage_dir}genome_index/tair10"
    PATH_TO_GTF = f"{root_storage_dir}genome_index/Arabidopsis_thaliana.TAIR10.56.gtf"
    REFERENCE_MAP = {
        "col-0": {
            "fasta": f"{root_storage_dir}files_for_rna_seq/tair12/tair12.fasta.gz",
            "gtf": f"{root_storage_dir}files_for_rna_seq/tair12/tair12_annotation.gff.gz",
            "salmon_index": f"{root_storage_dir}files_for_rna_seq/tair12/salmon_index",
        },
        "ler": {
            "fasta": f"{root_storage_dir}files_for_rna_seq/tair12/tair12.fasta.gz",
            "gtf": f"{root_storage_dir}files_for_rna_seq/tair12/tair12_annotation.gff.gz",
            "salmon_index": f"{root_storage_dir}files_for_rna_seq/tair12/salmon_index",
        },
        # Default fallback
        "unknown": {
            "fasta": f"{root_storage_dir}files_for_rna_seq/tair12/tair12.fasta.gz",
            "gtf": f"{root_storage_dir}files_for_rna_seq/tair12/tair12_annotation.gff.gz",
            "salmon_index": f"{root_storage_dir}files_for_rna_seq/tair12/salmon_index",
        },
    }
    processor = RNASeq_processor(threads=4, genome_index=PATH_TO_INDEX, gtf_annotation=PATH_TO_GTF, profile="singularity,slurm")
    tracker_save_path = os.path.join(output_dir, "rnaseq_tracker_stats.json")
    valid_gse_ids = []

    # Filter list for things already processed
    todos = [g for g in gse_list if (not tracker.is_processed(g) or metadata_only) and not tracker.is_ignored(g) and not tracker.is_error(g)]

    skipped_processed = [g for g in gse_list if tracker.is_processed(g) and not metadata_only]
    skipped_ignored = [g for g in gse_list if tracker.is_ignored(g)]
    skipped_error = [g for g in gse_list if tracker.is_error(g)]

    print(f"To process ({len(todos)}): {todos}")
    print(f"Skipped — already processed ({len(skipped_processed)}): {skipped_processed}")
    print(f"Skipped — ignored ({len(skipped_ignored)}): {skipped_ignored}")
    print(f"Skipped — error ({len(skipped_error)}): {skipped_error}")

    ecotype_groups: dict[str, list[str]] = defaultdict(list)

    print("Detecting ecotypes for all studies...")
    for gse_id in todos:
        ecotype = "col-0"  # safe default always set first
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
                    ecotype = "col-0"
                tracker.mark_ecotype(gse_id, ecotype)
                print(f"  {gse_id} -> {ecotype} (detected)")
        except Exception as e:
            print(f"  WARNING: unexpected error for {gse_id} ({e}), defaulting to col-0")
            ecotype = "col-0"
        finally:
            ecotype_groups[ecotype].append(gse_id)  # ALWAYS appended no matter what

    for ecotype, ids in ecotype_groups.items():
        print(f"  Ecotype '{ecotype}': {len(ids)} studies -> {ids}")

    for ecotype, gse_ids_for_ecotype in ecotype_groups.items():
        refs = REFERENCE_MAP.get(ecotype, REFERENCE_MAP["col-0"])
        print(f"\n{'=' * 60}")
        print(f"  ECOTYPE GROUP: {ecotype} ({len(gse_ids_for_ecotype)} studies)")
        print(f"  Reference FASTA:  {refs['fasta']}")
        print(f"  Salmon Index:     {refs['salmon_index']}")
        print(f"{'=' * 60}")

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
                    cluster_temp = os.environ.get("TMPDIR", "/tmp")
                    try:
                        gse = GEOparse.get_GEO(geo=gse_id, destdir=output_dir, silent=True)
                    except Exception:
                        print(f"Metadata error for {gse_id}")
                        tracker.mark_ignore(gse_id)
                        continue

                    if not check_metadata_for_sra_boolean(gse):
                        print(f"No SRA data for {gse_id}")
                        tracker.mark_ignore(gse_id)
                        continue
                    if len(gse.metadata.get("type", [])) > 1:
                        print(f"Skipping {gse.name}: Contains multiple experiment types.")
                        tracker.mark_ignore(gse_id)
                        continue
                    if len(gse.gsms) < 5:  # type: ignore  # noqa: PGH003
                        print("We are ignoring this study because it is < 5")
                        tracker.mark_ignore(gse_id)
                        continue

                    # --- technology compatibility check ---
                    compatible, reason = is_study_compatible(gse)
                    if not compatible:
                        print(f"  [!] Skipping {gse_id}: {reason}")
                        tracker.mark_ignore(gse_id)
                        continue

                    print("  - Processing metadata for all samples in study...")
                    try:
                        for _, gsm in gse.gsms.items():  # type: ignore  # noqa: PGH003
                            process_metadata(gse_id, gse, gsm, save_path=os.path.join(output_dir, "metadata", gse_id))
                    except Exception as e:
                        print(f"    > Metadata processing failed: {e}")

                    if metadata_only:
                        continue
                    # save_rnaseq_sample_metadata(gse_id, gse, output_dir)
                    # 3. Download (FIX 2: retry on corruption baked into download_fastq)
                    if download_raw and not tracker.is_downloaded(gse_id):
                        try:
                            processor.download_fastq(gse, fastq_folder, cluster_temp, container,old)
                            tracker.mark_downloaded(gse_id)
                            print(f"Download completed for {gse_id}")
                        except Exception as e:
                            print(f"Download failed for {gse_id}: {e}")
                            tracker.mark_error(gse_id)  # transient failure — allow retry
                            shutil.rmtree(fastq_folder, ignore_errors=True)
                            continue

                    if os.path.exists(fastq_folder) and os.listdir(fastq_folder):
                        print(f"Generating sample sheet rows for: {gse_id}")
                        rows = processor.get_samplesheet_rows(gse_id, fastq_folder)
                        if rows:
                            batch_samplesheet_rows.extend(rows)
                            batch_study_map[gse_id] = [r[0] for r in rows]
                            batch_fastq_dirs.append(fastq_folder)
                            print(f"DONE generating sample sheet for: {gse_id}")
                        else:
                            print(f"No valid FASTQ pairs found for {gse_id}")
                            tracker.mark_ignore(gse_id)
                    else:
                        print(f"fastq_folder={fastq_folder} | exists={os.path.exists(fastq_folder)}")
                        tracker.mark_ignore(gse_id)
                except Exception as e:
                    print(f"Error preparing {gse_id}: {e}")
                    tracker.mark_error(gse_id)  # transient failure — allow retry
            if metadata_only:
                continue
            # --- PHASE 2: EXECUTE BATCH ---
            if not batch_samplesheet_rows:
                print("Skipping batch execution (no valid samples).")
                continue

            # Build a collision-proof batch ID.
            # Format: batch_{ecotype}_{slurm_array_task_id}_{utc_timestamp}_{batch_index}
            # - SLURM_ARRAY_TASK_ID differentiates parallel array jobs
            # - UTC timestamp (to-the-second) makes it unique even if the same
            #   array task re-runs (e.g. after a failure)
            # - batch_index (position within this task's ecotype loop) handles
            #   the case where one task processes multiple ecotype groups
            slurm_task = os.environ.get("SLURM_ARRAY_TASK_ID", "local")
            utc_ts = time.strftime("%Y%m%dT%H%M%SZ", time.gmtime())
            batch_index = sum(len(v) for v in ecotype_groups.values() if v is not gse_ids_for_ecotype)
            batch_id = f"batch_{ecotype}_{slurm_task}_{utc_ts}_{batch_index:02d}"
            batch_dir = os.path.join(output_dir, "batch_processing", batch_id)
            os.makedirs(batch_dir, exist_ok=True)

            samplesheet_path = os.path.join(batch_dir, "samplesheet.csv")
            with open(samplesheet_path, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["sample", "fastq_1", "fastq_2", "strandedness"])
                writer.writerows(batch_samplesheet_rows)
                f.flush()
                os.fsync(f.fileno())

            print("Waiting 10 seconds for umbrella drive to sync samplesheet...")
            time.sleep(10)

            # FIX 3: run_pipeline_batch now returns (success, bad_samples)
            success, bad_samples = processor.run_pipeline_batch(samplesheet_path, batch_dir, refs,old)

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

                split_success = split_merged_counts(batch_dir, effective_study_map, output_dir, batch_id=batch_id)
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
                            metadata_dir = os.path.join(output_dir, "metadata", gse_id)
                            if os.path.exists(soft_file):
                                print(f"Cleaning metadata files for {gse_id}")
                                shutil.rmtree(metadata_dir)

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
                                keep = (
                                    (name == "deseq2.plots.pdf" and "deseq2_qc" in root)
                                    or (name == "meta_info.json" and "aux_info" in root)
                                    or (name == "multiqc_report.html")
                                    or ("software_versions" in name)
                                    or (name in ["salmon.merged.transcript_counts.tsv", "salmon.merged.transcript_tpm.tsv"])
                                )
                                if not keep:
                                    with contextlib.suppress(OSError):
                                        os.remove(filepath)
                            for name in dirs:
                                with contextlib.suppress(OSError):
                                    os.rmdir(os.path.join(root, name))
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
                    for gse_id in batch_study_map:
                        if bad_sample.startswith(gse_id):
                            guilty_gses.add(gse_id)
                            break

                # If we couldn't identify any guilty GSE (e.g. no bad_samples
                # were extracted from the log), mark everything as error to
                # avoid an infinite retry loop
                if not guilty_gses:
                    print("  [!] Could not identify guilty GSEs — marking all as error.")
                    guilty_gses = set(batch_study_map.keys())

                for gse_id in batch_study_map:
                    if gse_id in guilty_gses:
                        print(f"  Marking {gse_id} as error (caused pipeline failure).")
                        tracker.mark_error(gse_id)
                    else:
                        # Leave at STATUS_DOWNLOADED — will be retried next run
                        print(f"  Leaving {gse_id} as downloaded (collateral, will retry).")

    return valid_gse_ids


def combine_files_rnaseq(
    folder: str,
    new_file_name: str,
    new_file_location: str,
    combination_method: str | Callable = "sum",  # Changed default to 'sum' for TPMs
) -> tuple[pd.DataFrame, dict[str, str]]:
    """
    Combines individual GSE RNA-seq TSV files (from star_salmon) into a single large expression matrix
    and generates a mapping of Samples to Study IDs.

    Args:
        folder (str): Root directory containing GSE subfolders (e.g., processed_rnaseq).
        new_file_name (str): Output filename (e.g. 'RMA_RNAseq_Combined.csv').
        new_file_location (str): Output directory.
        combine_genes (bool): If True, merges 'Gene.1', 'Gene.2' into 'Gene'.
        combination_method (str or callable): Method to merge overlapping values (should be 'sum' for TPM).

    Returns:
        Tuple[pd.DataFrame, Dict]: (The combined expression dataframe, Dictionary {SampleID: StudyID})
    """

    # 1. Setup Output Directory
    if not os.path.exists(new_file_location):
        try:
            os.makedirs(new_file_location)
        except OSError as e:
            print(f"Error creating directory {new_file_location}: {e}")
            raise e

    dataframes = []
    sample_to_study_map = {}

    print(f"Scanning '{folder}' for processed RNA-seq files...")

    # 2. Iterate over the folder structure
    for d in os.listdir(folder):
        study_folder = os.path.join(folder, d)

        if os.path.isdir(study_folder):
            study_id = d

            expected_tsv_name = "salmon.merged.gene_tpm.tsv"
            file_path = os.path.join(study_folder, "star_salmon", expected_tsv_name)

            if os.path.exists(file_path):
                print(f"  - Found counts for Study: {study_id}")
                try:
                    df = pd.read_csv(file_path, sep="\t", index_col=0)

                    if "gene_name" in df.columns:
                        df = df.drop(columns=["gene_name"])

                    df.index = df.index.astype(str).str.upper()

                    # --- ARABIDOPSIS-SPECIFIC NCBI CLEANUP ---
                    # 1. Extract base AGI code (removes NCBI prefixes and isoform suffixes)
                    df.index = df.index.str.extract(r"(AT[1-5CM]G\d+)", expand=False)

                    # 2. Drop non-coding RNAs / unmapped features (which become NaN)
                    df = df[df.index.notna()]

                    # 3. Sum isoforms into gene-level TPMs (or use user's combination method)
                    df = df.groupby(df.index).agg(combination_method)
                    # -----------------------------------------

                    for sample_id in df.columns:
                        sample_to_study_map[sample_id] = study_id

                    dataframes.append(df)

                except Exception as e:
                    print(f"    ! Error reading {expected_tsv_name} in {study_id}: {e}")

    if not dataframes:
        msg = "No valid TSV files found to combine."
        print(msg)
        raise ValueError(msg)

    # 3. Combine Dataframes
    print(f"\nMerging {len(dataframes)} datasets... (This may take a moment)")

    # Because we already resolved duplicates in step 2, this merge is perfectly clean
    combined_df = pd.concat(dataframes, axis=1, join="outer", sort=True)

    # Fallback to handle any cross-dataset duplicates if they somehow exist
    if combined_df.index.duplicated().any():
        print(f"  - resolving duplicate indices in final merge (using {combination_method})...")
        combined_df = combined_df.groupby(combined_df.index).agg(combination_method)

    # 5. Save Data Matrix
    output_path = os.path.join(new_file_location, new_file_name)
    print(f"Saving combined matrix to: {output_path}")

    try:
        combined_df.to_csv(output_path)
        print("SUCCESS: Data matrix saved.")
        print(f"Final Dimensions: {combined_df.shape[0]} Genes x {combined_df.shape[1]} Samples")
    except Exception as e:
        print(f"Error saving data file: {e}")

    # 6. Save Sample Map
    base_name = os.path.splitext(new_file_name)[0]
    map_filename = f"{base_name}_sample_map.csv"
    map_output_path = os.path.join(new_file_location, map_filename)

    print(f"Saving sample map to: {map_output_path}")
    try:
        map_df = pd.DataFrame.from_dict(sample_to_study_map, orient="index", columns=["StudyID"])
        map_df.index.name = "SampleID"
        map_df.to_csv(map_output_path)
        print("SUCCESS: Sample map saved.")
    except Exception as e:
        print(f"Error saving map file: {e}")

    return combined_df, sample_to_study_map
