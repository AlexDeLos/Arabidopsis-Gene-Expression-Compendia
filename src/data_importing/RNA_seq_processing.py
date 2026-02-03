import os
import subprocess
import shutil
import csv  # Added for samplesheet generation
from Bio import Entrez
from tqdm import tqdm
import GEOparse
from GEOparse.GEOTypes import GSE,GSM, GPL, GDS
import sys
module_dir = './'
sys.path.append(module_dir)
from src.data_importing.helpers.download_helper import check_metadata_for_sra_boolean
from src.constants import *
from src.data_importing.helpers.file_tracker import FileTracker




class RNASeq_processor:
    def __init__(self, threads=4, genome_index=None, gtf_annotation=None, profile='docker'):
        self.threads = str(threads)
        # self.genome_index = genome_index # Path to HISAT2 index prefix
        # self.gtf_annotation = gtf_annotation # Path to .gtf file
        self.profile = profile
        
        # Verify Tools
        required = ['fasterq-dump', 'trimmomatic', 'hisat2', 'samtools', 'featureCounts']
        for tool in required:
            if not shutil.which(tool):
                print(f"WARNING: {tool} not found in PATH. Pipeline will fail.")

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

    def download_fastq(self, gse, output_folder, temp_files):
        """
        Downloads using fasterq-dump (fast) and immediately compresses using pigz (fast).
        """
        if not os.path.exists(output_folder): os.makedirs(output_folder)
        if not os.path.exists(temp_files): os.makedirs(temp_files)

        sra_map = {gsm: self.get_srr_ids(gsm) for gsm in gse.gsms.keys()}
        sra_map = {k:v for k,v in sra_map.items() if v}
        
        fasterq_path = os.path.expanduser("~/fasterq-dump")

        # Check if 'pigz' is installed (It is standard on most clusters)
        # If not, fallback to 'gzip'
        zip_cmd = "pigz" if shutil.which("pigz") else "gzip"
        print(f"Using {zip_cmd} for compression.")

        for gsm, srrs in tqdm(sra_map.items(), desc="Downloading SRRs", leave=False):
            for srr in srrs:
                # 1. Check if the ZIPPED file already exists
                # We look for files starting with SRR... and ending with .gz
                existing_gz = [f for f in os.listdir(output_folder) if f.startswith(srr) and f.endswith('.gz')]
                if existing_gz: 
                    continue

                # 2. Run fasterq-dump (Produces raw .fastq)
                if CLUSTER_RUN:
                    cmd = [fasterq_path, "--split-files", "--outdir", output_folder, "--temp", temp_files , "--threads", self.threads, srr]
                else:
                    cmd = ["fasterq-dump", "--split-files", "--outdir", output_folder, "--temp", temp_files , "--threads", self.threads, srr]
                
                try:
                    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
                    
                    # 3. Compression Step
                    # Find the raw FASTQ files for this SRR
                    raw_fastqs = [f for f in os.listdir(output_folder) if f.startswith(srr) and f.endswith('.fastq')]
                    
                    if not raw_fastqs:
                        print(f"Warning: Download ran but no .fastq files found for {srr}")
                        continue

                    for fq in raw_fastqs:
                        full_path = os.path.join(output_folder, fq)
                        try:
                            # -f forces overwrite, -p 8 uses 8 cores (if using pigz)
                            # We adjust args based on which tool we found
                            if zip_cmd == "pigz":
                                subprocess.run(["pigz", "-f", "-p", self.threads, full_path], check=True)
                            else:
                                subprocess.run(["gzip", "-f", full_path], check=True)
                                
                        except subprocess.CalledProcessError:
                            print(f"Error compressing {fq}")
                            # Optional: cleanup failed zip?
                            
                except subprocess.CalledProcessError as e:
                    print(f"\nCRITICAL ERROR downloading {srr}:")
                    print(f"Command tried: {' '.join(cmd)}")
                    raise e

    def generate_samplesheet(self, gse_id, fastq_folder, output_csv):
        """
        Scans the FASTQ folder and creates a valid nf-core samplesheet.csv.
        Columns: sample,fastq_1,fastq_2,strandedness
        """
        samples = {}
        # Identify samples and pair files
        files = [f for f in os.listdir(fastq_folder) if f.endswith('.fastq') or f.endswith('.fq') or f.endswith('.gz')]
        
        for f in files:
            path = os.path.join(fastq_folder, f)
            # Basic parsing logic for SRR/ERR files
            if '_1' in f:
                srr = f.split('_1')[0]
                samples.setdefault(srr, {'1': None, '2': None})['1'] = path
            elif '_2' in f:
                srr = f.split('_2')[0]
                samples.setdefault(srr, {'1': None, '2': None})['2'] = path
            else:
                # Single end assumption if no _1/_2
                srr = f.split('.')[0]
                samples.setdefault(srr, {'1': None, '2': None})['1'] = path

        rows = []
        for srr, paths in samples.items():
            unique_sample_name = f"{gse_id}_{srr}"
            fq1 = paths['1']    
            fq2 = paths['2']
            if not fq1: continue # Skip if no R1
            
            # Strandedness 'auto' lets Salmon/Nextflow decide
            if fq2:
                rows.append([unique_sample_name, fq1, fq2, 'auto'])
            else:
                rows.append([unique_sample_name, fq1, '', 'auto'])
        
        if not rows:
            return False

        with open(output_csv, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['sample', 'fastq_1', 'fastq_2', 'strandedness'])
            writer.writerows(rows)
            
        return True

    def run_pipeline_on_study(self, gse_id, fastq_folder, output_base_dir):
        """
        Runs nf-core/rnaseq using Kallisto-style pseudo-alignment (via Salmon).
        """
        study_out_dir = output_base_dir #os.path.join(output_base_dir, gse_id)
        os.makedirs(study_out_dir, exist_ok=True)

        print(f"Checking for uncompressed FASTQs in {fastq_folder}...")
        files_to_compress = [f for f in os.listdir(fastq_folder) if f.endswith('.fastq') or f.endswith('.fq')]
        
        if files_to_compress:
            print(f"Compressing {len(files_to_compress)} files (Required for nf-core)...")
            # Use parallel gzip (pigz) if available, else gzip
            # We use standard gzip here for compatibility
            for f in tqdm(files_to_compress, desc="Gzipping"):
                full_path = os.path.join(fastq_folder, f)
                try:
                    subprocess.run(["gzip", "-f", full_path], check=True)
                except subprocess.CalledProcessError:
                    print(f"Failed to compress {f}")
                    return False
                
        
        # 1. Generate Samplesheet
        samplesheet_path = os.path.join(study_out_dir, "samplesheet.csv")
        has_samples = self.generate_samplesheet(gse_id, fastq_folder, samplesheet_path)
        
        if not has_samples:
            print(f"No samples found for {gse_id}")
            return False

        print(f"Running nf-core/rnaseq (Salmon mode) on {gse_id}...")

        # 2. Construct Nextflow Command
        # Note: --pseudo_aligner salmon activates the fast "Kallisto-like" mode
        # We skip alignment (STAR/HISAT2) and other QC to focus on counts.
        cmd = [
            "nextflow", "run", "nf-core/rnaseq",
            "-profile", self.profile,
            "-revision", "3.14.0",
            "--input", samplesheet_path,
            "--outdir", study_out_dir,
            "--pseudo_aligner", "kallisto",
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
            # Run Nextflow
            subprocess.run(cmd, check=True)
            
            # Check if expected output exists
            # Standard Salmon output in nf-core/rnaseq
            expected_counts = os.path.join(study_out_dir, "star_salmon", "salmon.merged.gene_counts.tsv")
            if not os.path.exists(expected_counts):
                # Fallback location depending on version
                expected_counts = os.path.join(study_out_dir, "salmon", "salmon.merged.gene_counts.tsv")
            
            if os.path.exists(expected_counts):
                print(f"Success: Count matrix found at {expected_counts}")
                return True
            else:
                print(f"Pipeline finished but count matrix not found in standard paths.")
                return False
                
        except subprocess.CalledProcessError as e:
            print(f"Nextflow Pipeline Error for {gse_id}: {e}")
            return False
    


#outside functions

def download_experiments_RNA_seq_nf_core(gse_list:list[str],root_storage_dir:str,output_dir:str, tracker:FileTracker, download_raw:bool=True, scan:bool=True,run_and_delete:bool=True):
    """
    Orchestrates the download and processing of RNA-Seq studies using nf-core/rnaseq.
    """
    #TODO: WE MIGHT NOT NEED THS ANYMORE
    PATH_TO_INDEX = f"{root_storage_dir}genome_index/tair10"
    PATH_TO_GTF = f"{root_storage_dir}genome_index/Arabidopsis_thaliana.TAIR10.56.gtf"
    
    processor = RNASeq_processor(threads=1, genome_index=PATH_TO_INDEX, gtf_annotation=PATH_TO_GTF,profile='singularity')
    tracker_save_path = os.path.join(output_dir, "rnaseq_tracker_stats.json")
    
    # if not os.path.exists(output_dir): os.makedirs(output_dir)
    valid_gse_ids = []

    for gse_id in tqdm(gse_list, desc="Processing RNA-Seq", unit="study"):
        
        if tracker.is_ignored(gse_id): continue
        if tracker.is_processed(gse_id):
            valid_gse_ids.append(gse_id); continue
            
        # Folders
        soft_path = os.path.join(output_dir, f"{gse_id}_family.soft.gz")
        fastq_folder = os.path.join(output_dir, "fastq_storage", gse_id)
        results_folder = os.path.join(output_dir, "processed_rnaseq", gse_id)
        cluster_temp = os.environ.get('TMPDIR', '/tmp')
        # temp_files = os.path.join(output_dir,'temp')
        # if not os.path.exists(temp_files): os.makedirs(temp_files)
        if not os.path.exists(results_folder): os.makedirs(results_folder)
        try:
            # 1. Metadata
            try:
                gse:GSM | GSE | GPL | GDS = GEOparse.get_GEO(geo=gse_id, destdir=output_dir, silent=True)
                if isinstance(gse,GSE):
                    pass
                else:
                    raise ValueError('Not the correct type of ID')
            except:
                tracker.mark_ignore(gse_id); continue

            # Tracker Update
            # platform = gse.metadata.get('platform_id', ['Unknown'])[0]
            # num_samples = len(gse.gsms)
            has_sra = check_metadata_for_sra_boolean(gse) # Ensure this helper is imported
            
            if scan:
                # tracker.update_platform(platform, num_samples, has_raw=has_sra)
                if has_sra: valid_gse_ids.append(gse_id)
                if os.path.exists(soft_path): os.remove(soft_path)
                continue
            
            if not has_sra:
                tracker.mark_ignore(gse_id)
                continue

            # 2. Pipeline Execution
            if download_raw:
                # tracker.update_platform(platform, num_samples, has_raw=True)
                
                # A. Download
                if not tracker.is_downloaded(gse_id):
                    try:
                        # processor.download_fastq(gse, fastq_folder,cluster_temp)
                        tracker.mark_downloaded(gse_id)
                        tracker.save_to_json(tracker_save_path)
                    except Exception as e:
                        tqdm.write(f"Download Error {gse_id}: {e}")
                        tracker.mark_ignore(gse_id)
                        if os.path.exists(fastq_folder): shutil.rmtree(fastq_folder)
                        continue
                if not os.path.exists(fastq_folder) or not os.listdir(fastq_folder):
                    print(f"FASTQs not found for {gse_id}, skipping...")
                    tracker.mark_ignore(gse_id)
                    # continue

                # 3. Process with Nextflow
                if not tracker.is_processed(gse_id):
                    try:
                        # success = processor.run_pipeline_on_study(gse_id, fastq_folder, results_folder)
                        success = True
                        if success:
                            # 4. Cleanup
                            if run_and_delete:
                                print(f"Cleaning raw FASTQs for {gse_id}...")
                                if os.path.exists(fastq_folder): shutil.rmtree(fastq_folder)
                                
                            tracker.mark_processed(gse_id)
                            tracker.save_to_json(tracker_save_path)
                            valid_gse_ids.append(gse_id)
                        else:
                            tracker.mark_ignore(gse_id)
                            
                    except Exception as e:
                        tqdm.write(f"Pipeline Error {gse_id}: {e}")
                        tracker.mark_ignore(gse_id)
                        continue

        except Exception as e:
            print(f"General Error {gse_id}: {e}")
            tracker.mark_ignore(gse_id)

    return valid_gse_ids