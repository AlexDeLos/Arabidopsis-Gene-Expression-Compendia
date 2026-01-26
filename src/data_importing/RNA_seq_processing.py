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
from src.constants import *

class RNASeq_tracker:
    # ... (Keep the exact same Tracker class from previous iterations) ...
    def __init__(self) -> None:
        self.platform_counts: dict = {}
        self.totals: dict = {
            'total_studies_seen': 0, 'total_sample_seen': 0,
            'total_samples_used': 0, 'total_studies_used': 0
        }
        self.states: dict = {'ignore': set(), 'downloaded': set(), 'processed': set()}

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

    def mark_ignore(self, gse_id):
        self.states['ignore'].add(gse_id); self.states['downloaded'].discard(gse_id); self.states['processed'].discard(gse_id)
    def mark_downloaded(self, gse_id):
        self.states['downloaded'].add(gse_id); self.states['ignore'].discard(gse_id); self.states['processed'].discard(gse_id)
    def mark_processed(self, gse_id):
        self.states['processed'].add(gse_id); self.states['downloaded'].discard(gse_id); self.states['ignore'].discard(gse_id)
    def is_ignored(self, gse_id): return gse_id in self.states['ignore']
    def is_processed(self, gse_id): return gse_id in self.states['processed']
    def is_downloaded(self, gse_id): return gse_id in self.states['downloaded']
    
    def save_to_json(self, filename="rnaseq_tracker_results.json"):
        serializable_states = {k: list(v) for k, v in self.states.items()}
        data = {"totals": self.totals, "platform_counts": self.platform_counts, "states": serializable_states}
        with open(filename, 'w') as f: json.dump(data, f, indent=4)
    
    @classmethod
    def load_from_json(cls, filename="rnaseq_tracker_results.json"):
        if not os.path.exists(filename): return cls()
        with open(filename, 'r') as f: data = json.load(f)
        tracker = cls()
        tracker.totals = data.get("totals", tracker.totals)
        tracker.platform_counts = data.get("platform_counts", {})
        loaded_states = data.get("states", {})
        tracker.states = {k: set(loaded_states.get(k, [])) for k in ['ignore', 'downloaded', 'processed']}
        return tracker


class RNASeq_processor:
    def __init__(self, threads=4, genome_index=None, gtf_annotation=None):
        self.threads = str(threads)
        self.genome_index = genome_index # Path to HISAT2 index prefix
        self.gtf_annotation = gtf_annotation # Path to .gtf file
        
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
        """Downloads FASTQ files."""
        if not os.path.exists(output_folder): os.makedirs(output_folder)
        
        # Ensure temp directory exists (Critical for fasterq-dump)
        if not os.path.exists(temp_files): os.makedirs(temp_files)

        sra_map = {gsm: self.get_srr_ids(gsm) for gsm in gse.gsms.keys()}
        
        # Filter empty
        sra_map = {k:v for k,v in sra_map.items() if v}
        
        # --- FIX 1: Point to the correct binary name ---
        # Was: .../bin/sratools.3.0.10
        # Now: .../bin/fasterq-dump
        fasterq_path = '/tudelft.net/staff-umbrella/AT GE Datasets/sratoolkit.3.0.10-ubuntu64/bin/fasterq-dump'
        # FASTERQ_PATH = "/tudelft.net/staff-umbrella/AT GE Datasets/sratoolkit.3.0.10-ubuntu64/bin/fasterq-dump"
        for gsm, srrs in tqdm(sra_map.items(), desc="Downloading SRRs", leave=False):
            for srr in srrs:
                # Check existence
                if any(f.startswith(srr) for f in os.listdir(output_folder)): continue
                
                if CLUSTER_RUN:
                    cmd = [fasterq_path, "--split-files", "--outdir", output_folder, "--temp", temp_files , "--threads", self.threads, srr]
                else:
                    cmd = ["fasterq-dump", "--split-files", "--outdir", output_folder, "--temp", temp_files , "--threads", self.threads, srr]
                
                # --- FIX 2: Remove DEVNULL so you can see errors in the log ---
                # stdout=subprocess.DEVNULL is fine, but keep stderr visible
                try:
                    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
                except subprocess.CalledProcessError as e:
                    print(f"\nCRITICAL ERROR downloading {srr}:")
                    print(f"Command tried: {' '.join(cmd)}")
                    # This will print the actual error from the tool (e.g. 'disk full' or 'not found')
                    raise e

    def _find_adapter_file(self, filename="TruSeq3-PE.fa"):
        """
        Dynamically finds the absolute path of the adapter file 
        within the current Conda environment.
        """
        conda_prefix = os.environ.get("CONDA_PREFIX")
        
        if not conda_prefix:
            raise EnvironmentError("CONDA_PREFIX not found. Are you running this inside a Conda environment?")

        # Construct the find command: find $CONDA_PREFIX -name filename
        try:
            # check_output runs the command and returns the stdout result as bytes
            # text=True converts bytes to string
            output = subprocess.check_output(
                ["find", conda_prefix, "-name", filename], 
                text=True
            ).strip()
            
            # If multiple are found, 'output' might contain newlines. Take the first one.
            if "\n" in output:
                output = output.split("\n")[0]
                
            if not output or not os.path.exists(output):
                raise FileNotFoundError(f"Could not find {filename} in {conda_prefix}")
                
            return output
            
        except subprocess.CalledProcessError as e:
            raise FileNotFoundError(f"Error searching for adapters: {e}")
        
    def run_pipeline_on_study(self, gse_id, fastq_folder, output_folder):
        """
        Executes the Bio-protocol PDF workflow:
        Trim (Trimmomatic) -> Align (HISAT2) -> Sort (Samtools) -> Count (featureCounts)
        """
        print(f"Running pipeline on study {gse_id}")
        if not os.path.exists(output_folder): os.makedirs(output_folder)
        
        # 1. Identify Samples
        files = [f for f in os.listdir(fastq_folder) if f.endswith('.fastq')]
        samples = {} # SRR -> [file1, file2]
        for f in files:
            srr = f.split('_')[0].split('.')[0]
            if srr not in samples: samples[srr] = []
            samples[srr].append(os.path.join(fastq_folder, f))
            samples[srr].sort()

        bam_files = []

        # PROCESS PER SAMPLE (Trim & Align)
        for srr, fqs in tqdm(samples.items(), desc=f"Aligning {gse_id}", leave=False):
            
            # Define Filenames
            trimmed_1 = os.path.join(output_folder, f"{srr}_1_paired.fq")
            trimmed_2 = os.path.join(output_folder, f"{srr}_2_paired.fq")
            # Unpaired output (usually discarded in standard pipelines)
            unpaired_1 = os.path.join(output_folder, f"{srr}_1_unpaired.fq") 
            unpaired_2 = os.path.join(output_folder, f"{srr}_2_unpaired.fq")
            
            bam_out = os.path.join(output_folder, f"{srr}.sorted.bam")
            
            # --- STEP 1: TRIMMOMATIC (Cleaning) ---
            # Command structure based on Bio-protocol PDF recommendations
            adapter_path = self._find_adapter_file("TruSeq3-PE.fa")

            # 2. Build the command
            is_paired = len(fqs) == 2
            
            trim_cmd = ["trimmomatic", "PE" if is_paired else "SE", "-threads", self.threads]
            
            if is_paired:
                trim_cmd.extend([fqs[0], fqs[1], trimmed_1, unpaired_1, trimmed_2, unpaired_2])
            else:
                trim_cmd.extend([fqs[0], trimmed_1]) 

            # 3. Add the Clipping parameter using the dynamic path
            # Note the f"..." string formatting
            trim_cmd.extend([
                f"ILLUMINACLIP:{adapter_path}:2:30:10", 
                "LEADING:3", 
                "TRAILING:3", 
                "SLIDINGWINDOW:4:15", 
                "MINLEN:36"
            ])
            
            subprocess.run(trim_cmd, check=True, stderr=subprocess.DEVNULL)            
            # --- STEP 2: HISAT2 (Alignment) ---
            # We pipe HISAT2 directly to Samtools Sort to avoid huge SAM files
            
            
            hisat_cmd = ["hisat2", "-p", self.threads, "-x", self.genome_index]
            if is_paired:
                hisat_cmd.extend(["-1", trimmed_1, "-2", trimmed_2])
            else:
                hisat_cmd.extend(["-U", trimmed_1])
                
            # Pipeline: HISAT2 -> Samtools View (BAM) -> Samtools Sort
            p1 = subprocess.Popen(hisat_cmd, stdout=subprocess.PIPE)
            p2 = subprocess.Popen(["samtools", "sort", "-@", self.threads, "-o", bam_out], stdin=p1.stdout)
            p1.stdout.close()
            p2.communicate()
            
            bam_files.append(bam_out)

            # CLEANUP INTERMEDIATE FASTQS (Crucial for space)
            for f in [trimmed_1, trimmed_2, unpaired_1, unpaired_2]:
                if os.path.exists(f): os.remove(f)

        # --- STEP 3: FEATURECOUNTS (Quantification) ---
        
        
        count_file = os.path.join(output_folder, f"{gse_id}_counts.txt")
        
        fc_cmd = [
            "featureCounts", 
            "-T", self.threads, 
            "-a", self.gtf_annotation, 
            "-o", count_file,
            "-p", # Count fragments instead of reads (for paired-end)
            "-t", "exon", # Feature type
            "-g", "gene_id" # Attribute type
        ]
        fc_cmd.extend(bam_files)
        
        print(f"  - Counting features for {len(bam_files)} BAMs...")
        subprocess.run(fc_cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        # CLEANUP BAMS (We only need the count matrix)
        for b in bam_files:
            if os.path.exists(b): os.remove(b)

        # Clean up the featureCounts format to match Microarray style (GeneID index)
        # featureCounts output has a header and metadata columns we don't need
        df = pd.read_csv(count_file, sep="\t", comment="#")
        # Keep Geneid and sample columns (samples start from col 7)
        df = df.set_index("Geneid")
        df = df.iloc[:, 5:] # Drop Chr, Start, End, Strand, Length
        
        # Rename columns (currently paths like 'folder/SRR123.sorted.bam') to just SRR
        df.columns = [c.split('/')[-1].split('.')[0] for c in df.columns]
        
        final_csv = os.path.join(output_folder, f"{gse_id}_counts.csv")
        df.to_csv(final_csv)
        
        # Remove original raw output
        os.remove(count_file)
        os.remove(count_file + ".summary")
        
        return True
    
def download_experiments_RNA_seq(gse_list,root_storage_dir ,output_dir, tracker, download_raw=True, scan=True,run_and_delete:bool=True):
    
    PATH_TO_INDEX = f"{root_storage_dir}genome_index/tair10"
    PATH_TO_GTF = f"{root_storage_dir}genome_index/Arabidopsis_thaliana.TAIR10.56.gtf"
    
    processor = RNASeq_processor(threads=8, genome_index=PATH_TO_INDEX, gtf_annotation=PATH_TO_GTF)
    tracker_save_path = os.path.join(output_dir, "rnaseq_tracker_stats.json")
    
    if not os.path.exists(output_dir): os.makedirs(output_dir)
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
                gse = GEOparse.get_GEO(geo=gse_id, destdir=output_dir, silent=True)
            except:
                tracker.mark_ignore(gse_id); continue

            # Tracker Update
            platform = gse.metadata.get('platform_id', ['Unknown'])[0]
            num_samples = len(gse.gsms)
            has_sra = check_metadata_for_sra(gse) # Ensure this helper is imported
            
            if scan:
                tracker.update_platform(platform, num_samples, has_raw=has_sra)
                if has_sra: valid_gse_ids.append(gse_id)
                if os.path.exists(soft_path): os.remove(soft_path)
                continue
            
            if not has_sra:
                tracker.mark_ignore(gse_id)
                continue

            # 2. Pipeline Execution
            if download_raw:
                tracker.update_platform(platform, num_samples, has_raw=True)
                
                # A. Download
                if not tracker.is_downloaded(gse_id):
                    try:
                        processor.download_fastq(gse, fastq_folder,cluster_temp)
                        tracker.mark_downloaded(gse_id)
                        tracker.save_to_json(tracker_save_path)
                    except Exception as e:
                        tqdm.write(f"Download Error {gse_id}: {e}")
                        tracker.mark_ignore(gse_id)
                        if os.path.exists(fastq_folder): shutil.rmtree(fastq_folder)
                        continue

                # B. Process (Trim -> Align -> Count)
                if not tracker.is_processed(gse_id):
                    try:
                        success = processor.run_pipeline_on_study(gse_id, fastq_folder, results_folder)
                        
                        if success:
                            # C. DELETE FASTQS (Critical Step)
                            if run_and_delete:
                                if os.path.exists(fastq_folder): shutil.rmtree(fastq_folder)
                                if os.path.exists(soft_path): os.remove(soft_path)
                            
                            tracker.mark_processed(gse_id)
                            tracker.save_to_json(tracker_save_path)
                            valid_gse_ids.append(gse_id)
                            
                    except Exception as e:
                        tqdm.write(f"Pipeline Error {gse_id}: {e}")
                        # Keep FASTQ in case of error for debugging, or uncomment next line to force clean
                        # if os.path.exists(fastq_folder): shutil.rmtree(fastq_folder)
                        tracker.mark_ignore(gse_id)
                        continue

        except Exception as e:
            tracker.mark_ignore(gse_id)

    return valid_gse_ids

