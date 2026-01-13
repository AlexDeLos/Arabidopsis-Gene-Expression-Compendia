import os
import pandas as pd
import numpy as np
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
import sys
import json

module_dir = './'
sys.path.append(module_dir)

# Import the processing function from your other file
from src.data_importing.helpers import process_metadata

class Microarray_tracker:
    def __init__(self) -> None:
        self.platform_counts: dict = {}
        self.totals: dict = {
            'total_studies_seen': 0,
            'total_sample_seen': 0,
            'total_samples_used': 0,
            'total_studies_used': 0
        }
        self.study_tracker:dict = {}

    def update_platform(self, platform: str, samples: int, used: bool):
        """
        Updates the tracker with a new study.
        """
        # 1. Ensure platform entry exists
        if platform not in self.platform_counts:
            self.platform_counts[platform] = {
                'studies_seen': 0,
                'samples_seen': 0,
                'studies_used': 0,
                'samples_used': 0
            }
        
        # 2. Update 'SEEN' counts
        self.totals['total_studies_seen'] += 1
        self.totals['total_sample_seen'] += samples
        
        self.platform_counts[platform]['studies_seen'] += 1
        self.platform_counts[platform]['samples_seen'] += samples

        # 3. Update 'USED' counts
        if used:
            self.totals['total_studies_used'] += 1
            self.totals['total_samples_used'] += samples
            
            self.platform_counts[platform]['studies_used'] += 1
            self.platform_counts[platform]['samples_used'] += samples

    def classify_study(self,study_id,classification):
        if study_id in self.study_tracker.keys():
            raise ValueError('Study already classified')
        self.study_tracker[study_id] = classification

    def save_to_json(self, filename="tracker_results.json"):
        """Saves current state to a JSON file."""
        data = {
            "totals": self.totals,
            "platform_counts": self.platform_counts,
            "study_tracker": self.study_tracker
        }
        with open(filename, 'w') as f:
            json.dump(data, f, indent=4)
        print(f"Tracker data saved to {filename}")

    @classmethod
    def load_from_json(cls, filename="tracker_results.json"):
        """Creates a new tracker instance loaded from a JSON file."""
        if not os.path.exists(filename):
            print(f"File {filename} not found. Returning empty tracker.")
            return cls()
            
        with open(filename, 'r') as f:
            data = json.load(f)
        
        tracker = cls()
        tracker.totals = data.get("totals", tracker.totals)
        tracker.platform_counts = data.get("platform_counts", {})
        tracker.study_tracker = data.get("study_tracker", {})
        return tracker

    def print_summary(self):
        print("\n" + "="*40)
        print("       MICROARRAY PROCESSING SUMMARY")
        print("="*40)
        print(f"Total Studies: {self.totals['total_studies_seen']} Seen | {self.totals['total_studies_used']} Used")
        print(f"Total Samples: {self.totals['total_sample_seen']} Seen | {self.totals['total_samples_used']} Used")
        print("-" * 40)
        print("BREAKDOWN BY PLATFORM:")
        print(f"{'Platform':<15} | {'St(Seen)':<8} {'Sa(Seen)':<8} | {'St(Used)':<8} {'Sa(Used)':<8}")
        print("-" * 60)
        for plat, stats in self.platform_counts.items():
            print(f"{plat:<15} | {stats['studies_seen']:<8} {stats['samples_seen']:<8} | {stats['studies_used']:<8} {stats['samples_used']:<8}")
        print("="*40 + "\n")

class RNASeq_tracker:
    def __init__(self) -> None:
        self.platform_counts: dict = {}
        self.totals: dict = {
            'total_studies_seen': 0,
            'total_sample_seen': 0,
            'total_samples_used': 0, # Used here means "Has SRA/Raw Data"
            'total_studies_used': 0
        }

    def update_platform(self, platform: str, samples: int, has_raw: bool):
        """
        Updates the tracker. 
        'has_raw' = True implies the study has accessible SRA data (FASTQ).
        """
        # 1. Ensure platform entry exists
        if platform not in self.platform_counts:
            self.platform_counts[platform] = {
                'studies_seen': 0, 'samples_seen': 0,
                'studies_with_raw': 0, 'samples_with_raw': 0
            }
        
        # 2. Update 'SEEN' counts
        self.totals['total_studies_seen'] += 1
        self.totals['total_sample_seen'] += samples
        self.platform_counts[platform]['studies_seen'] += 1
        self.platform_counts[platform]['samples_seen'] += samples

        # 3. Update 'USED' (Has Raw) counts
        if has_raw:
            self.totals['total_studies_used'] += 1
            self.totals['total_samples_used'] += samples
            self.platform_counts[platform]['studies_with_raw'] += 1
            self.platform_counts[platform]['samples_with_raw'] += samples

    def save_to_json(self, filename="rnaseq_tracker_results.json"):
        with open(filename, 'w') as f:
            json.dump({"totals": self.totals, "platform_counts": self.platform_counts}, f, indent=4)
        print(f"RNA-Seq Tracker data saved to {filename}")

    def print_summary(self):
        print("\n" + "="*40)
        print("       RNA-SEQ PROCESSING SUMMARY")
        print("="*40)
        print(f"Total Studies: {self.totals['total_studies_seen']} Seen | {self.totals['total_studies_used']} Has Raw (SRA)")
        print(f"Total Samples: {self.totals['total_sample_seen']} Seen | {self.totals['total_samples_used']} Has Raw (SRA)")
        print("-" * 60)
        print(f"{'Platform':<15} | {'St(Seen)':<8} {'Sa(Seen)':<8} | {'St(Raw)':<8} {'Sa(Raw)':<8}")
        print("-" * 60)
        for plat, stats in self.platform_counts.items():
            if stats['studies_seen'] > 0:
                print(f"{plat:<15} | {stats['studies_seen']:<8} {stats['samples_seen']:<8} | {stats['studies_with_raw']:<8} {stats['samples_with_raw']:<8}")
        print("="*40 + "\n")
    
    @classmethod
    def load_from_json(cls, filename="rnaseq_tracker_results.json"):
        """Creates a new tracker instance loaded from a JSON file."""
        if not os.path.exists(filename):
            print(f"File {filename} not found. Returning empty tracker.")
            return cls()
            
        with open(filename, 'r') as f:
            data = json.load(f)
        
        tracker = cls()
        tracker.totals = data.get("totals", tracker.totals)
        tracker.platform_counts = data.get("platform_counts", {})
        return tracker


class Data_processing:
    def __init__(self):
        
        # --- CHECK R DEPENDENCIES ---
        try:
            self.oligo = importr('oligo')
            self.base = importr('base')
            self.biobase = importr('Biobase')
        except Exception as e:
            print(f"Error: R packages (oligo, Biobase) not installed.:{e}")
            return
        
    def _install_missing_r_package(self, package_name):
        """
        Attempts to install a missing Bioconductor package on the fly.
        """
        print(f"    > ATTEMPTING TO AUTO-INSTALL: {package_name}...")
        try:
            # Import R utils
            utils = importr('utils')
            utils.chooseCRANmirror(ind=1)
            
            # Use BiocManager to install
            biocmanager = importr('BiocManager')
            biocmanager.install(ro.StrVector([package_name]), ask=False, update=False)
            
            # Check if it worked
            if rpackages.isinstalled(package_name):
                print(f"    > SUCCESS: {package_name} installed.")
                return True
            else:
                print(f"    > FAILURE: Could not install {package_name}.")
                return False
        except Exception as e:
            print(f"    > FAILURE during installation: {e}")
            return False
    

    def process_microarray_metadata_and_rma(self, experiment_id, raw_data_folder, output_folder, gse):
        """
        1. Recursively searches for .CEL files.
        2. Uses R 'oligo' to read and RMA-normalize them.
        3. Maps Probe IDs to Gene Symbols using the GSE platform data.
        4. Exports a Pandas DataFrame (Rows=Genes, Cols=Samples).
        """
        print(f"\n[RMA Pipeline] Processing {experiment_id}...")

        # --- SAVE ---
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
            
        # --- RECURSIVE FILE DISCOVERY ---
        cel_files = []
        for root, dirs, files in os.walk(raw_data_folder):
            for f in files:
                if f.lower().endswith('.cel') or f.lower().endswith('.cel.gz'):
                    full_path = os.path.join(root, f)
                    cel_files.append(full_path)
                        
        if not cel_files:
            print(f"  - ABORTING: No CEL files found in {raw_data_folder}.")
            return

        print(f"  - Found {len(cel_files)} CEL files.")

        # --- METADATA PROCESSING ---
        print("  - Processing metadata for all samples in study...")
        try:
            # process_metadata_study(experiment_id, gse, save_path=output_folder)
            for gsm_name, gsm in gse.gsms.items():
                # Passing save_path as requested
                process_metadata(experiment_id, gse, gsm, save_path=output_folder)
        except Exception as e:
            print(f"    > Metadata processing failed: {e}")

        # --- R MAGIC (RMA NORMALIZATION) ---
        try:
            cel_files_r = ro.StrVector(cel_files)
            
            print("  - Reading CEL files (oligo)...")
            try:
                raw_data = self.oligo.read_celfiles(cel_files_r)
            except Exception as e_r:
                error_msg = str(e_r)
                
                # Check for the specific "package not loaded" error
                if "could not be loaded" in error_msg and "pd." in error_msg:
                    import re
                    # Extract package name (e.g., pd.atdschip.expr)
                    match = re.search(r"(pd\.[a-zA-Z0-9\.]+)", error_msg)
                    
                    if match:
                        missing_pkg = match.group(1)
                        # CALL THE HELPER FUNCTION TO INSTALL
                        installed = self._install_missing_r_package(missing_pkg)
                        
                        if installed:
                            print("    > Retrying read_celfiles...")
                            try:
                                # RETRY READING THE FILES
                                raw_data = self.oligo.read_celfiles(cel_files_r)
                            except Exception as e_retry:
                                print(f"  - SKIPPING {experiment_id}: Still failed after install. {e_retry}")
                                raise e_r
                        else:
                            print(f"  - SKIPPING {experiment_id}: Cannot install required package {missing_pkg}.")
                            raise e_r
                    else:
                        print(f"  - SKIPPING {experiment_id}: Annotation package error (Unknown name).")
                        raise e_r
                else:
                    raise e_r

            print("  - Running RMA normalization...")
            eset = self.oligo.rma(raw_data)
            
            # Extract expression matrix
            exprs_matrix_r = self.biobase.exprs(eset)
            
            # --- FIX STEP 1: Handle Data & Index Separately ---
            data_np = np.array(exprs_matrix_r)
            
            # 2. Extract Row Names (Probe IDs) explicitly from the R object
            try:
                probe_ids = list(exprs_matrix_r.rownames)
            except:
                # Fallback: Ask the ExpressionSet for feature names directly
                probe_ids = list(self.biobase.featureNames(eset))
                
            # 3. Create DataFrame and Assign the Index
            df = pd.DataFrame(data_np)
            df.index = probe_ids  
            df.index.name = "ProbeID"
            
            # --- CLEAN COLUMN NAMES (Samples) ---
            clean_names = []
            for c in cel_files:
                filename = os.path.basename(c)
                if filename.lower().endswith('.gz'):
                    filename = filename[:-3]
                if filename.lower().endswith('.cel'):
                    filename = filename[:-4]
                clean_names.append(filename.split('_')[0])
                
            df.columns = clean_names

            # --- ANNOTATION: PROBE -> GENE ---
            print("  - Mapping Probe IDs to Genes...")
            try:
                # 1. Get Platform ID (GPL)
                platform_id = gse.metadata.get('platform_id', [''])[0]
                
                if platform_id and platform_id in gse.gpls:
                    gpl = gse.gpls[platform_id]
                    annot = gpl.table
                    
                    # 2. Find Gene Column (Heuristic)
                    candidates = ['Gene Symbol', 'GENE_SYMBOL', 'Gene_Symbol', 'SYMBOL', 'ORF', 'GeneSymbol']
                    gene_col = None
                    for col in candidates:
                        if col in annot.columns:
                            gene_col = col
                            break
                    
                    if gene_col:
                        print(f"    > Found gene column: '{gene_col}'")
                        
                        # 3. Create Mapping Dictionary
                        def clean_id(val):
                            try:
                                return str(int(float(val)))
                            except:
                                return str(val)

                        annot['ID'] = annot['ID'].apply(clean_id)
                        
                        # Clean the gene symbol
                        clean_genes = annot[gene_col].astype(str).apply(lambda x: x.split('///')[0].split('//')[0].strip())
                        
                        probe_map = dict(zip(annot['ID'], clean_genes))
                        
                        # 4. Map the Index
                        df.index = df.index.map(str)
                        df['GeneSymbol'] = df.index.map(probe_map)
                        df['GeneSymbol'] = df['GeneSymbol'].fillna(df.index.to_series())
                        
                        # 5. Aggregate Duplicates
                        df = df.groupby('GeneSymbol').mean()
                        
                        print(f"    > Annotation complete. Collapsed {len(probe_map)} probes into {len(df)} genes.")
                    else:
                        print(f"    > WARNING: Could not find Gene Symbol column in GPL ({platform_id}). Keeping Probe IDs.")
                else:
                    print("    > Platform data not found in GSE object. Keeping Probe IDs.")

            except Exception as e_annot:
                print(f"    > Annotation step failed: {e_annot}. Keeping Probe IDs.")
                raise e_annot
                
            output_file = os.path.join(output_folder, f"{experiment_id}_RMA_Genes.csv")
            df.to_csv(output_file)
            print(f"  - SUCCESS: Saved to {output_file}")
            
        except Exception as e:
            print(f"  - RMA Failed for {experiment_id}: {e}")
            raise e