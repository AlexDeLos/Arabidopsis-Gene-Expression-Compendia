import os
import pandas as pd
import numpy as np
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
import rpy2.robjects.packages as rpackages
from rpy2.robjects.vectors import StrVector
import sys
import json
import GEOparse
from tqdm import tqdm
import re
import mygene
module_dir = './'
sys.path.append(module_dir)

# Import the processing function from your other file
from src.data_importing.helpers.helpers import process_metadata
from src.data_importing.helpers.download_helper import clean_files,check_metadata_for_cel, count_cel_files


LOCUS_REGEX = re.compile(r'^AT[1-5MC]G\d{5}$', re.IGNORECASE)
class Microarray_tracker:
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

    def update_platform(self, platform: str, samples: int, used: bool):
        """
        Updates the statistical counters for platforms.
        """
        if platform not in self.platform_counts:
            self.platform_counts[platform] = {
                'studies_seen': 0,
                'studies_used': 0,
                'samples_seen': 0,
                'samples_used': 0
            }
        
        self.totals['total_studies_seen'] += 1
        self.totals['total_sample_seen'] += samples
        self.platform_counts[platform]['studies_seen'] += 1
        self.platform_counts[platform]['samples_seen'] += samples

        if used:
            self.totals['total_studies_used'] += 1
            self.totals['total_samples_used'] += samples
            self.platform_counts[platform]['studies_used'] += 1
            self.platform_counts[platform]['samples_used'] += samples

    # --- STATE MANAGEMENT METHODS ---
    def sync_with_filesystem(self, download_dir: str, processed_dir: str):
        """
        Scans the actual directories on disk and updates the tracker states 
        to match reality.
        
        Logic:
        1. Files on disk always override the current memory state.
        2. If a processed CSV exists (checked inside subfolders) -> State = PROCESSED.
           Structure checked: processed_dir/GSExxx/GSExxx_RMA_Genes.csv
        3. If a raw folder exists (and no processed CSV) -> State = DOWNLOADED.
        4. If it was in 'ignore' but files exist, it is removed from 'ignore'.
        """
        print(f"Syncing tracker with files...")
        print(f"  - Scannning Processed: {processed_dir}")
        print(f"  - Scannning Downloads: {download_dir}")

        # 1. Identify what is actually processed on disk
        actual_processed_ids = set()
        if os.path.exists(processed_dir):
            # Iterate over folders in the root processed directory
            for gse_id in os.listdir(processed_dir):
                subfolder_path = os.path.join(processed_dir, gse_id)
                
                # Ensure it is a directory (e.g., GSE12345)
                if os.path.isdir(subfolder_path):
                    # Construct the expected file path inside that subfolder
                    expected_csv = os.path.join(subfolder_path, f"{gse_id}_RMA_Genes.csv")
                    
                    if os.path.exists(expected_csv):
                        actual_processed_ids.add(gse_id)

        # 2. Identify what is actually downloaded on disk
        actual_downloaded_ids = set()
        if os.path.exists(download_dir):
            for d in os.listdir(download_dir):
                full_path = os.path.join(download_dir, d)
                # Check for raw folders starting with GSE
                if d.startswith("GSE") and os.path.isdir(full_path):
                    actual_downloaded_ids.add(d)

        # 3. Update Tracker State based on Disk Reality
        changes_made = 0

        # Sync Processed
        for gse in actual_processed_ids:
            if not self.is_processed(gse):
                self.mark_processed(gse) # Helper handles mutual exclusivity
                changes_made += 1

        # Sync Downloaded
        # (Only mark as downloaded if it is NOT already marked as processed)
        for gse in actual_downloaded_ids:
            if gse not in actual_processed_ids:
                if not self.is_downloaded(gse):
                    self.mark_downloaded(gse) # Helper handles mutual exclusivity
                    changes_made += 1
        
        print(f"  - Sync complete. Updated {changes_made} records based on file existence.")

    def compare_states(self, other_tracker) -> bool:
        """
        Compares the state sets (ignore, downloaded, processed) of this tracker 
        against another tracker.
        
        Returns:
            True: If the sets in both trackers are identical.
            False: If there is any discrepancy in the sets.
        """
        # We compare the 'states' dictionary which contains the 3 sets.
        # Sets are unordered, so '==' checks for content equality, which is what we want.
        
        if not isinstance(other_tracker, Microarray_tracker):
            print("Warning: Comparing Microarray_tracker with a different object type.")
            return False

        are_equal = self.states == other_tracker.states
        
        if not are_equal:
            # Optional: Debug print to show where they differ
            keys = ['downloaded', 'processed']
            for k in keys:
                if self.states[k] != other_tracker.states[k]:
                    diff_1 = self.states[k] - other_tracker.states[k]
                    diff_2 = other_tracker.states[k] - self.states[k]
                    if diff_1: print(f"  Difference in '{k}': Self has extra {len(diff_1)} items (e.g., {list(diff_1)[:3]}...)")
                    if diff_2: print(f"  Difference in '{k}': Other has extra {len(diff_2)} items (e.g., {list(diff_2)[:3]}...)")
        
        return are_equal
    
    def mark_ignore(self, gse_id):
        """Moves ID to ignore set; removes from others."""
        self.states['ignore'].add(gse_id)
        self.states['downloaded'].discard(gse_id)
        self.states['processed'].discard(gse_id)

    def mark_downloaded(self, gse_id):
        """Moves ID to downloaded set; removes from others."""
        self.states['downloaded'].add(gse_id)
        self.states['ignore'].discard(gse_id)
        self.states['processed'].discard(gse_id)

    def mark_processed(self, gse_id):
        """Moves ID to processed set; removes from others."""
        self.states['processed'].add(gse_id)
        self.states['downloaded'].discard(gse_id)
        self.states['ignore'].discard(gse_id)

    def is_ignored(self, gse_id):
        return gse_id in self.states['ignore']

    def is_processed(self, gse_id):
        return gse_id in self.states['processed']
    
    def is_downloaded(self, gse_id):
        return gse_id in self.states['downloaded']

    def save_to_json(self, filename="tracker_results.json"):
        """Saves current state to a JSON file. Converts sets to lists for JSON serialization."""
        # Convert sets to lists for JSON serialization
        serializable_states = {
            k: list(v) for k, v in self.states.items()
        }
        
        data = {
            "totals": self.totals,
            "platform_counts": self.platform_counts,
            "study_tracker": self.study_tracker,
            "states": serializable_states
        }
        with open(filename, 'w') as f:
            json.dump(data, f, indent=4)
        # print(f"Tracker saved to {filename}") # Optional: Comment out to reduce console spam

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
        
        # Load states and convert lists back to sets
        loaded_states = data.get("states", {})
        tracker.states = {
            'ignore': set(loaded_states.get('ignore', [])),
            'downloaded': set(loaded_states.get('downloaded', [])),
            'processed': set(loaded_states.get('processed', []))
        }
        
        return tracker

    def print_summary(self):
        print("\n" + "="*40)
        print("       MICROARRAY PROCESSING SUMMARY")
        print("="*40)
        print(f"Total Studies: {self.totals['total_studies_seen']} Seen | {self.totals['total_studies_used']} Used")
        print(f"Ignored: {len(self.states['ignore'])} | Downloaded: {len(self.states['downloaded'])} | Processed: {len(self.states['processed'])}")
        print("-" * 40)
        print("BREAKDOWN BY PLATFORM:")
        print(f"{'Platform':<15} | {'St(Seen)':<8} {'Sa(Seen)':<8} | {'St(Used)':<8} {'Sa(Used)':<8}")
        print("-" * 60)
        for plat, stats in self.platform_counts.items():
            print(f"{plat:<15} | {stats['studies_seen']:<8} {stats['samples_seen']:<8} | {stats['studies_used']:<8} {stats['samples_used']:<8}")
        print("="*40 + "\n")



class Microarray_data_processing:
    def __init__(self):
        # Initialize R environment immediately
        self._initialize_r_environment()
        
    def _initialize_r_environment(self):
        """
        Installs (if missing) and loads CORE R packages at startup.
        Assigns them to self for reuse.
        """
        print("\n[Init] Setting up R environment...")
        
        # 1. Select a mirror to avoid interactive prompts
        utils = rpackages.importr('utils')
        utils.chooseCRANmirror(ind=1) 

        # 2. Define Core Packages (The tools needed for every file)
        # We DO NOT include specific chip annotations here (e.g. pd.at...) 
        # as those must still be handled dynamically.
        core_packages = ['BiocManager', 'oligo', 'Biobase']

        # 3. Check and Install Missing Core Packages
        names_to_install = [x for x in core_packages if not rpackages.isinstalled(x)]
        if len(names_to_install) > 0:
            print(f"  - Installing missing core packages: {names_to_install}")
            
            # Ensure BiocManager is ready first
            if not rpackages.isinstalled('BiocManager'):
                utils.install_packages(StrVector(['BiocManager']))
            
            bioc_manager = rpackages.importr('BiocManager')
            bioc_manager.install(StrVector(names_to_install), update=False, ask=False)

        # 4. Load Packages and Attach to Self
        # This fixes the "object has no attribute 'oligo'" error
        try:
            print("  - Loading R libraries...")
            self.oligo = rpackages.importr('oligo')
            self.biobase = rpackages.importr('Biobase')
            
            # Keep a reference to utilities we might need later
            self.r_utils = utils
            # self.bioc_manager = rpackages.importr('BiocManager')
            
            print("  - [Success] R environment ready.")
            
        except Exception as e:
            print(f"  - [CRITICAL] Failed to load R packages: {e}")
            raise e

    # def _install_missing_r_package(self, package_name):
    #     """
    #     Helper for DYNAMIC annotation packages (e.g. pd.hugene...)
    #     Only used when a specific chip type is encountered.
    #     """
    #     try:
    #         print(f"    > Attempting to install missing annotation package: {package_name}...")
    #         self.bioc_manager.install(StrVector([package_name]), update=False, ask=False)
    #         return True
    #     except Exception as e:
    #         print(f"    > Failed to install {package_name}: {e}")
    #         return False

    def _map_symbols_to_locus_ids(self, identifiers):
        """
        Helper: Takes a list of strings (mix of Locus IDs and Symbols).
        Returns a dictionary { 'original_string': 'AT_LOCUS_ID' }.
        Drops items that cannot be resolved.
        """
        mg = mygene.MyGeneInfo()
        mapping = {}
        to_query = []

        # 1. First pass: Identify what is already a Locus ID
        for item in identifiers:
            item_clean = str(item).strip().upper()
            if LOCUS_REGEX.match(item_clean):
                mapping[item] = item_clean
            else:
                to_query.append(item)

        # 2. Second pass: Query external library for the rest
        if to_query:
            print(f"    > Querying MyGene.info for {len(to_query)} non-standard IDs...")
            # Query for symbol, alias, or name and ask for the 'locus_tag' or 'ensembl.gene'
            try:
                results = mg.querymany(to_query, scopes='symbol,alias,name,reporter', 
                                       fields='locus_tag,ensembl.gene', species='3702', 
                                       returnall=False, as_dataframe=False, verbose=False)
                
                for res in results:
                    query = res['query']
                    # Check if we got a hit
                    locus = None
                    if 'locus_tag' in res:
                        locus = res['locus_tag']
                    elif 'ensembl' in res:
                        # Sometimes ensembl is a list, sometimes a dict
                        ens = res['ensembl']
                        if isinstance(ens, list) and len(ens) > 0:
                            locus = ens[0].get('gene')
                        elif isinstance(ens, dict):
                            locus = ens.get('gene')
                    
                    # If we found a locus and it looks valid, save it
                    if locus and LOCUS_REGEX.match(locus):
                        mapping[query] = locus.upper()
                        
            except Exception as e:
                print(f"    > Warning: MyGene query failed ({e}). Dropping these IDs.")

        return mapping
    
    def process_microarray_metadata_and_rma(self, experiment_id, raw_data_folder, output_folder, gse) -> int:
        print(f"\n[RMA Pipeline] Processing {experiment_id}...")
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
            
        # --- 1. File Discovery ---
        cel_files = []
        for root, _, files in os.walk(raw_data_folder):
            for f in files:
                if f.lower().endswith('.cel') or f.lower().endswith('.cel.gz'):
                    cel_files.append(os.path.join(root, f))
        
        if not cel_files:
            print(f"  - ABORTING: No CEL files found in {raw_data_folder}.")
            raise ValueError('No cel files found')

        print(f"  - Found {len(cel_files)} CEL files.")

        # --- 2. Metadata Processing ---
        # (Kept your original logic here)
        print("  - Processing metadata for all samples in study...")
        try:
            for _, gsm in gse.gsms.items():
                process_metadata(experiment_id, gse, gsm, save_path=output_folder)
        except Exception as e:
            print(f"    > Metadata processing failed: {e}")

        # --- 3. RMA Normalization ---
        try:
            # RPY2 Imports assumed to be available in 'self' or globally
            cel_files_r = ro.StrVector(cel_files) # adjusted self.ro based on context
            print("  - Reading CEL files (oligo)...")
            
            # (Kept your original Package installation retry logic)
            try:
                raw_data = self.oligo.read_celfiles(cel_files_r)
            except Exception as e_r:
                error_msg = str(e_r)
                if "could not be loaded" in error_msg and "pd." in error_msg:
                    import re
                    match = re.search(r"(pd\.[a-zA-Z0-9\.]+)", error_msg)
                    if match:
                        missing_pkg = match.group(1)
                        print(f"    > Missing R package detected: {missing_pkg}")
                        installed = self._install_missing_r_package(missing_pkg)
                        if installed:
                            print("    > Retrying read_celfiles...")
                            raw_data = self.oligo.read_celfiles(cel_files_r)
                        else:
                            raise e_r
                    else:
                        raise e_r
                else:
                    raise e_r

            print("  - Running RMA normalization...")
            eset = self.oligo.rma(raw_data)
            exprs_matrix_r = self.biobase.exprs(eset)
            data_np = np.array(exprs_matrix_r)
            
            # Get Probe IDs
            try:
                probe_ids = list(exprs_matrix_r.rownames)
            except:
                probe_ids = list(self.biobase.featureNames(eset))
                
            df = pd.DataFrame(data_np)
            df.index = probe_ids  
            df.index.name = "ProbeID"
            
            # Clean Column Names
            clean_names = []
            for c in cel_files:
                filename = os.path.basename(c)
                filename = filename.replace('.gz', '').replace('.cel', '').replace('.CEL', '')
                clean_names.append(filename.split('_')[0])
            df.columns = clean_names

            # --- 4. Mapping to Locus ID (Modified Section) ---
            print("  - Mapping Probe IDs to Arabidopsis Locus IDs...")
            try:
                platform_id = gse.metadata.get('platform_id', [''])[0]
                if platform_id and platform_id in gse.gpls:
                    gpl = gse.gpls[platform_id]
                    annot = gpl.table
                    
                    # Priority List: Look for explicit Locus columns first, then Symbols
                    # 'ORF' and 'AGI' usually contain the AT code directly
                    candidates = ['ORF', 'AGI', 'Locus', 'Locus_id', 'Gene Symbol', 'GENE_SYMBOL', 'SYMBOL']
                    
                    target_col = None
                    for col in candidates:
                        if col in annot.columns:
                            target_col = col
                            break
                    
                    if target_col:
                        print(f"    > Found annotation column: '{target_col}'")
                        
                        # 1. Create a dictionary: ProbeID -> Raw Annotation (e.g., "FLC" or "AT5G10140")
                        # Handle multiple entries like "AT1G01010 /// AT1G01020" by taking the first one
                        annot['Cleaned_Val'] = annot[target_col].astype(str).apply(lambda x: x.split('///')[0].split('//')[0].strip())
                        
                        # Filter out empty strings or 'nan'
                        valid_annot = annot[annot['Cleaned_Val'].replace('nan', '') != '']
                        probe_to_raw_map = dict(zip(valid_annot['ID'], valid_annot['Cleaned_Val']))
                        
                        # 2. Get unique raw values to convert (optimization)
                        unique_raw_values = list(set(probe_to_raw_map.values()))
                        
                        # 3. Resolve these raw values to Locus IDs using the helper
                        raw_to_locus_map = self._map_symbols_to_locus_ids(unique_raw_values)
                        
                        # 4. Map the DataFrame
                        # Map Probe -> Raw Value -> Locus ID
                        df['Raw_Val'] = df.index.map(probe_to_raw_map)
                        df['LocusID'] = df['Raw_Val'].map(raw_to_locus_map)
                        
                        # 5. Drop rows where LocusID is NaN (Failed mapping or invalid ID)
                        initial_count = len(df)
                        df = df.dropna(subset=['LocusID'])
                        
                        # 6. Group by LocusID (average duplicates if multiple probes map to same gene)
                        df = df.drop(columns=['Raw_Val'])
                        df = df.groupby('LocusID').mean()
                        
                        print(f"    > Annotation complete. Reduced {initial_count} probes to {len(df)} validated Locus IDs.")
                    else:
                        raise ValueError(f"    > WARNING: Could not find suitable gene column in GPL {platform_id}.")
                else:
                    raise ValueError("    > Platform data not found in GSE object.")
                    
            except Exception as e_annot:
                print(f"    > Annotation step failed: {e_annot}. Aborting this study.")
                raise e_annot

            # --- 5. Save ---
            output_file = os.path.join(output_folder, f"{experiment_id}_RMA_LocusID.csv")
            df.to_csv(output_file)
            print(f"  - SUCCESS: Saved to {output_file}")
            return 0
            
        except Exception as e:
            print(f"  - RMA Failed for {experiment_id}: {e}")
            raise e
        
    def process_microarray_metadata_and_rma_old(self, experiment_id, raw_data_folder, output_folder, gse)-> int:
        print(f"\n[RMA Pipeline] Processing {experiment_id}...")
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        cel_files = []
        for root, dirs, files in os.walk(raw_data_folder):
            for f in files:
                if f.lower().endswith('.cel') or f.lower().endswith('.cel.gz'):
                    full_path = os.path.join(root, f)
                    cel_files.append(full_path)       
        if not cel_files:
            print(f"  - ABORTING: No CEL files found in {raw_data_folder}.")
            raise ValueError('No cel files found')

        print(f"  - Found {len(cel_files)} CEL files.")
        print("  - Processing metadata for all samples in study...")
        try:
            for gsm_name, gsm in gse.gsms.items():
                process_metadata(experiment_id, gse, gsm, save_path=output_folder)
        except Exception as e:
            print(f"    > Metadata processing failed: {e}")

        try:
            cel_files_r = ro.StrVector(cel_files)
            print("  - Reading CEL files (oligo)...")
            try:
                raw_data = self.oligo.read_celfiles(cel_files_r)
            except Exception as e_r:
                error_msg = str(e_r)
                if "could not be loaded" in error_msg and "pd." in error_msg:
                    import re
                    match = re.search(r"(pd\.[a-zA-Z0-9\.]+)", error_msg)
                    if match:
                        missing_pkg = match.group(1)
                        installed = self._install_missing_r_package(missing_pkg)
                        if installed:
                            print("    > Retrying read_celfiles...")
                            try:
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
            exprs_matrix_r = self.biobase.exprs(eset)
            data_np = np.array(exprs_matrix_r)
            try:
                probe_ids = list(exprs_matrix_r.rownames)
            except:
                probe_ids = list(self.biobase.featureNames(eset))
            df = pd.DataFrame(data_np)
            df.index = probe_ids  
            df.index.name = "ProbeID"
            clean_names = []
            for c in cel_files:
                filename = os.path.basename(c)
                if filename.lower().endswith('.gz'):
                    filename = filename[:-3]
                if filename.lower().endswith('.cel'):
                    filename = filename[:-4]
                clean_names.append(filename.split('_')[0])
            df.columns = clean_names
            print("  - Mapping Probe IDs to Genes...")
            try:
                platform_id = gse.metadata.get('platform_id', [''])[0]
                if platform_id and platform_id in gse.gpls:
                    gpl = gse.gpls[platform_id]
                    annot = gpl.table
                    candidates = ['Gene Symbol', 'GENE_SYMBOL', 'Gene_Symbol', 'SYMBOL', 'ORF', 'GeneSymbol']
                    gene_col = None
                    for col in candidates:
                        if col in annot.columns:
                            gene_col = col
                            break
                    if gene_col:
                        print(f"    > Found gene column: '{gene_col}'")
                        def clean_id(val):
                            try:
                                return str(int(float(val)))
                            except:
                                return str(val)
                        annot['ID'] = annot['ID'].apply(clean_id)
                        clean_genes = annot[gene_col].astype(str).apply(lambda x: x.split('///')[0].split('//')[0].strip())
                        probe_map = dict(zip(annot['ID'], clean_genes))
                        df.index = df.index.map(str)
                        df['GeneSymbol'] = df.index.map(probe_map)
                        df['GeneSymbol'] = df['GeneSymbol'].fillna(df.index.to_series())
                        df = df.groupby('GeneSymbol').mean()
                        print(f"    > Annotation complete. Collapsed {len(probe_map)} probes into {len(df)} genes.")
                    else:
                        raise ValueError(f"    > WARNING: Could not find Gene Symbol column in GPL ({platform_id}). Keeping Probe IDs.")
                else:
                    raise ValueError("    > Platform data not found in GSE object. Keeping Probe IDs.")
            except Exception as e_annot:
                print(f"    > Annotation step failed: {e_annot}. Keeping Probe IDs.")
                raise e_annot
            output_file = os.path.join(output_folder, f"{experiment_id}_RMA_Genes.csv")
            df.to_csv(output_file)
            print(f"  - SUCCESS: Saved to {output_file}")
            return 0
        except Exception as e:
            print(f"  - RMA Failed for {experiment_id}: {e}")
            raise e

def download_experiments_microarray(gse_list, download_dir, tracker:Microarray_tracker, download_raw=True, scan=True, output_folder='./microarray_processed_data/'):
    '''
    Downloads experiments using the tracker object to persist state (ignore/downloaded/processed).
    '''
    data_processor = Microarray_data_processing()
    if not os.path.exists(download_dir):
        os.makedirs(download_dir)

    valid_gse_ids = []
    
    # Path where the tracker will be auto-saved
    tracker_save_path = os.path.join(download_dir, "tracker_stats.json")

    # --- PROGRESS BAR ---
    for gse_id in tqdm(gse_list, desc='Processing Datasets', unit='study'):
        
        # 1. CHECK TRACKER: Skip if ignored or already processed (unless we want to re-verify)
        if tracker.is_ignored(gse_id):
            # tqdm.write(f"  [{gse_id}] Skipping (Ignored).")
            continue
            
        if tracker.is_processed(gse_id):
            # Already done, just add to valid list if needed
            valid_gse_ids.append(gse_id)
            continue

        # Define expected file paths
        soft_file_path = os.path.join(download_dir, f'{gse_id}_family.soft.gz')
        exp_folder = os.path.join(download_dir, gse_id)
        
        try:
            gse = None
            
            # --- PHASE 1: DOWNLOAD & CHECK METADATA ---
            if not tracker.is_downloaded(gse_id):
                try:
                    # Attempt download or load from disk if partial
                    gse = GEOparse.get_GEO(geo=gse_id, destdir=download_dir, silent=True)
                except Exception as e:
                    tqdm.write(f"  [{gse_id}] Metadata download failed: {e}")
                    tracker.mark_ignore(gse_id)
                    clean_files([soft_file_path])
                    tracker.save_to_json(tracker_save_path)
                    continue
                
                # Check Platform
                try:
                    platform = gse.metadata.get('platform_id', ['Unknown'])[0]
                    num_samples = len(gse.gsms)
                    
                    if platform != 'GPL198':
                        # Ignore this study forever
                        tracker.mark_ignore(gse_id)
                        clean_files([soft_file_path])
                        tracker.save_to_json(tracker_save_path)
                        continue
                except:
                    platform = 'Unknown'
                    num_samples = 0
                
                # Check for CEL files in Metadata
                has_cel = check_metadata_for_cel(gse) if download_raw else False

                # SCAN MODE
                if scan:
                    tracker.update_platform(platform, num_samples, used=has_cel)
                    if has_cel:
                        valid_gse_ids.append(gse_id)
                    
                    # In scan mode, we don't 'download' the raw files, so we clean up metadata
                    clean_files([soft_file_path])
                    # We do NOT mark as 'downloaded' or 'ignored' in scan mode usually,
                    # unless we want to remember we scanned it. 
                    # For now, let's assuming scanning doesn't change persistent state 
                    # OR you might want to mark it as ignored if no CEL files found.
                    continue

                # NORMAL MODE
                if download_raw:
                    if not has_cel:
                        tqdm.write(f'  [{gse_id}] SKIPPING: No .CEL files in metadata.')
                        tracker.update_platform(platform, num_samples, used=False)
                        tracker.mark_ignore(gse_id)
                        clean_files([soft_file_path])
                        tracker.save_to_json(tracker_save_path)
                        continue
                    
                    # Download Raw Files
                    if not os.path.exists(exp_folder):
                        os.makedirs(exp_folder)
                    
                    gse.download_supplementary_files(directory=exp_folder, download_sra=False)

                    if count_cel_files(exp_folder) > 0:
                        valid_gse_ids.append(gse_id)
                        tracker.update_platform(platform, num_samples, used=True)
                        
                        # SUCCESSFUL DOWNLOAD -> Update State
                        tracker.mark_downloaded(gse_id)
                        tracker.save_to_json(tracker_save_path)
                    else:
                        tqdm.write(f'  [{gse_id}] WARNING: Downloaded but 0 CEL files found.')
                        tracker.update_platform(platform, num_samples, used=False)
                        tracker.mark_ignore(gse_id)
                        clean_files([exp_folder, soft_file_path])
                        tracker.save_to_json(tracker_save_path)
                        continue
            
            # --- PHASE 2: PROCESSING ---
            # If we are here, it is downloaded (or we just downloaded it), and we are not in scan mode.
            if download_raw and not scan:
                # If we don't have the GSE object yet (e.g., resumed from 'downloaded' state)
                if gse is None:
                    # Try loading from local SOFT file first
                    try:
                        gse = GEOparse.get_GEO(filepath=soft_file_path, destdir=download_dir, silent=True)
                    except:
                        # Fallback to fetching ID if file missing (unlikely if marked downloaded)
                         gse = GEOparse.get_GEO(geo=gse_id, destdir=download_dir, silent=True)

                tqdm.write(f'  [{gse_id}] Processing RMA...')
                try:
                    data_processor.process_microarray_metadata_and_rma(gse_id, exp_folder, f'{output_folder}/{gse_id}', gse)
                    
                    # SUCCESSFUL PROCESSING -> Update State
                    tracker.mark_processed(gse_id)
                    tracker.save_to_json(tracker_save_path)
                    
                except Exception as e_proc:
                    tqdm.write(f'  [{gse_id}] PROCESSING FAILED: {e_proc}')
                    # If processing fails, we mark as ignored so we don't retry endlessly
                    tracker.mark_ignore(gse_id)
                    clean_files([soft_file_path, exp_folder])
                    tracker.save_to_json(tracker_save_path)
                    raise e_proc

        except Exception as e:
            tqdm.write(f'  [{gse_id}] FAILED: {e}')
            tracker.mark_ignore(gse_id)
            tracker.save_to_json(tracker_save_path)
            clean_files([soft_file_path, exp_folder])

    return valid_gse_ids
