import os
import GEOparse
from Bio import Entrez
import sys
import shutil
from tqdm import tqdm
import json
# Ensure we can import the local modules
module_dir = './'
sys.path.append(module_dir)

# Import the processing function from your other file
from src.data_importing.microarray_data_processing import Data_processing

def search_geo_accessions(query, max_results=20):
    """
    Searches NCBI GEO for datasets matching the query.
    """
    print(f"Searching GEO with query: {query}...")
    try:
        handle = Entrez.esearch(db="gds", term=query, retmax=max_results)
        record = Entrez.read(handle)
        handle.close()
        
        id_list = record["IdList"]
        print(f"Found {len(id_list)} datasets (limit set to {max_results}).")
        
        if not id_list:
            return []
            
        gse_ids = []
        batch_size = 100 
        
        # We can also add a progress bar here for the search batches if you like
        for i in range(0, len(id_list), batch_size):
            batch = id_list[i:i+batch_size]
            handle = Entrez.esummary(db="gds", id=",".join(batch))
            summaries = Entrez.read(handle)
            handle.close()
            
            for summary in summaries:
                acc = summary.get("Accession", "")
                if acc.startswith("GSE"):
                    gse_ids.append(acc)
        
        return list(set(gse_ids))
        
    except Exception as e:
        print(f"Error during search: {e}")
        return []

def count_cel_files(directory):
    """
    Recursively counts unique CEL files in the directory.
    """
    cel_count = 0
    for root, dirs, files in os.walk(directory):
        for f in files:
            if f.lower().endswith('.cel') or f.lower().endswith('.cel.gz'):
                cel_count += 1
    return cel_count

def check_metadata_for_cel(gse):
    """
    Inspects the GSM metadata to see if CEL files are even listed.
    """
    # print("  - Inspecting metadata for CEL files...") # Commented out to keep progress bar clean
    for gsm_name, gsm in gse.gsms.items():
        for key, value in gsm.metadata.items():
            val_str = str(value).lower()
            if "supplementary_file" in key and ".cel" in val_str:
                return True
    return False

def clean_files(files:list[str]):
    for file in files:
        if os.path.exists(file):
            shutil.rmtree(file)


def download_experiments_microarray(gse_list, download_dir, tracker, download_raw=True, scan=True,output_folder='./microarray_processed_data/'):
    """
    Downloads experiments. 
    If a study is skipped or fails (in scan=False mode), its files are deleted to save space.
    """
    data_processor = Data_processing()
    if not os.path.exists(download_dir):
        os.makedirs(download_dir)

    valid_gse_ids:list = []
    seen_path:str = f'{download_dir}seen.json'
    try:
        with open(seen_path, 'r') as file:
            seen_studies = json.load(file)
        #turn lists to sets
        for key in seen_studies:
            seen_studies[key] = set(seen_studies[key])
            
    except:
        seen_studies:dict = {
            "ignore":set(),
            "downloaded": set(),
            "processed":set()
        }
    # --- PROGRESS BAR ---
    for gse_id in tqdm(gse_list, desc="Processing Datasets", unit="study"):
        # Define expected file paths early for cleanup logic
        soft_file_path = os.path.join(download_dir, f"{gse_id}_family.soft.gz")
        exp_folder = os.path.join(download_dir, gse_id)
        downloaded = os.path.isdir(exp_folder)
        processed = os.path.isfile(f'{output_folder}{gse_id}/{gse_id}_RMA_Genes.csv')
        if processed:
            seen_studies['processed'].add(gse_id)
        if downloaded:
            seen_studies['downloaded'].add(gse_id)
        try:
            if (not downloaded) and  not(gse_id in seen_studies["ignore"]):
                # 1. Download Metadata (Strictly necessary for stats)
                try:
                    gse = GEOparse.get_GEO(destdir=download_dir, silent=True,filepath=soft_file_path)
                except:
                    gse = GEOparse.get_GEO(geo=gse_id, destdir=download_dir, silent=True)
                
                # --- EXTRACT STATS ---
                try:
                    platform = gse.metadata.get('platform_id', ['Unknown'])[0]
                    num_samples = len(gse.gsms)
                    if platform != 'GPL198':
                        seen_studies['ignore'].add(gse_id)
                        clean_files([soft_file_path])
                        continue
                except Exception as e:
                    platform = 'Unknown'
                    num_samples = 0
                
                # Check Metadata for CEL availability
                has_cel_reference = False
                if download_raw:
                    has_cel_reference = check_metadata_for_cel(gse)

                # --- SCAN MODE ---
                if scan:
                    tracker.update_platform(platform, num_samples, used=has_cel_reference)
                    if has_cel_reference:
                        valid_gse_ids.append(gse_id)
                    # CLEANUP: Always remove metadata in scan mode
                    clean_files([soft_file_path])
                    continue

                # --- NORMAL PROCESSING MODE (scan=False) ---
                
                # 2. Pre-Check (Metadata)
                if download_raw:
                    if not has_cel_reference:
                        tqdm.write(f"  [{gse_id}] SKIPPING: No .CEL files in metadata.")
                        tracker.update_platform(platform, num_samples, used=False)
                        # CLEANUP: Remove the SOFT file since we won't use this study
                        clean_files([soft_file_path])
                        seen_studies["ignore"].add(gse_id)
                        continue

                # 3. Download Supplementary Files
                found_raw = 0
                
                if download_raw:
                    if not os.path.exists(exp_folder):
                        os.makedirs(exp_folder)
                    
                    gse.download_supplementary_files(directory=exp_folder, download_sra=False)

                    # 4. Verify disk content
                    found_raw = count_cel_files(exp_folder)
                    
                    if found_raw > 0:
                        valid_gse_ids.append(gse_id)
                        tracker.update_platform(platform, num_samples, used=True)
                        seen_studies['downloaded'].add(gse_id)
                    else:
                        tqdm.write(f"  [{gse_id}] WARNING: Downloaded but 0 CEL files found.")
                        tracker.update_platform(platform, num_samples, used=False)
                        
                        # CLEANUP: Remove the useless folder AND the SOFT file
                        clean_files([exp_folder,soft_file_path])
                        seen_studies['ignore'].add(gse_id)
                        continue
            else:
                #what fo do if we did not download?
                pass
            # 5. Process
            if download_raw and not (gse_id in seen_studies['ignore']):
                tqdm.write(f"  [{gse_id}] Processing RMA...")
                if not processed:
                    try:
                        data_processor.process_microarray_metadata_and_rma(gse_id, exp_folder, f"{output_folder}/{gse_id}", gse)
                        seen_studies["processed"].add(gse_id)
                    except Exception as e:
                        seen_studies["ignore"].add(gse_id)
                        clean_files([soft_file_path,exp_folder])
                        raise e

            # if we already have it downloaded we should check that it is in the seen
            # seen_studies["seen"].add(gse_id)
            # if processed:
            #     seen_studies["used"].add(gse_id)
        except Exception as e:
            tqdm.write(f"  [{gse_id}] FAILED: {e}")
            # CLEANUP ON FAILURE
            # If it crashed, we don't want partial files taking up space
            clean_files([soft_file_path,exp_folder])
    #turn sets to lits
    for key in seen_studies:
        seen_studies[key] = list(seen_studies[key])
    json_str = json.dumps(seen_studies, indent=4)
    with open(seen_path, "w") as f:
        f.write(json_str)
    return valid_gse_ids


def check_metadata_for_sra(gse):
    """
    Checks if the study has SRA (Sequence Read Archive) links.
    This indicates 'Raw Data' (FASTQ) is available in SRA, even if not in GEO.
    """
    for gsm_name, gsm in gse.gsms.items():
        # Method 1: Check Relation field
        for relation in gsm.metadata.get('relation', []):
            if "SRA:" in relation or "BioProject:" in relation:
                return True
        # Method 2: Check explicit SRA field
        if 'sra_id' in gsm.metadata:
            return True
            
        # Method 3: Check supplementary files for FASTQ (Rare in GEO, but possible)
        for key, value in gsm.metadata.items():
            val_str = str(value).lower()
            if "supplementary_file" in key and (".fastq" in val_str or ".fq" in val_str):
                return True
    return False

def download_experiments_RNA_seq(gse_list, output_dir, tracker, download_raw=True, scan=True):
    """
    Specific downloader/scanner for RNA-Seq data.
    - Checks for SRA availability (Raw Data).
    - Can download supplementary counts if available.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    valid_gse_ids = []

    for gse_id in tqdm(gse_list, desc="Processing RNA-Seq", unit="study"):
        try:
            # 1. Download Metadata
            gse = GEOparse.get_GEO(geo=gse_id, destdir=output_dir, silent=True)
            
            # --- EXTRACT STATS ---
            try:
                platform = gse.metadata.get('platform_id', ['Unknown'])[0]
                num_samples = len(gse.gsms)
            except:
                platform = 'Unknown'
                num_samples = 0
            
            # 2. Check for RAW Data (SRA Availability)
            has_raw_sra = check_metadata_for_sra(gse)

            # --- SCAN MODE ---
            if scan:
                # Update tracker with findings
                tracker.update_platform(platform, num_samples, has_raw=has_raw_sra)
                
                if has_raw_sra:
                    valid_gse_ids.append(gse_id)

                # CLEANUP: Delete the soft file to save space
                soft_file_path = os.path.join(output_dir, f"{gse_id}_family.soft.gz")
                if os.path.exists(soft_file_path):
                    os.remove(soft_file_path)
                
                continue

            # --- NORMAL MODE (Scan=False) ---
            # If we are here, we might want to download supplementary files (Count matrices)
            if download_raw:
                if not has_raw_sra:
                    tqdm.write(f"  [{gse_id}] SKIPPING: No SRA link found.")
                    tracker.update_platform(platform, num_samples, has_raw=False)
                    continue
                
                # If it has raw data, we update tracker as successful
                tracker.update_platform(platform, num_samples, has_raw=True)
                valid_gse_ids.append(gse_id)

                # NOTE: For RNA-Seq, we usually do NOT download the SRA files via this script
                # because they are massive (TB). We usually just download the metadata/supp files here.
                # If you really want to download supplementary files (like counts.txt):
                
                exp_folder = os.path.join(output_dir, gse_id)
                if not os.path.exists(exp_folder):
                    os.makedirs(exp_folder)
                
                # tqdm.write(f"  [{gse_id}] Downloading supplementary files (counts/tables)...")
                gse.download_supplementary_files(directory=exp_folder, download_sra=False)

        except Exception as e:
            tqdm.write(f"  [{gse_id}] FAILED: {e}")

    return valid_gse_ids
