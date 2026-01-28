import os
from Bio import Entrez
import sys
import shutil

# Ensure we can import the local modules
module_dir = './'
sys.path.append(module_dir)

# Import the processing classes
def search_geo_accessions(query, max_results=20, filter_organism="Arabidopsis thaliana"):
    '''
    Searches NCBI GEO for datasets matching the query and strict-filters the organism.
    '''
    print(f'Searching GEO with query: {query}...')
    try:
        handle = Entrez.esearch(db='gds', term=query, retmax=max_results)
        record = Entrez.read(handle)
        handle.close()
        
        id_list = record['IdList']
        print(f'Found {len(id_list)} datasets (limit set to {max_results}).')
        
        if not id_list:
            return []
            
        gse_ids = []
        batch_size = 100 
        
        for i in range(0, len(id_list), batch_size):
            batch = id_list[i:i+batch_size]
            handle = Entrez.esummary(db='gds', id=','.join(batch))
            summaries = Entrez.read(handle)
            handle.close()
            
            for summary in summaries:
                acc = summary.get('Accession', '')
                
                # --- FIX: Verify Organism ---
                # 'taxon' field contains the organism name (e.g., "Homo sapiens")
                taxon = summary.get('taxon', '')
                
                # If we requested a filter, checking it matches the metadata
                if filter_organism and filter_organism.lower() not in taxon.lower():
                    # print(f"  Skipping {acc} (Organism: {taxon})") # Optional debug
                    continue

                if acc.startswith('GSE'):
                    gse_ids.append(acc)
        
        print(f"Returned {len(gse_ids)} validated {filter_organism} IDs.")
        return list(set(gse_ids))
        
    except Exception as e:
        print(f'Error during search: {e}')
        return []

def count_cel_files(directory):
    '''
    Recursively counts unique CEL files in the directory.
    '''
    cel_count = 0
    for root, dirs, files in os.walk(directory):
        for f in files:
            if f.lower().endswith('.cel') or f.lower().endswith('.cel.gz'):
                cel_count += 1
    return cel_count

def check_metadata_for_cel(gse):
    '''
    Inspects the GSM metadata to see if CEL files are even listed.
    '''
    for gsm_name, gsm in gse.gsms.items():
        for key, value in gsm.metadata.items():
            val_str = str(value).lower()
            if 'supplementary_file' in key and '.cel' in val_str:
                return True
    return False

def clean_files(files:list[str]):
    for file in files:
        if os.path.exists(file):
            if os.path.isdir(file):
                shutil.rmtree(file)
            elif os.path.isfile(file):
                os.remove(file)

def check_metadata_for_counts(gse):
    """
    Checks if the study has supplementary files that look like count matrices.
    """
    # ADDED '.tar' to valid extensions
    valid_exts = ('.txt', '.csv', '.tsv', '.tab', '.xlsx', 'counts', 'abundance', 'fpkm', 'tpm', '.tar')
    
    # Check GSE level supplements
    for name, metadata in gse.metadata.items():
        if 'supplementary_file' in name:
            for url in metadata:
                url_lower = str(url).lower()
                if any(ext in url_lower for ext in valid_exts):
                    return True
                        
    # Check GSM level
    for gsm_name, gsm in gse.gsms.items():
        for key, value in gsm.metadata.items():
            val_str = str(value).lower()
            if 'supplementary_file' in key and any(ext in val_str for ext in valid_exts):
                return True
    return False

def check_metadata_for_sra_boolean(gse):
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
def check_metadata_for_sra(gse, output_dir):
    """
    Checks if a given GSE ID contains SRA (RNA-Seq) data.
    Returns metadata dict if valid, False otherwise.
    """
    try:
        # 1. Download/Parse Metadata
        # We silence the output to keep the console clean
        
        valid_sra_samples = 0
        detected_platform = "Unknown"
        
        # 2. Iterate through Samples (GSMs)
        for gsm_name, gsm in gse.gsms.items():
            is_sra = False
            
            # Check A: Relation field (Standard for SRA)
            # Looks for "SRA: SRX..." or "BioProject: PRJNA..."
            for relation in gsm.metadata.get('relation', []):
                if "SRA:" in relation or "BioProject:" in relation:
                    is_sra = True
                    break
            
            # Check B: Explicit SRA field
            if not is_sra and 'sra_id' in gsm.metadata:
                is_sra = True

            # Check C: Supplementary Check (Fallback)
            if not is_sra:
                for key, value in gsm.metadata.items():
                    val_str = str(value).lower()
                    if "supplementary_file" in key and (".fastq" in val_str or ".fq" in val_str):
                        is_sra = True
                        break
            
            if is_sra:
                valid_sra_samples += 1
                # Grab the platform (GPL) from the first valid sample
                if detected_platform == "Unknown":
                    # 'platform_id' is usually a list ["GPL1234"]
                    detected_platform = gsm.metadata.get('platform_id', ['Unknown'])[0]

        # 3. Return Logic
        if valid_sra_samples > 0:
            return {
                'platform': detected_platform,
                'n_samples': valid_sra_samples
            }
        else:
            return False

    except Exception as e:
        # If GEOparse fails (network issue, deleted record), return False to skip
        print(f"Metadata check failed for {gse}: {e}")
        return False