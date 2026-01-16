import os
from Bio import Entrez
import sys
import shutil

# Ensure we can import the local modules
module_dir = './'
sys.path.append(module_dir)

# Import the processing classes
def search_geo_accessions(query, max_results=20):
    '''
    Searches NCBI GEO for datasets matching the query.
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
                if acc.startswith('GSE'):
                    gse_ids.append(acc)
        
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