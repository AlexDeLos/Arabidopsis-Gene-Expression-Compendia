import GEOparse
import pandas as pd
import logging
import sys, os
module_dir = './'
sys.path.append(module_dir)
from src.constants import *
# It's assumed your 'helpers.py' file with 'get_geo_list' exists in the same directory.
from data_importing.helpers.helpers import *

# --- Configuration ---
# Set up basic logging to track progress and errors
logging.basicConfig(level=logging.WARNING, format='%(asctime)s - %(levelname)s - %(message)s')

# Define paths for input and output
# GEO_DOWNLOAD_DIR = './downloads/geo_downloads/'
# METADATA_OUTPUT_DIR = './downloads/metadata/'
# COMBINED_DATA_OUTPUT_FILE = './downloads/combined_expression_data.csv'
# SOFT_PATH = './downloads/old_geo_downloads/'
def import_data(track:tracker):
    """
    Main function to download, process, and combine GEO data into a single DataFrame.
    """
    setup_directories()
    
    # Load the list of GEO studies and the master gene list
    try:
        geo_list = Studies#get_geo_list('core_lists/data_addresses.csv')
        df_index = pd.read_csv('core_lists/genes_list.csv', index_col=0)
        master_gene_list_df = pd.read_csv('core_lists/genes_list.csv', index_col=0)
    except FileNotFoundError as e:
        logging.error(f"Core list file not found: {e}. Please ensure 'core_lists' directory is correct.")
        return

    # Initialize the final DataFrame with the index from the master gene list
    final_df = pd.DataFrame(index=master_gene_list_df.index)
    
    for geo_accession in geo_list:
        logging.info(f"--- Processing study: {geo_accession} ---")
        try:
            gse = GEOparse.get_GEO(geo=geo_accession, destdir=GEO_DOWNLOAD_DIR, silent=True)
        except Exception as e:
            logging.warning(f"Could not download or parse {geo_accession}. Error: {e}")
            try:
                gse = GEOparse.get_GEO(filepath=f"{SOFT_PATH}{geo_accession}_family.soft.gz")
            except Exception as ex:
                logging.error(f"Error whilt getting the local version of {geo_accession}. Error: {ex}")
                continue

        gpl_names:list = list(gse.gpls.keys())# this run this for both gpls[0]
        for gpl_name in gpl_names:
            try:
                gpl_table = gse.gpls[gpl_name].table
                # probe_to_gene_map = create_probe_to_gene_map(gpl)
                probe_to_gene_map = dict(zip(gpl_table['ID'], gpl_table['ORF'].map(mapping)))
                if probe_to_gene_map is None:
                    logging.error(f"Skipping study {geo_accession} due to missing gene map, it is RNA seq.")
                    continue
            except Exception as e:
                logging.error(f"Error processing platform for {geo_accession}: {e}")
                continue

            for gsm_name, gsm in gse.gsms.items():
                logging.info(f"Processing sample: {gsm_name}")
                try:
                    # Goal 1: Store filtered metadata (unchanged)
                    process_metadata(geo_accession, gse, gsm)
                    
                    # Goal 2: Process data and get a DataFrame to append
                    #TODO: check if GSM680342 is in the final df
                    processed_sample_df = process_sample_data(geo_accession, gsm, probe_to_gene_map)
                    
                    # If data was processed successfully, add it to the final DataFrame
                    if processed_sample_df is not None and not processed_sample_df.empty:
                        track.update_counter(gsm_name,gse.gpls[gpl_name],geo_accession,list(gsm.table['VALUE']))
                        processed_sample_df.index = processed_sample_df.index.str.upper()
                        processed_sample_df = handle_duplicates(df_index,processed_sample_df)
                        complete_in = pd.concat([df_index, processed_sample_df], axis=1)
                        # complete_in = complete_in.transform(lambda x: x.fillna(x.mean()))
                        final_df = pd.concat([final_df, complete_in], axis=1)
                        # processed_sample_df = processed_sample_df.loc[processed_sample_df.index.isin(final_df.index)]
                        # final_df = pd.concat([final_df, processed_sample_df], axis=1)
                        # final_df = final_df.join(processed_sample_df)
                        logging.info(f"Appended data for {gsm_name} to the main DataFrame.")

                except Exception as e:
                    logging.error(f"An unexpected error occurred while processing sample {gsm_name}: {e}")

    logging.info("--- All studies processed. Finalizing the combined DataFrame. ---")
    
    # After merging, some columns might be all NaN if they had no overlapping genes.
    # It's good practice to drop them.
    # final_df.dropna(axis=1, how='all', inplace=True)

    # Save the final combined DataFrame to a single CSV file
    final_df.to_csv(COMBINED_DATA_OUTPUT_FILE)
    track.store()
    logging.info(f"✅ Success! Combined DataFrame saved to '{COMBINED_DATA_OUTPUT_FILE}'")

if __name__ == "__main__":
    track = tracker(STORAGE_DIR)
    import_data(track)
    track.load()
    track.log_transform_high_intensity_studies()
    track.plot()