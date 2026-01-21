import json
import os
from tqdm import tqdm
# from extraction_test_core import *
import sys
module_dir = './'
sys.path.append(module_dir)
from src.meta_data_processing.utils.extractors_full import *
from src.meta_data_processing.utils.llm_utils import get_condensed_labels,get_metadata_script
from src.meta_data_processing.utils.classes import LabelMap
from src.constants import *
from src.meta_data_processing.utils.llm_utils import get_condensed_labels
from src.meta_data_processing.utils.classes import LabelMap
import dotenv

dotenv.load_dotenv()

def load_labels_study(path):
    labels = {}
    for file in os.listdir(path):
        if file.endswith('.json'):
            file_path = os.path.join(path, file)
            study = file_path.split('/')[-1].split('.')[0]
            with open(file_path, 'r') as file:
                data = json.load(file)
            labels[study] = data
    return labels

def load_json(path:str):
    with open(path, 'r') as file:
        object = json.load(file)
    return object

import os
import json
import sys
from tqdm import tqdm

def condense_labels(Studies=Studies, in_folder='new_storage/processed_microarray_data/', saving_path=LABELS_PATH, llm_grounding:bool = True, model=None):
    os.makedirs(saving_path, exist_ok=True)
    labels = {}
    seen = LabelMap('./data/maps')
    
    # Check if input directory exists
    if not os.path.exists(in_folder):
        print(f"Error: Input folder '{in_folder}' not found.")
        return

    # 1. Iterate over Study Directories
    study_dirs = [d for d in os.listdir(in_folder) if os.path.isdir(os.path.join(in_folder, d))]
    
    for study_id in tqdm(study_dirs, desc="Processing Studies"):
        # Filter by Studies list if provided
        if Studies and (study_id not in Studies):
            continue

        study_path = os.path.join(in_folder, study_id)
        
        # Initialize study entry in labels if not present
        if study_id not in labels:
            labels[study_id] = {}

        # 2. Iterate over Sample JSON files inside the study folder
        sample_files = [f for f in os.listdir(study_path) if f.endswith('.json')]
        
        for sample_file in sample_files:
            file_path = os.path.join(study_path, sample_file)
            
            try:
                # Load the JSON containing both study and sample metadata
                with open(file_path, 'r') as f:
                    data = json.load(f)
                
                # Extract relevant sections
                study_meta = data.get('study_metadata', {})
                sample_meta = data.get('sample_metadata', {})
                sample_id = data.get('sample_id', sample_file.replace('.json', ''))
                
            except Exception as e:
                print(f"Error loading {sample_file}: {e}")
                continue

            # 3. Dynamic Extractor Check & Generation
            extractor_name = f'{study_id}_extractor'
            
            # Check if the extractor function exists in the global scope
            if extractor_name not in globals():
                tqdm.write(f"Warning: Extractor function '{extractor_name}' not defined.")
                tqdm.write(f"Generating new extractor for {study_id}...")
                
                try:
                    # Generate the script using the LLM utility
                    # We pass the sample metadata as a reference for the schema generation
                    new_script = get_metadata_script(sample_meta, study_id=study_id, model=model)
                    
                    # Append the new function to the extractors file
                    with open('extractors_full.py', 'a') as f:
                        f.write('\n\n')
                        f.write(new_script)
                    
                    # Execute the new script string to load the function into the current runtime
                    exec(new_script, globals())
                    tqdm.write(f"  -> Successfully generated and loaded {extractor_name}.")
                    
                except Exception as e:
                    tqdm.write(f"  -> Failed to generate extractor for {study_id}: {e}")
                    continue # Skip this sample if we can't extract

            # 4. Run the Extractor
            try:
                # Retrieve the function from globals and execute
                extractor_func = globals()[extractor_name]
                python_object = extractor_func(sample_meta)
                
                # --- Post-Processing Logic (Same as original) ---
                
                # Handle missing treatments
                if not python_object.get('treatment'):
                    python_object['treatment'] = ['no treatment/control']
                
                # Debug print for specific condition
                if "Light 24h" in python_object['treatment']:
                    print(f'study: {study_id} sample: {sample_id}')
                
                # LLM Grounding / Mapping
                if llm_grounding and seen.check_past(python_object):
                    # Note: Using study_meta extracted from the JSON
                    condensed = dict(get_condensed_labels(study_info=study_meta, sample_info=python_object))
                    seen.add_mapping(python_object, condensed)
                    python_object = condensed
                else:
                    condensed = seen.apply_mappings(python_object)
                
                # Helper to flatten lists
                def flatten(lst, ret=None):
                    if ret is None: ret = []
                    for el in lst:
                        if isinstance(el, list):
                            flatten(el, ret)
                        else:
                            ret.append(el)
                    return ret

                # Flatten treatments and save
                condensed['treatment'] = list(set(flatten(condensed['treatment'])))
                labels[study_id][sample_id] = dict(condensed)
                seen.save_map()

            except Exception as e:
                print(f"Error processing sample {sample_id} in {study_id}: {e}")
                continue

    # 5. Save Results
    print("Saving condensed labels...")
    for study in labels:
        if labels[study]: # Only save if we actually processed data for this study
            with open(f'{saving_path}/{study}.json', 'w') as handle:
                json.dump(labels[study], handle, indent=4)

if __name__ == '__main__':
    condense_labels(model='gemini-1.5-pro')

    labels_1 = load_labels_study(LABELS_PATH)
    res = []
    for st in labels_1:
        for sam in labels_1[st]:
            save = labels_1[st][sam]
            save['id'] = sam
            res.append(labels_1[st][sam])
            
    with open(f'llm_condensed_labels.json', 'w') as handle:
        json.dump(res, handle)
