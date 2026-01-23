import json
import os
from tqdm import tqdm
import dotenv
import sys
module_dir = './'
sys.path.append(module_dir)
import src.meta_data_processing.utils.extractors_full as extractors_full
from src.constants import *
from src.meta_data_processing.utils.classes import *
from src.meta_data_processing.utils.llm_utils import get_batch_labels_treatment,get_batch_labels_tissues,get_metadata_script
from src.meta_data_processing.utils.universal_extractor import UniversalExtractor
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

def save_labels(labels,saving_path):
    # 5. Save Results
    print("Saving condensed labels...")
    for study in labels:
        if labels[study]: # Only save if we actually processed data for this study
            with open(f'{saving_path}/{study}.json', 'w') as handle:
                json.dump(labels[study], handle, indent=4)


def condense_labels(in_folder, saving_path, Studies=None,extractors_path = 'src/meta_data_processing/utils/extractors_full.py'):
    os.makedirs(saving_path, exist_ok=True)
    uni_extractor = UniversalExtractor()
    # --- COMPONENT 1: Restore LabelMap ---
    # This loads your historical mappings (map.json, map_treatment.json, etc.)
    seen = LabelMap('./new_storage/maps')
    
    # Initialize Optimizer
    optimizer = GroundingOptimizer()
    
    study_dirs = [d for d in os.listdir(in_folder) if os.path.isdir(os.path.join(in_folder, d))]
    
    for study_id in tqdm(study_dirs, desc="Processing Studies"):
        # 1. Skip if filtered
        if Studies and study_id not in Studies:
            continue
            
        # 2. Check if already done
        output_file = os.path.join(saving_path, f"{study_id}.json")
        if os.path.exists(output_file):
            continue

        study_path = os.path.join(in_folder, study_id)
        raw_samples = []

        # --- PHASE 1: EXTRACTION ---
        sample_files = [f for f in os.listdir(study_path) if f.endswith('.json')]
        for sample_file in sample_files:
            try:
                with open(os.path.join(study_path, sample_file), 'r') as f:
                    data = json.load(f)
                
                # Verify extractor exists before processing files
                if False:
                    extractor_name = f'{study_id}_extractor'
                    if hasattr(extractors_full, extractor_name):
                        extractor = getattr(extractors_full, extractor_name)
                    elif extractor_name in globals():
                        extractor = globals()[extractor_name]
                    else:
                        # You might want to generate it here if missing, or skip
                        new_extractor = get_metadata_script(sample_metadata=data,study_id=study_id)
                        
                        # Append to file
                        with open(extractors_path, 'a') as f:
                            f.write('\n\n')
                            f.write(new_extractor)
                        
                        # Load into current runtime
                        exec(new_extractor, globals())
                        extractor = globals()[extractor_name] # Now we can grab it
                    
                    
                    extracted = extractor(data.get('sample_metadata', {}))
                else:
                    extracted = uni_extractor.extract(
                        sample_metadata=data.get('sample_metadata', {}),
                        study_metadata=data.get('study_metadata', {})
                       )
                
                # IMPORTANT: Keep 'sample_id' for tracking, but we will remove it before saving
                extracted['sample_id'] = data.get('sample_id', sample_file.replace('.json',''))
                raw_samples.append(extracted)
            except Exception as e:
                print(f"Extraction error {study_id}: {e}")

        if not raw_samples:
            continue

        # --- PHASE 2: BATCH GROUNDING ---
        try:
            # CHANGE: Pass the whole 'seen' object as 'label_map_obj'
            grounded_samples = optimizer.batch_process_study(
                data=data, # type: ignore
                extracted_samples=raw_samples,
                llm_func=get_batch_labels_treatment,
                llm_func_tis= get_batch_labels_tissues,
                label_map_obj=seen
            )
            
            # --- COMPONENT 2: Sync and Save Progress ---
            final_output = {}
            
            for raw, grounded in zip(raw_samples, grounded_samples):
                
                # Update the LabelMap (seen.add_mapping handles the complexity of splitting dicts)
                seen.add_mapping(raw, grounded)
                
                s_id = grounded.pop('sample_id')
                final_output[s_id] = grounded

            seen.save_map()
            
            # 4. Save the Study Data
            with open(output_file, 'w') as f:
                json.dump(final_output, f, indent=4)
                
        except Exception as e:
            print(f"Grounding error {study_id}: {e}")
            # Optional: Save partial progress or handle error gracefully



if __name__ == '__main__':
    condense_labels(in_folder='new_storage/processed_microarray_data/',saving_path=LABELS_PATH)#, Studies=['GSE5622','GSE119383','GSE10670'])


    labels_1 = load_labels_study(LABELS_PATH)
    res = []
    for st in labels_1:
        for sam in labels_1[st]:
            save = labels_1[st][sam]
            save['id'] = sam
            res.append(labels_1[st][sam])
            
    with open(f'llm_condensed_labels.json', 'w') as handle:
        json.dump(res, handle)
