import json
import os
from tqdm import tqdm
import sys
module_dir = './'
sys.path.append(module_dir)
from src.constants import *


def batch_metadata_by_study():
    lables = {}
    study_metdata = dict.fromkeys(Studies,None)
    for file in tqdm(os.listdir(METADATA_OUTPUT_DIR)):
        if file.endswith('.json'):
            file_path = os.path.join(METADATA_OUTPUT_DIR, file)
            study_id,sample_id = file_path.split('/')[-1].split('.')[0].split('_')
            if not (study_id in Studies):
                continue
            if  study_id in lables:
                if sample_id in lables[study_id]:
                    continue
            with open(file_path, 'r') as file:
                data = json.load(file)
            if  not(study_id in lables):
                lables[study_id] = {}
            meta_data = data.copy()
            meta_data['sample_title'] = data['sample_metadata']['title']
            meta_data['sample_source_name'] = data['sample_metadata']['source_name_ch1']
            try:
                meta_data['sample_characteristicts'] = data['sample_metadata']['characteristics_ch1']
            except KeyError:
                meta_data['sample_characteristicts'] = ''
            try:
                meta_data['sample_extraction_protocol'] = data['sample_metadata']['extract_protocol_ch1']
            except KeyError:
                meta_data['sample_extraction_protocol'] = ''
            del meta_data['sample_metadata']
            del meta_data['study_metadata']
            if study_metdata[study_id] is None:
                study_metdata[study_id] = {'study_metadata': data['study_metadata']}
            study_metdata[study_id][sample_id] = meta_data # pyright: ignore[reportOptionalSubscript]
    os.makedirs('{METADATA_OUTPUT_DIR}/study_batch_metadata/',exist_ok=True)
    for study in study_metdata:
        with open(f'{METADATA_OUTPUT_DIR}/study_batch_metadata/{study}.json', 'w') as handle:
            json.dump(study_metdata[study], handle)

if __name__ == '__main__':
    batch_metadata_by_study()