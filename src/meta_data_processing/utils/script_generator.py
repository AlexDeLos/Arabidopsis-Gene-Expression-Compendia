import json
import os
from tqdm import tqdm
from llm_utils import get_metadata, get_metadata_script
def load_labels_study(path):
    labels = {}
    for file in tqdm(os.listdir(path)):
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


def condense_labels_script(Studies,in_folder,model, saving_path):
    os.makedirs(saving_path,exist_ok=True)
    try:
        labeling_script = load_labels_study(saving_path)
    except:
        labeling_script= {}

    count = 0
    for file in tqdm(os.listdir(in_folder)):
        if file.endswith('.json'):
            file_path = os.path.join(in_folder, file)
            
            
            study_id = file_path.split('/')[-1].split('.')[0].split('_')[0]
            # if study_id != 'GSE41935' or _ !='GSM1027685':
            #     continue

            if not (study_id in Studies):
                continue
            if study_id in labeling_script:
                continue

            dict = load_json(file_path)
            study_info = dict['Study metadata']
            del dict['Study metadata']
            samples_info = dict
            if not (study_id in labeling_script):
                python_code = get_metadata_script(study_info, samples_info,study_id, model = model, temp=0)
                python_code = python_code.replace('```python', '')
                python_code = python_code.replace('```', '')
                print(python_code, file=open('extractors_full.py', 'a'))
                labeling_script[study_id] = str(python_code)
            else:
                pass
            
            
            count = count +1
            if count%10 == 0:
                for study in labeling_script:
                    with open(f'{saving_path}/{study}.json', 'w') as handle:
                        json.dump(labeling_script[study], handle)
                print("saved")

    for study in labeling_script:
        with open(f'{saving_path}/{study}.json', 'w') as handle:
            json.dump(labeling_script[study], handle)

#! TEMP TODO: delete
# Studies = ['GSE201609','GSE46205','GSE27548','GSE41935','GSE34188','GSE27550','GSE76827','GSE5620','GSE110079','GSE60960','GSE5624','GSE5622','GSE126373','GSE5628','GSE63128']
Studies = ['GSE44053', 'GSE77815', 'GSE16474', 'GSE4062', 'GSE9415', 'GSE18624', 'GSE40061', 'GSE112161', 'GSE20494', 'GSE79997', 'GSE27552', 'GSE16222', 'GSE110857', 'GSE4760', 'GSE46205', 'GSE71001', 'GSE58616', 'GSE162310', 'GSE22107', 'GSE62163', 'GSE51897', 'GSE72949', 'GSE90562', 'GSE5628', 'GSE26266', 'GSE34188', 'GSE34595', 'GSE76827', 'GSE119383', 'GSE65046', 'GSE11758', 'GSE65414', 'GSE37408', 'GSE5624', 'GSE10643', 'GSE15577', 'GSE11538', 'GSE70861', 'GSE6583', 'GSE27551', 'GSE12619', 'GSE121225', 'GSE108070', 'GSE78713', 'GSE110079', 'GSE63128', 'GSE60960', 'GSE37118', 'GSE79681', 'GSE63372', 'GSE5622', 'GSE26983', 'GSE27550', 'GSE19603', 'GSE95202', 'GSE53308', 'GSE16765', 'GSE71855', 'GSE58620', 'GSE24177', 'GSE35258', 'GSE10670', 'GSE49418', 'GSE18666', 'GSE83136', 'GSE44655', 'GSE27549', 'GSE19700', 'GSE103398', 'GSE63522', 'GSE201609', 'GSE5620', 'GSE66369', 'GSE2268', 'GSE71237', 'GSE48474', 'GSE41935', 'GSE27548', 'GSE5623', 'GSE72050', 'GSE126373']

in_folder = 'study_batch_metadata'
experiment = '1.1'
model = 'gemini-2.5-flash'
saving_path = f'scripts/{model}/{experiment}'
condense_labels_script(Studies,in_folder,model, saving_path)