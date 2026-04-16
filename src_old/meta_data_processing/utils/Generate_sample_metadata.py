import json
import os
from tqdm import tqdm

def generate_metadata(folder_path = '/home/alex/Documents/GitHub/Data_collection/metadata/output_3',Studies = ['GSE44053', 'GSE77815', 'GSE16474', 'GSE4062', 'GSE9415', 'GSE18624', 'GSE40061', 'GSE112161', 'GSE20494', 'GSE79997', 'GSE27552', 'GSE16222', 'GSE110857', 'GSE4760', 'GSE46205', 'GSE71001', 'GSE58616', 'GSE162310', 'GSE22107', 'GSE62163', 'GSE51897', 'GSE72949', 'GSE90562', 'GSE5628', 'GSE26266', 'GSE34188', 'GSE34595', 'GSE76827', 'GSE119383', 'GSE65046', 'GSE11758', 'GSE65414', 'GSE37408', 'GSE5624', 'GSE10643', 'GSE15577', 'GSE11538', 'GSE70861', 'GSE6583', 'GSE27551', 'GSE12619', 'GSE121225', 'GSE108070', 'GSE78713', 'GSE110079', 'GSE63128', 'GSE60960', 'GSE37118', 'GSE79681', 'GSE63372', 'GSE5622', 'GSE26983', 'GSE27550', 'GSE19603', 'GSE95202', 'GSE53308', 'GSE16765', 'GSE71855', 'GSE58620', 'GSE24177', 'GSE35258', 'GSE10670', 'GSE49418', 'GSE18666', 'GSE83136', 'GSE44655', 'GSE27549', 'GSE19700', 'GSE103398', 'GSE63522', 'GSE201609', 'GSE5620', 'GSE66369', 'GSE2268', 'GSE71237', 'GSE48474', 'GSE41935', 'GSE27548', 'GSE5623', 'GSE72050', 'GSE126373']):
    lables = {}
    
    study_metdata = dict.fromkeys(Studies,None)
    for file in tqdm(os.listdir(folder_path)):
        if file.endswith('.json'):
            file_path = os.path.join(folder_path, file)
            try:
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
                meta_data['Sample title:'] = data['Sample metadata:']['title']
                meta_data['Sample source_name:'] = data['Sample metadata:']['source_name_ch1']
                try:
                    meta_data['Sample characteristicts:'] = data['Sample metadata:']['characteristics_ch1']
                except KeyError:
                    meta_data['Sample characteristicts:'] = ''
                try:
                    meta_data['Sample extraction protocol:'] = data['Sample metadata:']['extract_protocol_ch1']
                except KeyError:
                    meta_data['Sample extraction protocol:'] = ''
                # meta_data['Samples extra info:'] = meta_data['Sample metadata:']
                # meta_data['Study metadata:'] = meta_data['Study metadata:']
                del meta_data['Sample metadata:']
                del meta_data['Study metadata:']
                if study_metdata[study_id] is None:
                    study_metdata[study_id] = {'Study metadata': data['Study metadata:']}
                study_metdata[study_id][sample_id] = meta_data
            except Exception as e:
                print(f"Error loading {file}: {str(e)}")


    x = 0
    for study in study_metdata:
        for sample in study_metdata[study]:
            with open(f'sample_metadata/{study}_{sample}.json', 'w') as handle:
                json.dump(study_metdata[study][sample], handle)

generate_metadata()