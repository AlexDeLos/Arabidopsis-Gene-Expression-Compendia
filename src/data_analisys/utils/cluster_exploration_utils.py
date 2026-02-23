import json
import pandas as pd
import numpy as np
import os
import math
import matplotlib.pyplot as plt
import re
from collections import defaultdict
import sys
module_dir = './'
sys.path.append(module_dir)

from src.constants import *
def to_int(x:list,name:str,path:str):
    new_x = []
    seen = {}
    floor = 1
    for i in x:
        if type(i) == list:
            # we need to find a way to deal with lists
            i.sort()
            i = str(i)
        if type(i) == dict:
            i = '?'
        if i in seen:
            new_x.append(seen[i])
        elif i is None:
            new_x.append(0)
        else:
            new_x.append(floor)
            seen[i]= floor
            floor = floor+1
    save_json(object=seen,name=f'{name}_dict',path=path)
    return new_x

def to_int_(x:list):
    new_x = []
    seen = {}
    floor = 1
    for i in x:
        if type(i) == list:
            # we need to find a way to deal with lists
            i.sort()
            i = str(i)
        if type(i) == dict:
            i = '?'
        if i in seen:
            new_x.append(seen[i])
        elif i is None:
            new_x.append(0)
        else:
            new_x.append(floor)
            seen[i]= floor
            floor = floor+1
    return new_x,seen


def get_study(sample: str):
    """Extracts StudyID from sample name."""
    try:
        sample_key = sample.split('.')[0]
        return SAMPLE_STUDY_MAP.loc[sample_key, 'StudyID']
    except KeyError:
        # Fallback if sample not found
        return "Unknown_Study"
def get_studies(df):
    studies=[]
    for el in df.columns:
        studies.append(get_study(el))
    return studies


def get_label_map_new(data_df,labels_df):
    labels_map = {}
    for sample,study in zip(data_df.columns,get_studies(data_df)):
        sample=sample.split('.')[0].upper() #TODO: where does this lower case come from? gsm77182
        temp = dict(labels_df.loc[sample])
        for i in temp:
            temp[i] = str(temp[i])
        labels_map[sample] = temp
        labels_map[sample]['study'] = study
    return labels_map

def add_map(maps:dict,map:list,name:str):
    assert(len(maps)==len(map))
    for i,samples in enumerate(maps):
        maps[samples][name] = str(map[i])
    return maps
def get_map(maps:dict,name:str):
    ret = []
    for i,samples in enumerate(maps):
        ret.append(maps[samples][name])
    return ret


def flatten_extend(matrix):
    flat_list = []
    for row in matrix:
        if type(row) is list:
            flat_list.extend(row)
        else:
            flat_list = matrix
    return flat_list

def sanatize_labels(labels):
    for study in labels:
        for sample in labels[study]:
            if labels[study][sample]['treatment'] is None:
                labels[study][sample]['treatment'] = [None]
            
            labels[study][sample]['treatment'] = flatten_extend(labels[study][sample]['treatment'])
            for index in range(len(labels[study][sample]['treatment'])):
                st = False
                if type(labels[study][sample]['treatment']) is str:
                    st= True
                    labels[study][sample]['treatment'] = [labels[study][sample]['treatment']]
                treatment = labels[study][sample]['treatment'][index]
                if treatment is None:
                    if st:
                        break
                    continue
                if treatment == 'null':
                    labels[study][sample]['treatment'][index] = None
                    if st:
                        break
                    continue
                labels[study][sample]['treatment'][index] = labels[study][sample]['treatment'][index].replace(' ','_')
                labels[study][sample]['treatment'][index] = labels[study][sample]['treatment'][index].replace('_stress','')
                if 'nacl' in treatment:
                    labels[study][sample]['treatment'][index] = 'salt'
                if 'salt' in treatment:
                    labels[study][sample]['treatment'][index] = 'salt'
                if st:
                    break
    return labels

def save_json(object,name, path):
    try:
        with open(f'{path}/{name}.json', 'w') as handle:
            json.dump(object, handle)
    except:
        pass

#cluster only the data of the most active genes

def make_df_from_labels(data_dict,col):
    records = []
# Iterate through the dictionary structure
    for study_id, samples in data_dict.items():
        for sample_id, sample_data in samples.items():
            # Create a new record with study_id included
            record = {
                'study_id': study_id.upper(),
                'sample_id': sample_id.upper(),
                **sample_data  # Unpack all sample data
            }
            records.append(record)

    # Create DataFrame
    df = pd.DataFrame(records)

    # Set sample_id as index
    df.set_index('sample_id', inplace=True)

    # Select and order the columns you want
    # df = df[['study_id', 'control', 'treatment', 'tissue', 'agar']]
    return df[col]


def get_incomplete_studies(maps,sample_index,label):
    mask = list(map(lambda x: x is None,maps[label]))
    studies = np.unique(np.array(maps['study'])[list(mask)])
    samples = []
    for s in studies:
        index = maps['study'].index(s)
        samples.append(sample_index[index])
    return samples

def load_labels_study(path):
    labels = {}
    for file in os.listdir(path):
        if file.endswith('.json'):
            file_path = os.path.join(path, file)
            study = file_path.split('/')[-1].split('.')[0]
            if re.search('GSM',study):
                continue
            with open(file_path, 'r') as file:
                data = json.load(file)
            labels[study] = data
    return labels


def plot_pie_chart(data,path):
    sqrt = math.ceil(math.sqrt(len(data)))
    fig, axes = plt.subplots(sqrt, sqrt, figsize=(12, 10))
    # if sqrt*sqrt > len(data):
    #     to_remove = sqrt*sqrt- len(data)
    #     for i in range(to_remove):
    #         fig.delaxes(axes[sqrt-1,(sqrt-1)-i]) # The indexing is zero-based here
    axes = axes.ravel()
    for i, (key, values) in enumerate(data.items()):
        values = list(map(lambda x: str(x), values))
        unique, counts = np.unique(values, return_counts=True)
        axes[i].pie(counts, labels=unique, autopct='%1.1f%%', startangle=90)
        axes[i].set_title(f'Label Distribution - {key}')
        axes[i].axis('equal')

    plt.tight_layout()
    plt.savefig(f'{path}/pie.svg')
    plt.close()

def fuse_columns_by_sample(df):
    """
    Identifies columns with patterns and fuses them based on common name2:
    - Fuses {name1}_{name2}.1 with {name1}_{name2}
    - Also fuses columns that share the same name2 part
    
    Parameters:
    df (pd.DataFrame): Input dataframe
    
    Returns:
    pd.DataFrame: Dataframe with fused columns
    """
    
    result_df = df.copy()
    all_columns = result_df.columns.tolist()
    columns_to_drop = []
    fused_pairs = []
    
    # Pattern for .1 duplicates
    pattern_duplicate = re.compile(r'^(.+)\.1$')
    
    # Pattern to extract name2 (the part after last underscore)
    def extract_name_parts(column_name):
        """Extract name1 and name2 from column name"""
        if '_' in column_name:
            # Split by underscore and get the last part as name2
            parts = column_name.split('_')
            name1 = '_'.join(parts[:-1])  # Everything before last underscore
            name2 = parts[-1]  # Last part after last underscore
            # Remove .1 suffix if present
            if name2.endswith('.1'):
                name2 = name2[:-2]
            return name1, name2
        return None, None
    
    # First pass: handle .1 duplicates (original behavior)
    for col in all_columns:
        match = pattern_duplicate.match(col)
        if match:
            base_col_name = match.group(1)
            if base_col_name in all_columns and base_col_name not in columns_to_drop:
                result_df[base_col_name] = result_df[[base_col_name, col]].mean(axis=1, skipna=True)
                columns_to_drop.append(col)
                fused_pairs.append((base_col_name, col))
    
    # Update columns list after first pass
    current_columns = [col for col in result_df.columns if col not in columns_to_drop]
    
    # Second pass: group columns by name2 and fuse them
    name2_groups = defaultdict(list)
    
    # Group columns by their name2 part
    for col in current_columns:
        name1, name2 = extract_name_parts(col)
        if name2 is not None:
            name2_groups[name2].append(col)
    
    # Fuse columns that share the same name2
    for name2, columns in name2_groups.items():
        if len(columns) > 1:
            print(f"Fusing columns with same name2 '{name2}': {columns}")
            
            # Use the first column as the base for fusion
            base_col = columns[0]
            cols_to_fuse = columns[1:]
            
            # Calculate average of all columns with same name2
            result_df[base_col] = result_df[columns].mean(axis=1, skipna=True)
            
            # Mark other columns for removal
            columns_to_drop.extend(cols_to_fuse)
            
            for col in cols_to_fuse:
                fused_pairs.append((base_col, col))
    
    # Drop all marked columns
    result_df = result_df.drop(columns=columns_to_drop)
    
    # Print summary
    if fused_pairs:
        print(f"\nFused {len(fused_pairs)} column pairs:")
        for base, duplicate in fused_pairs:
            print(f"  '{base}' + '{duplicate}' -> '{base}'")
    else:
        print("No column pairs found for fusion")
    
    return result_df

def apply_mask_to_maps(maps,mask):
    maps_ = maps.copy()
    for map in maps:
        map_ = list(compress(maps_[map], mask))
        maps_[map] = map_
    return maps_
def pca_variance(robust_df,path):
    pca = PCA()
    #
    # Determine transformed features
    #
    X_train_pca = pca.fit_transform(robust_df.T)
    #
    # Determine explained variance using explained_variance_ration_ attribute
    #
    exp_var_pca = pca.explained_variance_ratio_
    #
    # Cumulative sum of eigenvalues; This will be used to create step plot
    # for visualizing the variance explained by each principal component.
    #
    cum_sum_eigenvalues = np.cumsum(exp_var_pca)
    #
    # Create the visualization plot
    #
    plt.bar(range(0,len(exp_var_pca)), exp_var_pca, alpha=0.5, align='center', label='Individual explained variance')
    plt.step(range(0,len(cum_sum_eigenvalues)), cum_sum_eigenvalues, where='mid',label='Cumulative explained variance')
    plt.ylabel('Explained variance ratio')
    plt.xlabel('Principal component index')
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig(f'{path}/test.svg')
    plt.close()


def plot_var(df,path):
    plt.plot(np.sort(np.array(df.var())))
    # plt.plot(np.array(df.var()))
    plt.savefig(f'{path}/gene_var.svg')
    plt.close()
    
def get_df_and_maps(df:pd.DataFrame,maps:dict,key:str,label:str):
    df_copy = df.copy()
    cols_to_keep = []
    samples = []
    maps_copy = maps.copy()
    for col in df.columns:
        col_sample_id = col.split('_')[-1]
        if maps[col_sample_id][key] == label:
            cols_to_keep.append(col)
            samples.append(col_sample_id)

    maps_copy = {key: maps_copy[key] for key in samples}
    df_copy = df_copy[cols_to_keep]
    return df_copy, maps_copy
