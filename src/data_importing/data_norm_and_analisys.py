import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from inmoose.pycombat import pycombat_norm
import os
from sklearn.preprocessing import RobustScaler
import sys
module_dir = './'
sys.path.append(module_dir)
from src.constants import *
from src.data_importing.helpers import get_first_indexs, apply_KNN_impute,box_plot
from src.data_analisys.utils.cluster_exploration_utils import *
from src.data_importing.helpers import find_and_plot_missing_genes

def get_study(sample: str):
    return sample.split('_')[0]

def run_preprocessing(
    plot_nan = True,
    plot_boxPlots = False,
    no_change = False):

    path = PROCESSED_DATA_FOLDER

    out_path = FILTERING_FIGURES

    # Create directory if it doesn't exist
    os.makedirs(out_path, exist_ok=True)

    try:
        # raise FileNotFoundError
        filtered_df = pd.read_csv(path+'filter.csv', index_col=0)
        # filtered_df_og = pd.read_csv('./data/downloads/old_plocessed_data/'+'filter.csv', index_col=0)
        # filtered_df_og.columns = list(map(lambda x: x.split('_')[2]+'_'+x.split('_')[0],filtered_df_og.columns))
        print('succesfully loaded data')
    except FileNotFoundError:
        if no_change:
            raise FileNotFoundError
        
        big_df = pd.read_csv(COMBINED_DATA_OUTPUT_FILE,index_col=0)
        big_df = fuse_columns_by_sample(big_df)
        # filter the data on 20%
        nan_percentage_genes = big_df.isna().mean(axis=1) * 100

        # Filter rows where NaN percentage <= 20%
        filtered_genes = nan_percentage_genes[nan_percentage_genes <= 20].index

        # Keep only the filtered rows
        filtered_df = big_df.loc[filtered_genes]


        # Calculate the percentage of NaN values in each column
        nan_percentage_samples = filtered_df.isna().mean() * 100

        # Filter columns where NaN percentage <= 20%
        filtered_columns = nan_percentage_samples[nan_percentage_samples <= 20].index

        # Keep only the filtered columns
        filtered_df = filtered_df[filtered_columns]

        filtered_df.to_csv(path+'filter.csv')
    
    print('data loaded')
    big_df = filtered_df
    find_and_plot_missing_genes(list(big_df.index), out_opath=FIGURES_DIR)

    matrix = big_df.to_numpy()
    matrix_nan = big_df.isna().to_numpy()

    chromosomes = ['1','2','3','4','5']
    indices:list[int] = get_first_indexs(big_df.index,chromosomes)

    if plot_nan:
        row_nan_count = big_df.isna().sum(axis=1)
        col_nan_count = big_df.isna().sum(axis=0)
        
        plt.imshow(matrix_nan, cmap='hot', interpolation='nearest')
        plt.savefig(out_path+'matrix.svg')
        plt.close()
        # Create directory if it doesn't exist
        os.makedirs(out_path+'row', exist_ok=True)
        os.makedirs(out_path+'col', exist_ok=True)
        for i,_ in enumerate(indices):
            min_var = indices[i]
            try:
                max_var = indices[i+1]
            except:
                max_var = len(matrix)
            plt.bar(range(len(row_nan_count.index[min_var:max_var])),row_nan_count.values[min_var:max_var])
            plt.ylim(0,1850)
            plt.xlabel('Genes')
            plt.ylabel('Missing number of data points')
            plt.savefig(out_path+'row/row_dis'+chromosomes[i]+'.svg')
            plt.close()

        plt.bar(range(len(col_nan_count.index)),col_nan_count.values)
        plt.xlabel('Experiments')
        plt.ylabel('Missing number of data points')
        plt.savefig(out_path+'col/col_dis.svg')
        plt.close()


    # With simple replacement
    print('plotted Nans', plot_nan)
    np.nan_to_num(matrix,copy=False)
    print('nans filled with 0')

    for column in big_df:
        # print(big_df[column])
        if (big_df[column]>300).any():
            big_df[column] = np.log2(big_df[column]+1)

    print('applied log2')

    study_map = list(map(get_study,big_df.columns))

    # KNN Impute
    try: 
        print('reading file from: ' + path+'imputed.csv')
        df_impute:pd.DataFrame = pd.read_csv(path+'/imputed.csv', index_col=0)
    except FileNotFoundError:
        if no_change:
            raise FileNotFoundError
        #! is this a good imputing?
        #? it might not work due to the other samples in the study probably also are missing the data we are interested in...
        # print('Could not read file, running KNN impute')
        # df_impute:pd.DataFrame = apply_KNN_impute(big_df,5)
        # print('KNN impute ran, saving file')
        # df_impute.to_csv(path+'imputed.csv')
        # print('file saved at: ' + path+'imputed.csv')
        
        ### BETTER IMPUTE FUNCTION?
        # 1. Get a unique list of all study IDs from the column names
        # We split 'StudyId_SampleId' by '_' and take the first part
        study_ids = big_df.columns.str.split('_').str[0].unique()
        
        print(f"Found {len(study_ids)} studies. Starting study-wise imputation...")

        # 2. Create a list to hold the imputed data for each study
        imputed_studies_list = []

        # 3. Loop through each study, impute its data, and add it to the list
        for study_id in study_ids:
            # Select columns belonging to the current study
            study_cols = [col for col in big_df.columns if col.startswith(f"{study_id}_")]
            study_df = big_df[study_cols]
            
            # Apply KNN imputation only on this subset of data
            print(f"  -> Imputing data for study: {study_id} ({len(study_cols)} samples)")
            imputed_study_df = apply_KNN_impute(study_df, 5)
            imputed_studies_list.append(imputed_study_df)

        # 4. Concatenate all the imputed dataframes back into a single dataframe
        # axis=1 combines them column-wise
        df_impute = pd.concat(imputed_studies_list, axis=1)

        # Optional but recommended: ensure the final column order matches the original
        df_impute = df_impute[big_df.columns]
        df_impute.to_csv(path+'imputed.csv')
        print('file saved at: ' + path+'imputed.csv')


    # NORMALIZE
    #? Apply batch correction
    try:
        study_corrected_df:pd.DataFrame = pd.read_csv(path+'study_corrected.csv',index_col=0)
    except FileNotFoundError:
        if no_change:
            raise FileNotFoundError
        study_map = list(map(get_study,df_impute.columns))
        # raise ValueError("Running in the cluster")
        d = dict([(y,x+1) for x,y in enumerate(sorted(set(study_map)))])
        batches = []
        for el in study_map:
            batches.append(d[el])
        study_corrected_df:pd.DataFrame = pycombat_norm(df_impute, batches)
        study_corrected_df.to_csv(path+'study_corrected.csv')

    try:
        standardized_df = pd.read_csv(path+'/standardized.csv', index_col=0)
    except FileNotFoundError:
        if no_change:
            raise FileNotFoundError
        standardized_df = ((df_impute.T - df_impute.T.mean()) / df_impute.T.std()).T
        standardized_df.to_csv(path+'/standardized.csv')
    try:
        robust_df = pd.read_csv(path+'/robust.csv', index_col=0)
    except FileNotFoundError:
        if no_change:
            raise FileNotFoundError
        scaler = RobustScaler()
        robust_df = pd.DataFrame(scaler.fit_transform(df_impute), columns=df_impute.columns,index = df_impute.index)
        robust_df.to_csv(path+'/robust.csv')
    
    try:
        double_norm = pd.read_csv(path+'/2_way_norm.csv', index_col=0)
    except FileNotFoundError:
        if no_change:
            raise FileNotFoundError            
        mat = study_corrected_df.to_numpy()
        q75, q25 = np.percentile(mat, [75 ,25],axis=1,keepdims=True)
        iqr = q75 - q25
        norm = (mat - np.median(mat,axis=1,keepdims=True))/iqr
        double_norm = pd.DataFrame(norm, columns=study_corrected_df.columns,index=study_corrected_df.index)
        double_norm.to_csv(path+'/2_way_norm.csv')

# double_norm
    try:
        standardized_df_ = pd.read_csv(path+'/standardized+.csv', index_col=0)
    except FileNotFoundError:
        if no_change:
            raise FileNotFoundError
        standardized_df_ = ((study_corrected_df.T - study_corrected_df.T.mean()) / study_corrected_df.T.std()).T
        standardized_df_.to_csv(path+'/standardized+.csv')
    try:
        robust_df_ = pd.read_csv(path+'/robust+.csv', index_col=0)
    except FileNotFoundError:
        if no_change:
            raise FileNotFoundError
        scaler_ = RobustScaler()
        robust_df_ = pd.DataFrame(scaler_.fit_transform(study_corrected_df), columns=study_corrected_df.columns,index=study_corrected_df.index)
        robust_df_.to_csv(path+'/robust+.csv')
    
    
    if plot_boxPlots:
        box_plot(df_impute,120, out_path+'/uncorrected_box_plots/')
        box_plot(study_corrected_df,120, out_path+'/study_corrected_box_plots/')
        box_plot(double_norm,120, out_path+'/2_way_norm_plots/')
        box_plot(robust_df,120, out_path+'/robust/')
        box_plot(standardized_df,120,out_path+'/standardized/')
        box_plot(robust_df_,120, out_path+'/robust+/')
        box_plot(standardized_df_,120,out_path+'/standardized+/')
    print('Done')

def _flatten_covariate(value):
    """
    Converts a potential list (or other non-hashable type) 
    into a single, hashable string.
    """
    if isinstance(value, list):
        return '_'.join(map(str, value))
    return value

def normalize_by_cov(df_small,cov = 'treatment'):
    labels = load_labels_study(LABELS_PATH)
    # df = pd.read_csv(PROCESSED_DATA_FOLDER+'imputed.csv', index_col=0)
    # df_small = pd.concat([df.filter(like='GSE34188', axis=1),df.filter(like='GSE58620', axis=1)],axis =1)
    df = df_small
    study_map = list(map(get_study,df.columns))
    
    d = dict([(y,x+1) for x,y in enumerate(sorted(set(study_map)))])
    batches = []
    for el in study_map:
        batches.append(d[el])

    covariate_data = []
    for sample_col in df.columns:
        try:
            study_id, sample_id = sample_col.split('_', 1)
        except ValueError:
            covariate_data.append({'tissue': None, 'treatment': None})
            continue

        try:
            info = labels[study_id][sample_id]
            covariate_data.append({
                'tissue': _flatten_covariate(info['tissue']),
                'treatment': _flatten_covariate(info['treatment'])
            })
        except KeyError:
            covariate_data.append({'tissue': None, 'treatment': None})

    covar_mod = pd.DataFrame(covariate_data, index=df.columns)

    # --- START: DIAGNOSTIC ---
    # This will print a table showing the overlap.
    # Look for any row (like 'whole_plant') that has a count in ONLY ONE
    # column (batch). This confirms it is confounded.
    print("--- Checking for Confounding Variables ---")
    check_df = covar_mod.copy()
    check_df['batch'] = batches # Add the batch list for comparison
    
    print("\n--- Tissue vs. Batch Crosstab ---")
    print(pd.crosstab(check_df['tissue'], check_df['batch']))
    
    print("\n--- Treatment vs. Batch Crosstab ---")
    print(pd.crosstab(check_df['treatment'], check_df['batch']))
    print("--------------------------------------------")
    # --- END: DIAGNOSTIC ---

    
    # By default, assume we have to drop both.
    confounded_vars = ['treatment','tissue']
    confounded_vars.remove(cov)
    covar_mod_fixed = covar_mod.drop(columns=confounded_vars, errors='ignore')

    # If dropping those columns results in an empty DataFrame,
    # we must pass 'None' to pycombat_norm.
    if covar_mod_fixed.empty:
        print("All covariates were confounded. Running ComBat without covariates.")
        covar_mod_fixed = None
    else:
        print(f"Removed confounded variables: {confounded_vars}")
    if len(covar_mod_fixed.value_counts()) ==1:
        confounded_vars = None
    # Pass the fixed covariate model (or None) to the function.
    try:
        df_corrected = pycombat_norm(df, 
                                    batches, 
                                    covar_mod=covar_mod_fixed, 
                                    na_cov_action='fill')
    except ValueError:
        df_corrected = pycombat_norm(df, 
                                    batches)
    print("ComBat normalization successful.")
    print(df_corrected.head())
    return df_corrected

def normalize_by_tissue_2():
    # get all samples grouped by tissue

    df = pd.read_csv(PROCESSED_DATA_FOLDER+'imputed.csv', index_col=0)
    labels = load_labels_study(LABELS_PATH)
    labels = keys_upper(labels)
    # Fuse the labels

    labels_types = ['TREATMENT','TISSUE','MEDIUM']
    labels_df = make_df_from_labels(labels, labels_types)
    maps = get_label_map_new(df,labels_df)
    df_normalized = pd.DataFrame()
    for tissue in ["root", "leaf", "shoot", "rosette", "whole_plant", "callus", "seedling"]: # "flower","bud","silique","seed", 
        df_copy, maps_copy = get_df_and_maps(df,maps,'TISSUE',tissue)
        df_small = normalize_by_cov(df_copy)
        df_normalized= pd.concat([df_normalized,df_small],axis=1)

    # combat norm using treatment as a covarience
    df_normalized.to_csv(PROCESSED_DATA_FOLDER+'/tissue_normalized.csv')


    study_map = list(map(get_study,df_normalized.columns))
    batches = []
    d = dict([(y,x+1) for x,y in enumerate(sorted(set(study_map)))])
    for el in study_map:
        batches.append(d[el])
    df_norm_2 = pycombat_norm(df_normalized,batches)
    df_norm_2.to_csv(PROCESSED_DATA_FOLDER+'/tissue_normalized_2.csv')
    # return df_normalized


if __name__ == '__main__':
    run_preprocessing()
    normalize_by_tissue_2()