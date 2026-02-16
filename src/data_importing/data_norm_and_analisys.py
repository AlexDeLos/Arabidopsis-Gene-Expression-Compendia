import pandas as pd
import numpy as np
from collections import Counter
import os
import sys

# --- R Integration Imports ---
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri, numpy2ri
from rpy2.robjects.packages import importr
# Scikit-learn
from sklearn.preprocessing import RobustScaler

# Local imports
module_dir = './'
sys.path.append(module_dir)
from src.constants import *
# Assuming these helpers exist in your src folder
from src.data_importing.helpers.helpers import get_first_indexs, apply_KNN_impute, box_plot
from src.data_analisys.utils.cluster_exploration_utils import *

# --- Constants & Configuration ---
SAMPLE_STUDY_MAP = pd.read_csv(STORAGE_DIR+'/final_data/RMA_Microarray_Combined_sample_map.csv', index_col=0)

def get_study(sample: str):
    """Extracts StudyID from sample name."""
    try:
        sample_key = sample.split('.')[0]
        return SAMPLE_STUDY_MAP.loc[sample_key, 'StudyID']
    except KeyError:
        # Fallback if sample not found
        return "Unknown_Study"

def run_r_combat(df, batch_list, covar_df=None):
    """
    Memory-optimized wrapper to run R's sva::ComBat from Python.
    """
    print("  [R-ComBat] Initializing R interface...")
    
    # 1. Activate converters locally to avoid global pollution
    try:
        sva = importr('sva')
        stats = importr('stats')
        base = importr('base') # Access to R's garbage collector
    except Exception as e:
        raise ImportError("Could not import R 'sva' package.") from e

    # 2. Optimization: Ensure input is strictly float (prevents object-type overhead)
    # and Fortran-contiguous (R expects column-major, this might speed up conversion)
    print("  [R-ComBat] preparing input data...")
    if not np.issubdtype(df.values.dtype, np.floating):
        df = df.astype(float)
    
    # Create R object
    dat_r = numpy2ri.py2rpy(df.values)
    
    # Prepare batch vector
    batch_r = ro.StrVector(batch_list)
    
    # Prepare Model Matrix
    mod_r = ro.r('NULL')
    if covar_df is not None and not covar_df.empty:
        print(f"  [R-ComBat] Creating model matrix with covariates: {list(covar_df.columns)}")
        covar_r = pandas2ri.py2rpy(covar_df)
        formula_str = "~ " + " + ".join(covar_df.columns)
        formula = ro.Formula(formula_str)
        mod_r = stats.model_matrix(formula, data=covar_r)

    # 3. Run ComBat
    print("  [R-ComBat] Calling sva::ComBat...")
    combat_data_r = sva.ComBat(dat=dat_r, batch=batch_r, mod=mod_r, par_prior=True, prior_plots=False)

    # --- CRITICAL OPTIMIZATION START ---
    # Delete the R INPUT data immediately to free memory *before* allocating Python output
    print("  [R-ComBat] Cleaning up R input objects...")
    del dat_r
    del batch_r
    del mod_r
    
    # Force R Garbage Collection
    base.gc()
    # --- CRITICAL OPTIMIZATION END ---

    # 4. Convert back to Pandas
    print("  [R-ComBat] Converting result back to Python...")
    # This creates the Python copy. Since dat_r is gone, we have more room.
    combat_data_np = np.array(combat_data_r)
    
    # --- CRITICAL OPTIMIZATION 2 ---
    # Delete the R OUTPUT data immediately
    del combat_data_r
    base.gc()
    # -------------------------------

    result_df = pd.DataFrame(combat_data_np, index=df.index, columns=df.columns)
    
    return result_df

def run_preprocessing(plot_boxPlots=False, no_change=False):
    path = PROCESSED_DATA_FOLDER
    out_path = FILTERING_FIGURES
    os.makedirs(out_path, exist_ok=True)

    # --- 1. Load Data ---
    try:
        filtered_df = pd.read_csv(path+'filter.csv', index_col=0)
        print('Successfully loaded filtered data.')
    except FileNotFoundError:
        if no_change: raise
        
        print("Creating filter.csv from raw data...")
        big_df = pd.read_csv(COMBINED_DATA_OUTPUT_FILE, index_col=0)
        big_df = fuse_columns_by_sample(big_df)
        
        # Filter Rows (Genes) > 20% NaN
        nan_genes = big_df.isna().mean(axis=1) * 100
        filtered_df = big_df.loc[nan_genes <= 20]

        # Filter Cols (Samples) > 20% NaN
        nan_samples = filtered_df.isna().mean() * 100
        filtered_df = filtered_df[filtered_df.columns[nan_samples <= 20]]

        filtered_df.to_csv(path+'filter.csv')
    
    df_impute = filtered_df

    # --- 2. Batch Correction (ComBat) ---
    try:
        study_corrected_df = pd.read_csv(path+'study_corrected.csv', index_col=0)
    except FileNotFoundError:
        if no_change: raise
        
        print("Running Batch Correction...")
        all_studies = list(map(get_study, df_impute.columns))

        # Filter out studies with < 2 samples (ComBat requirement)
        study_counts = Counter(all_studies)
        single_sample_studies = {s for s, c in study_counts.items() if c < 2}

        if single_sample_studies:
            print(f"Removing studies with <2 samples: {single_sample_studies}")
        
        valid_cols = []
        valid_batches = []
        for col, study in zip(df_impute.columns, all_studies):
            if study not in single_sample_studies:
                valid_cols.append(col)
                valid_batches.append(study)

        df_impute = df_impute[valid_cols]
        
        # Run R ComBat
        study_corrected_df = run_r_combat(df_impute, valid_batches)
        study_corrected_df.to_csv(path + 'study_corrected.csv')

    return    
    # --- 3. Standard Scaler ---
    try:
        standardized_df = pd.read_csv(path+'/standardized.csv', index_col=0)
    except FileNotFoundError:
        if no_change: raise
        standardized_df = ((df_impute.T - df_impute.T.mean()) / df_impute.T.std()).T
        standardized_df.to_csv(path+'/standardized.csv')

    # --- 4. Robust Scaler ---
    try:
        robust_df = pd.read_csv(path+'/robust.csv', index_col=0)
    except FileNotFoundError:
        if no_change: raise
        scaler = RobustScaler()
        robust_df = pd.DataFrame(scaler.fit_transform(df_impute.T).T, 
                                 columns=df_impute.columns, index=df_impute.index)
        robust_df.to_csv(path+'/robust.csv')
    
    # --- 5. Two-way Norm (Median/IQR) ---
    try:
        double_norm = pd.read_csv(path+'/2_way_norm.csv', index_col=0)
    except FileNotFoundError:
        if no_change: raise            
        mat = study_corrected_df.to_numpy()
        q75, q25 = np.percentile(mat, [75, 25], axis=1, keepdims=True)
        iqr = q75 - q25
        # Avoid div by zero
        iqr[iqr == 0] = 1.0 
        norm = (mat - np.median(mat, axis=1, keepdims=True)) / iqr
        double_norm = pd.DataFrame(norm, columns=study_corrected_df.columns, index=study_corrected_df.index)
        double_norm.to_csv(path+'/2_way_norm.csv')

    # --- 6. Scalers on Corrected Data ---
    try:
        standardized_df_ = pd.read_csv(path+'/standardized+.csv', index_col=0)
    except FileNotFoundError:
        if no_change: raise
        standardized_df_ = ((study_corrected_df.T - study_corrected_df.T.mean()) / study_corrected_df.T.std()).T
        standardized_df_.to_csv(path+'/standardized+.csv')
    
    try:
        robust_df_ = pd.read_csv(path+'/robust+.csv', index_col=0)
    except FileNotFoundError:
        if no_change: raise
        scaler_ = RobustScaler()
        robust_df_ = pd.DataFrame(scaler_.fit_transform(study_corrected_df.T).T, 
                                  columns=study_corrected_df.columns, index=study_corrected_df.index)
        robust_df_.to_csv(path+'/robust+.csv')
    
    if plot_boxPlots:
        print("Generating Boxplots...")
        # Add your box_plot calls here...
    
    print('Preprocessing Done')

def _flatten_covariate(value):
    if isinstance(value, list):
        return '_'.join(map(str, value))
    return str(value)

def normalize_by_cov(df_small, cov='treatment'):
    labels = load_labels_study(LABELS_PATH)
    df = df_small
    
    # 1. Prepare Batches
    batches = list(map(get_study, df.columns))
    
    # 2. Prepare Covariates
    covariate_data = []
    for sample_col in df.columns:
        try:
            study_id, sample_id = sample_col.split('_', 1)
            info = labels.get(study_id, {}).get(sample_id, {})
            covariate_data.append({
                'tissue': _flatten_covariate(info.get('tissue', 'unknown')),
                'treatment': _flatten_covariate(info.get('treatment', 'unknown'))
            })
        except ValueError:
            covariate_data.append({'tissue': 'unknown', 'treatment': 'unknown'})

    covar_df = pd.DataFrame(covariate_data, index=df.columns)

    # 3. Diagnostic Print
    print("\n--- Confounding Check ---")
    check_df = covar_df.copy()
    check_df['batch'] = batches
    print(pd.crosstab(check_df['tissue'], check_df['batch']))
    print(pd.crosstab(check_df['treatment'], check_df['batch']))

    # 4. Handle Confounding
    # We want to keep 'cov' (e.g., treatment) and drop 'tissue' if it's not the target
    # R handles 'biological' variables vs 'adjustment' variables.
    # Here we just pass the one we want to PRESERVE to the model matrix.
    
    # If we want to preserve 'cov', we include it in the dataframe passed to R.
    # The others are effectively ignored by the model (not adjusted for).
    covar_to_use = covar_df[[cov]].copy()
    
    # Check if variable has at least 2 levels (otherwise model matrix fails)
    if covar_to_use[cov].nunique() < 2:
        print(f"Warning: Covariate '{cov}' has only 1 level. Running ComBat without covariates.")
        covar_to_use = None
    
    # 5. Run ComBat
    try:
        df_corrected = run_r_combat(df, batches, covar_df=covar_to_use)
        print("ComBat normalization successful.")
        return df_corrected
    except Exception as e:
        print(f"ComBat failed with covariates: {e}")
        print("Retrying without covariates...")
        return run_r_combat(df, batches)

def normalize_by_tissue_2():
    print("Running normalize_by_tissue_2...")
    df = pd.read_csv(PROCESSED_DATA_FOLDER+'imputed.csv', index_col=0)
    labels = load_labels_study(LABELS_PATH)
    labels = keys_upper(labels)

    labels_types = ['TREATMENT', 'TISSUE', 'MEDIUM']
    labels_df = make_df_from_labels(labels, labels_types)
    maps = get_label_map_new(df, labels_df)
    
    df_normalized_parts = []
    tissues = ["root", "leaf", "shoot", "rosette", "whole_plant", "callus", "seedling"]
    
    for tissue in tissues:
        print(f"Processing tissue: {tissue}")
        try:
            df_copy, maps_copy = get_df_and_maps(df, maps, 'TISSUE', tissue)
            if df_copy.empty:
                continue
            df_small = normalize_by_cov(df_copy)
            df_normalized_parts.append(df_small)
        except Exception as e:
            print(f"Skipping tissue {tissue}: {e}")

    if not df_normalized_parts:
        print("No tissues processed.")
        return

    df_normalized = pd.concat(df_normalized_parts, axis=1)
    df_normalized.to_csv(PROCESSED_DATA_FOLDER+'/tissue_normalized.csv')

    # Final batch correction on the merged result
    batches = list(map(get_study, df_normalized.columns))
    df_norm_2 = run_r_combat(df_normalized, batches)
    df_norm_2.to_csv(PROCESSED_DATA_FOLDER+'/tissue_normalized_2.csv')

if __name__ == '__main__':
    run_preprocessing()
    # normalize_by_tissue_2()