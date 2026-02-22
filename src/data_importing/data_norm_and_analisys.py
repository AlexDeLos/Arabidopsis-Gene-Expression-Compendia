import pandas as pd
import numpy as np
from collections import Counter
from sklearn.impute import KNNImputer
import matplotlib.pyplot as plt
import os
import sys

# --- R Integration Imports ---
from rpy2.robjects.conversion import localconverter
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri, numpy2ri
from rpy2.robjects.packages import importr
# Scikit-learn

# Local imports
module_dir = './'
sys.path.append(module_dir)
from src.constants import *
from src.data_analisys.utils.cluster_exploration_utils import *

def plot_filtering_summary(df_before, df_after, output_path):
    """
    Creates a combined visualization showing the counts of genes/samples
    and a heatmap of missing values before and after filtering.
    """
    print(f"  [Plotting] Generating filtering summary visualization...")
    
    fig = plt.figure(figsize=(16, 12))
    fig.suptitle('Data Filtering Summary (0s converted to NaNs)', fontsize=20, fontweight='bold')

    # Create a grid layout: Top half for bar chart, Bottom half for heatmaps
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 2.5])

    # --- 1. Bar Chart: Gene and Sample Counts ---
    ax_bars = fig.add_subplot(gs[0, :])
    
    categories = ['Genes (Rows)', 'Samples (Columns)']
    before_counts = [df_before.shape[0], df_before.shape[1]]
    after_counts = [df_after.shape[0], df_after.shape[1]]

    x = np.arange(len(categories))
    width = 0.35

    ax_bars.bar(x - width/2, before_counts, width, label='Before Filtering', color='#ff9999', edgecolor='black')
    ax_bars.bar(x + width/2, after_counts, width, label='After Filtering', color='#66b3ff', edgecolor='black')

    ax_bars.set_ylabel('Count', fontsize=12)
    ax_bars.set_title('Total Genes and Samples Before vs After', fontsize=14)
    ax_bars.set_xticks(x)
    ax_bars.set_xticklabels(categories, fontsize=12)
    ax_bars.legend()

    # Add text labels on top of the bars
    for i, v in enumerate(before_counts):
        ax_bars.text(i - width/2, v + (max(before_counts)*0.02), f"{v:,}", ha='center', va='bottom', fontweight='bold')
    for i, v in enumerate(after_counts):
        ax_bars.text(i + width/2, v + (max(before_counts)*0.02), f"{v:,}", ha='center', va='bottom', fontweight='bold')

    # --- 2. Heatmap: Missingness Before ---
    ax_heat_before = fig.add_subplot(gs[1, 0])
    # Yellow/Light = Missing (True), Dark/Purple = Present (False)
    ax_heat_before.imshow(df_before.isna(), aspect='auto', cmap='viridis', interpolation='nearest')
    
    pct_missing_before = (df_before.isna().sum().sum() / df_before.size) * 100
    ax_heat_before.set_title(f'Missingness Before\n({pct_missing_before:.1f}% Total Missing)', fontsize=14)
    ax_heat_before.set_xlabel('Samples', fontsize=12)
    ax_heat_before.set_ylabel('Genes', fontsize=12)

    # --- 3. Heatmap: Missingness After ---
    ax_heat_after = fig.add_subplot(gs[1, 1])
    ax_heat_after.imshow(df_after.isna(), aspect='auto', cmap='viridis', interpolation='nearest')
    
    pct_missing_after = (df_after.isna().sum().sum() / df_after.size) * 100
    ax_heat_after.set_title(f'Missingness After\n({pct_missing_after:.1f}% Total Missing)', fontsize=14)
    ax_heat_after.set_xlabel('Samples', fontsize=12)
    ax_heat_after.set_ylabel('Genes', fontsize=12)

    plt.tight_layout(rect=[0, 0, 1, 0.96]) # Leave room for suptitle
    
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  [Plotting] Saved to {output_path}")

def apply_KNN_impute(df: pd.DataFrame, n_neighbors: int = 5) -> pd.DataFrame:
    """
    Imputes missing values (NaNs) in a pandas DataFrame using K-Nearest Neighbors.
    
    Parameters:
    - df: The input DataFrame containing missing values.
    - n_neighbors: The number of nearest neighbors to use for imputation.
    
    Returns:
    - A new DataFrame with the missing values imputed, preserving index and columns.
    """
    print(f"  [KNN Impute] Imputing missing values using {n_neighbors} neighbors...")
    
    # Initialize the imputer
    imputer = KNNImputer(n_neighbors=n_neighbors, weights="uniform")
    
    # Fit and transform the data
    # Note: KNNImputer returns a numpy array, so we must reconstruct the DataFrame
    imputed_array = imputer.fit_transform(df)
    
    # Reconstruct the DataFrame with original indices and columns
    df_imputed = pd.DataFrame(imputed_array, index=df.index, columns=df.columns)
    
    print("  [KNN Impute] Imputation complete.")
    return df_imputed

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
        
        # 1. Force all columns to be strings so R can interpret them as categorical factors
        covar_df = covar_df.astype(str)
        
        # 2. Use the localconverter context to safely translate Pandas -> R
        with localconverter(ro.default_converter + pandas2ri.converter):
            covar_r = ro.conversion.py2rpy(covar_df)
            
        # 3. Build the formula
        formula_str = "~ " + " + ".join(covar_df.columns)
        formula = ro.Formula(formula_str)
        
        # 4. Generate the model matrix
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

def run_preprocessing(no_change=False):
    path = PROCESSED_DATA_FOLDER
    out_path = FILTERING_FIGURES
    os.makedirs(out_path, exist_ok=True)

    # --- 1. Load and Filter Data ---
    try:
        filtered_df = pd.read_csv(path+'filter.csv', index_col=0) 
        print('Successfully loaded filtered data.')
    except FileNotFoundError:
        if no_change: raise
        
        print("Creating filter.csv from raw data...")
        big_df = pd.read_csv(COMBINED_DATA_OUTPUT_FILE, index_col=0)
        big_df = fuse_columns_by_sample(big_df)
        
        # --- NEW: Treat 0 as NaN ---
        # This converts all exactly 0.0 values to NaN so they get filtered and imputed
        big_df = big_df.replace(0, np.nan)
        
        # Filter Rows (Genes) > 20% NaN (which now includes 0s)
        nan_genes = big_df.isna().mean(axis=1) * 100
        filtered_df = big_df.loc[nan_genes <= 20]

        # Filter Cols (Samples) > 20% NaN (which now includes 0s)
        nan_samples = filtered_df.isna().mean() * 100
        filtered_df = filtered_df[filtered_df.columns[nan_samples <= 20]]
        
        plot_output = os.path.join(out_path, 'filtering_missingness_summary.png')
        plot_filtering_summary(big_df, filtered_df, plot_output)
        filtered_df.to_csv(path+'filter.csv')
    
    # --- 2. NEW: KNN Imputation ---
    try:
        df_impute = pd.read_csv(path+'imputed.csv', index_col=0)
        print('Successfully loaded imputed data.')
    except FileNotFoundError:
        print("Running KNN Imputation on missing/zero values...")
        # Using the apply_KNN_impute function imported at the top of your script
        df_impute = apply_KNN_impute(filtered_df)
        df_impute.to_csv(path+'imputed.csv')

    # --- 3. Batch Correction (ComBat) ---
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

        df_impute_combat = df_impute[valid_cols]
        
        # Run R ComBat
        study_corrected_df = run_r_combat(df_impute_combat, valid_batches)
        study_corrected_df.to_csv(path + 'study_corrected.csv')

    return

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

def normalize_all_with_covariates():
    print("Running ComBat on full dataset with biological covariates...")
    
    # 1. Load data
    df = pd.read_csv(PROCESSED_DATA_FOLDER+'filter.csv', index_col=0)
    labels = load_labels_study(LABELS_PATH)

    # 2. Get labels
    labels_types = LABELS
    labels_df = make_df_from_labels(labels, labels_types)
    
    def combine_label_lists(series):
        combined = set()
        for val in series:
            # If the value is a list (e.g., ['heat stress']), add its elements
            if isinstance(val, list):
                combined.update([str(v) for v in val if pd.notna(v)])
            # If it's a single string/value, add it directly
            elif pd.notna(val):
                combined.add(str(val))
        # Return as a unique, sorted list to satisfy the requirement
        return sorted(list(combined))

    # Apply the fusion. This guarantees strictly unique indices.
    labels_df = labels_df.groupby(level=0).agg(combine_label_lists)
    
    # 3. Align DataFrame and metadata maps perfectly
    cols_in_maps = [c for c in df.columns if c in labels_df.index]
    df = df[cols_in_maps]
    maps_aligned = labels_df.loc[cols_in_maps]
    
    # 4. Extract batches (Study IDs)
    batches = [get_study(c) for c in df.columns]
    
    # 5. Build the Covariate DataFrame containing your labels
    # R's ComBat requires strings for categorical factors, not Python lists.
    # We flatten the fused lists with an underscore for the final model matrix.
    covariates_df = pd.DataFrame({
        'treatment': maps_aligned['treatment'].apply(lambda x: '_'.join(x)).tolist(),
        'tissue': maps_aligned['tissue'].apply(lambda x: '_'.join(x)).tolist(),
        'medium': maps_aligned['medium'].apply(lambda x: '_'.join(x)).tolist()
    }, index=df.columns)
    
    # 6. Run ComBat
    df_normalized = run_r_combat(df, batch_list=batches, covar_df=covariates_df)
    
    # 7. Save results
    df_normalized.to_csv(PROCESSED_DATA_FOLDER+'/fully_normalized_with_covariates.csv')
    print("Finished! Saved to fully_normalized_with_covariates.csv")
    
    
if __name__ == '__main__':
    run_preprocessing()
    normalize_all_with_covariates()