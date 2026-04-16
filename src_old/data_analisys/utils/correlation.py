import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency

def cramers_v(contingency_table):
    """
    Calculates Cramér's V for a contingency table.
    Cramér's V is a measure of association between two categorical variables.
    
    Args:
        contingency_table (pd.DataFrame): A contingency table (e.g., from pd.crosstab).
        
    Returns:
        float: The Cramér's V value, between 0 and 1.
    """
    # Get the chi-squared statistic
    # chi2_contingency returns: chi2, p-value, degrees_of_freedom, expected_frequencies
    chi2 = chi2_contingency(contingency_table)[0]
    
    # Get the total number of observations
    n = contingency_table.sum().sum()
    
    # Get the dimensions of the table (rows, columns)
    r, k = contingency_table.shape
    
    # Handle edge case where n is zero
    if n == 0:
        return 0.0
        
    # Calculate phi-squared
    phi2 = chi2 / n
    
    # Calculate the correction for the minimum dimension
    min_dim = min(r, k) - 1
    
    # Handle edge case where a variable is constant (min_dim == 0)
    if min_dim == 0:
        # Correlation is undefined or 0, depending on interpretation.
        # Returning 0.0 as there is no co-variance.
        return 0.0
        
    # Calculate Cramér's V
    v = np.sqrt(phi2 / min_dim)
    return v

def calculate_study_correlations(data_dict):
    """
    Calculates the correlation (Cramér's V) between a sample's study_id
    and its 'tissue' and 'treatment' variables.
    
    Args:
        data_dict (dict): The nested dictionary of study data with the structure:
                          {"study_id": {"sample_id": {"tissue": ..., "treatment": [...]}}}
            
    Returns:
        dict: A dictionary containing the correlation values for
              'study_vs_tissue' and 'study_vs_treatment'.
    """
    records = []
    # Flatten the nested dictionary into a list of records
    for study_id, samples in data_dict.items():
        for sample_id, details in samples.items():
            
            # --- Handle the 'treatment' list ---
            # Get the treatment list, default to empty list if key is missing
            treatment_list = details.get('treatment', [])
            
            # Ensure all items in the list are strings, then sort them.
            # Sorting ensures that ['A', 'B'] and ['B', 'A'] are treated as
            # the same category.
            sorted_treatments = sorted([str(t) for t in treatment_list])
            
            # Join them with an underscore to create a single categorical value
            treatment_str = "_".join(sorted_treatments)
            
            records.append({
                'study_id': study_id,
                'sample_id': sample_id,
                'tissue': details.get('tissue'),
                'treatment': treatment_str
            })
    
    # If no records were found, return NaNs
    if not records:
        return {'study_vs_tissue': np.nan, 'study_vs_treatment': np.nan}
        
    # Create a DataFrame from the flat list
    df = pd.DataFrame(records)
    
    # --- Correlation for Study vs. Tissue ---
    
    # We can only calculate correlation if there is variance in both variables
    if df['tissue'].nunique() <= 1 or df['study_id'].nunique() <= 1:
        study_tissue_corr = np.nan # Not enough variance to correlate
    else:
        # Create the contingency table
        study_tissue_contingency = pd.crosstab(df['study_id'], df['tissue'])
        # Calculate Cramér's V
        study_tissue_corr = cramers_v(study_tissue_contingency)
        
    # --- Correlation for Study vs. Treatment ---
    
    # Check for variance again
    if df['treatment'].nunique() <= 1 or df['study_id'].nunique() <= 1:
        study_treatment_corr = np.nan # Not enough variance to correlate
    else:
        # Create the contingency table
        study_treatment_contingency = pd.crosstab(df['study_id'], df['treatment'])
        # Calculate Cramér's V
        study_treatment_corr = cramers_v(study_treatment_contingency)
        
    return {
        'study_vs_tissue': study_tissue_corr,
        'study_vs_treatment': study_treatment_corr
    }
