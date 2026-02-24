import numpy as np
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri, Formula
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr
import os
import sys
module_dir = './'
sys.path.append(module_dir)
from src.constants import *

def diff_exp_combine_tissues(treatments,save_dir,data_type,out_dir,samples=None, pure=False,tissue=None,filter_low_combination:int = 10):
    output_dir: str = f'{out_dir}'
    os.makedirs(output_dir,exist_ok=True)
    for treatment in treatments:
        try:
            output_filename = f"{tissue if tissue is not None else 'All-tissues'}_{treatment}"
            file_to_output = f"{output_dir}{output_filename}_genes.csv"
            print(f"--- Starting analysis for treatment: {treatment} across all tissues ---")
            # if os.path.isfile(file_to_output):
            #     print(f"+++ File already exists, skipping it to save time  +++")
            #     continue

            # 1. Load and prepare data and design files
            data = pd.read_csv(f'{save_dir}/{data_type}.csv', index_col=0)
            if data_type == '2_way_norm_og':
                data.columns = [col.split('_')[0] for col in data.columns]
            else:
                data.columns = [col.split('_')[-1] for col in data.columns]

            design = pd.read_csv(f'{CLUSTER_EXPLORATION_FIGURES_DIR}{EXPERIMENT_NAME}/labels.csv')

            # 2. Create masks for the target treatment and control samples
            if tissue:
                is_tissue = design['tissue'].str.contains(tissue, na=False)
            else:
                is_tissue = design['tissue'].str.contains('', na=False)
            is_treatment = design['treatment'].str.contains(treatment, na=False)
            is_only_treatment = design['treatment'].apply(lambda x: len(x) == len(treatment)+4 and treatment in x)
            if pure:
                is_treatment = is_only_treatment
            is_control = design['treatment'].str.contains("No stress", na=False)
            
            if samples is not None:
                is_study = design['sample_id'].apply(lambda x: x in samples)
                design_filtered = design[(is_treatment | is_control) & is_study & is_tissue].copy()
            else:
                design_filtered = design[(is_treatment | is_control) & is_tissue].copy()

            # 3. Synchronize samples between expression data and design metadata
            # Find the common set of sample IDs present in both dataframes
            common_samples = list(set(design_filtered['sample_id']) & set(data.columns))
            
            # Filter both dataframes to keep only the common samples
            design_filtered = design_filtered[design_filtered['sample_id'].isin(common_samples)]
            data_filtered = data[common_samples]
            
            # drop duplicate columns
            data_filtered= data_filtered.T.drop_duplicates().T
            design_filtered = design_filtered.drop_duplicates()

            design_filtered = design_filtered.sort_values(by='sample_id').reset_index(drop=True)
            design_filtered = design_filtered.groupby(['treatment','tissue']).filter(lambda x: len(x) >= filter_low_combination)
            data_filtered = data_filtered[design_filtered['sample_id']]
            

            print(f"Found and aligned {len(design_filtered)} samples for '{treatment}' vs. 'No stress'.")
            
            # 4. Create the final metadata DataFrame for the model
            metadata = design_filtered[['sample_id', 'treatment', 'tissue']].copy()
            del design_filtered
            metadata.rename(columns={'tissue': 'tissue'}, inplace=True) # Rename tissue
            
            # Explicitly create 'Control' and 'treatment' labels.
            # This ensures 'Control' will be the first level alphabetically.
            
            # Find the rows that are 'No stress'
            is_control_mask = metadata['treatment'].str.contains("No stress", na=False)
            
            # Default everything to 'treatment'
            metadata['Target'] = 'treatment'
            # Then, use the mask to label the 'Control' rows
            metadata.loc[is_control_mask, 'Target'] = 'Control'
            
            # Now we can drop the original treatment column if we want
            del metadata['treatment']
            
            single_tissue = False
            if len(set(metadata['tissue'])) == 1:
                single_tissue = True
                del metadata['tissue']

            if len(set(metadata['Target']))==1:
                print(f"+++ No enogh samples for this treatment need at least {filter_low_combination}, skiping it  +++")
                continue

            
            # 5. Import R libraries
            base = importr('base')
            stats = importr('stats')
            if CLUSTER_RUN:
                limma = importr('limma')#, lib_loc='/home/alex/R/x86_64-pc-linux-gnu-library/4.5/')
                writexl = importr('writexl')#, lib_loc='/home/alex/R/x86_64-pc-linux-gnu-library/4.5/')
            else:
                limma = importr('limma', lib_loc='/home/alex/R/x86_64-pc-linux-gnu-library/4.5/')
                writexl = importr('writexl', lib_loc='/home/alex/R/x86_64-pc-linux-gnu-library/4.5/')

            # 6. Convert pandas DataFrames to R objects
            with localconverter(ro.default_converter + pandas2ri.converter):
                r_data = ro.conversion.py2rpy(data_filtered)
                r_metadata = ro.conversion.py2rpy(metadata)
                genes = ro.StrVector(data_filtered.index.tolist())

                
            # Clean the 'Target' column in the R metadata DataFrame sto ensure valid R names
            print("Cleaning target levels to be valid R names...")
            # Get the Target column from the R dataframe
            target_col = r_metadata.rx2('Target') 
            # Clean it using R's make.names function
            clean_target_col = base.make_names(target_col)
            # Find the index of the 'Target' column and replace it with the cleaned version
            target_idx = r_metadata.names.index('Target')
            r_metadata[target_idx] = clean_target_col


            # 7. Create the model matrix
            if single_tissue:
                r_formula = Formula('~0 + Target')
            else:
                r_formula = Formula('~0 + Target + tissue')
            r_design = stats.model_matrix(r_formula, data=r_metadata)

            # 8. Fit the linear model — This should now work without errorssample_id
            print("Fitting linear model with tissue as a covariate...")
            fit = limma.lmFit(r_data, r_design)
            
            # 9. Define the contrast
            target_levels = base.levels(base.factor(r_metadata.rx2('Target')))
            # Make the contrast robust against different numbers of levels
            if len(target_levels) < 2:
                print(f"⚠️ Skipping contrast: Not enough treatment levels found for '{treatment}'.")
                continue
            # Use paste0 to create a syntactically valid contrast string for R
            contrast_str = f"Target{target_levels[-1]}-Target{target_levels[0]}"
            print(f"Making contrast: {contrast_str}")
            contrast_matrix = limma.makeContrasts(contrast_str, levels=r_design)
            
            # 10. Fit the contrast and compute stats
            fit2 = limma.contrasts_fit(fit, contrast_matrix)
            fit2 = limma.eBayes(fit2)

            # 11. Extract and save results
            print("Extracting results...")
            r_output = limma.topTreat(fit2, coef=1, genelist=genes, number=np.inf)
            writexl.write_xlsx(r_output, f"{output_dir}{output_filename}.xlsx")

            diff_exp_results = pd.read_excel(f"{output_dir}{output_filename}.xlsx")
            significant_genes = diff_exp_results.sort_values(by=['adj.P.Val'])
            # genes_to_save = significant_genes[significant_genes['adj.P.Val']<0.05][['ID','t','logFC','adj.P.Val']]
            genes_to_save = significant_genes[['ID','t','logFC','adj.P.Val']]
            genes_to_save.to_csv(f"{output_dir}{output_filename}_genes.csv", index=False)
            del genes_to_save
            del diff_exp_results
            print(f"Successfully completed analysis for '{treatment}'. Results saved.")

        except Exception as e:
            print(f"ERROR processing '{treatment}': {e}")