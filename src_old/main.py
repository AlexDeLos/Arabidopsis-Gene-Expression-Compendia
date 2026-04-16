from data_analisys.label_cluster_exploration import run_label_cluster_exploration
from data_analisys.diff_and_GSEA_pipeline import run_diff_exp_and_enrichment
from data_importing.import_GEOparse import import_data
from data_importing.data_norm_and_analisys import run_preprocessing
from meta_data_processing.label_generation import condense_labels
# RUN DATA IMPORTING
import sys
module_dir = './'
sys.path.append(module_dir)
from src.constants import *


Data_types = ['filter','study_corrected']# ,'tissue_normalized','tissue_normalized_2','robust', 'standardized', '2_way_norm',
#process and filter the data
run_preprocessing()
#RUN METADATA AND LABELING
condense_labels(in_folder='new_storage/processed_microarray_data/',saving_path=LABELS_PATH)

# RUN UMAP AND CLUSTER ANALISYS
run_label_cluster_exploration(0,Data_types)
# run_label_cluster_exploration(10,Data_types)
run_label_cluster_exploration(15,Data_types)

# RUN DIFF EXP and ENRICHMENT ANALISYS
run_diff_exp_and_enrichment(data_types =Data_types)