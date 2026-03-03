import pandas as pd
import os
import sys
import gc
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotly.express as px

from sklearn.metrics import silhouette_score, adjusted_rand_score
from sklearn.cluster import KMeans, MiniBatchKMeans
from sklearn.neighbors import KNeighborsClassifier
from sklearn.preprocessing import MultiLabelBinarizer
from sklearn.model_selection import cross_val_score
from scipy.stats import chi2_contingency
from sklearn.linear_model import LinearRegression

module_dir = './'
sys.path.append(module_dir)

from src.constants import *
from src.data_analisys.bulk_transformer import get_bulkformer_embeddings
from src.data_analisys.cluster_explotation import plot_metrics_comparison,variance_explained_by_label,calculate_asw_batch_within_biology
from src.data_analisys.utils.cluster_exploration_utils_2 import (
    prepare_data_structure, align_labels_to_data, 
    run_pca, run_umap, run_tsne
)
from src.data_analisys.utils.cluster_exploration_utils import *



# --- (plot_interactive_projection function remains exactly the same) ---

def run_exploration_on_dataframe(
    data_df: pd.DataFrame, 
    labels_dict: dict, 
    experiment_name: str,
    output_folder: str,
    use_bulkformer: bool = False,         # <-- NEW FLAG
    gene_length_dict: dict = None,        # <-- REQUIRED FOR BULKFORMER
    target_vocab: list = None,            # <-- REQUIRED FOR BULKFORMER
    ortholog_map: dict = None             # <-- REQUIRED FOR PLANT DATA
):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # 1. PREPROCESSING / EMBEDDING
    if use_bulkformer:
        print(f"  >>> Extracting BulkFormer Embeddings for {experiment_name}...")
        # Transpose data_df (Genes x Samples) to (Samples x Genes) for the model
        embeddings_df = get_bulkformer_embeddings(
            count_df=data_df.T, 
            gene_length_dict=gene_length_dict,
            target_vocab=target_vocab,
            ortholog_map=ortholog_map
        )
        # Transpose back to (Latent_Dims x Samples) for the existing alignment logic
        df_aligned = embeddings_df.T 
    else:
        print(f"  >>> Using standard PCA preprocessing for {experiment_name}...")
        df_aligned = prepare_data_structure(data_df)

    categories = ['treatment', 'tissue', 'medium','study_id']
    
    # X_base is returned as (Samples x Features)
    X_base, _, _, _ = align_labels_to_data(df_aligned, labels_dict, 'study_id')
    
    meta_df = pd.DataFrame({
        c: align_labels_to_data(df_aligned, labels_dict, c)[1] 
        for c in categories
    })

    results_summary = []

    # 2. METRIC CALCULATION
    for cat in categories:
        print(f"\n[Metrics: {cat.upper()}]")
        text_labels_np = np.array(meta_df[cat].tolist())
        valid_mask = ~np.isin(text_labels_np, ['unknown', 'unspecified', 'None', 'nan'])

        X_metric = X_base[valid_mask]
        text_labels_metric = text_labels_np[valid_mask]
        batch_text_labels_metric = meta_df['study_id'].values[valid_mask]

        unique_classes, num_labels_metric = np.unique(text_labels_metric, return_inverse=True)

        if X_metric.shape[0] < 5 or len(unique_classes) < 2:
            print(f"  Not enough valid samples/classes for {cat}.")
            sil_score, ari_score, knn_purity, var_explained, batch_asw = [np.nan] * 5
        else:
            # ---> KEY CHANGE: Skip PCA if we already have BulkFormer Embeddings <---
            if use_bulkformer:
                X_rep_metric = X_metric
            else:
                X_rep_metric, _ = run_pca(X_metric, n_components=min(50, X_metric.shape[0]-1))

            sil_score = silhouette_score(X_rep_metric, num_labels_metric, sample_size=min(5000, X_rep_metric.shape[0]))
            
            kmeans = MiniBatchKMeans(n_clusters=len(unique_classes), random_state=42).fit(X_rep_metric)
            ari_score = adjusted_rand_score(num_labels_metric, kmeans.labels_)

            knn = KNeighborsClassifier(n_neighbors=min(5, X_rep_metric.shape[0] - 1))
            knn_purity = cross_val_score(knn, X_rep_metric, num_labels_metric, cv=2).mean()

            var_explained = variance_explained_by_label(X_rep_metric, text_labels_metric)
            batch_asw = calculate_asw_batch_within_biology(X_rep_metric, batch_text_labels_metric, text_labels_metric)
            
            print(f"  Silhouette: {sil_score:.3f}, ARI: {ari_score:.3f}, KNN Purity: {knn_purity:.3f}, Var Exp: {var_explained:.3f}, Batch ASW: {batch_asw:.3f}")

        results_summary.append({
            'Category': cat, 'Silhouette': sil_score, 'ARI': ari_score, 
            'KNN_Purity': knn_purity, 'Variance_Explained': var_explained, 
            'Batch_ASW_within_Bio': batch_asw
        })

    # 3. VISUALIZATIONS
    print(f"\nGenerating standard UMAP & TSNE for {experiment_name}...")
    
    # ---> KEY CHANGE: Skip global PCA if using BulkFormer <---
    if use_bulkformer:
        X_rep_full = X_base 
    else:
        X_rep_full, _ = run_pca(X_base, n_components=min(50, X_base.shape[0]-1))
        
    for method, run_func in [("UMAP", run_umap), ("TSNE", run_tsne)]:
        emb = run_func(X_rep_full)
        output_path = f'{output_folder}/{experiment_name}_{method}.html'
        plot_interactive_projection(
            emb, 
            meta_df, 
            f'{experiment_name} - {method}', 
            output_path
        )

    res_df = pd.DataFrame(results_summary)
    res_df.to_csv(f'{output_folder}/{experiment_name}_metrics.csv', index=False)
    
    return res_df


if __name__ == "__main__":
    all_metrics = {}
    
    print("Loading Labels Map...")
    labels_map = make_df_from_labels(load_labels_study(LABELS_PATH), LABELS).to_dict() 
    
    # Provide the necessary dictionaries for the transformer
    # (These need to be loaded from your metadata files)
    mock_gene_lengths = pd.read_csv("metadata/gene_lengths.csv", index_col=0).to_dict()['length']
    mock_target_vocab = [line.strip() for line in open("metadata/bulkformer_gene_list.txt")]
    mock_ortholog_map = json.load(open("metadata/arabidopsis_to_human_hgnc.json"))
    
    stages = ['filter', 'imputed', 'study_corrected', 'rankin']
    
    for file in stages:
        data_path = f'{STORAGE_DIR}/final_data/{file}.csv'
        
        if os.path.exists(data_path):
            print(f"\n{'='*50}\nProcessing {file}\n{'='*50}")
            df = pd.read_csv(data_path, index_col=0)
            
            print("  Cleaning sample IDs...")
            df.columns = [c.split('.')[0].upper() for c in df.columns]

            print("  Backfilling missing study_ids using get_study()...")
            if 'study_id' not in labels_map:
                labels_map['study_id'] = {}
                
            count_filled = 0
            for sample in df.columns:
                if sample not in labels_map['study_id']:
                    study_val = get_study(sample)
                    labels_map['study_id'][sample.upper()] = study_val
                    count_filled += 1
            print(f"  -> Added study_id labels for {count_filled} samples.")

            output_dir = f"{CLUSTER_EXPLORATION_FIGURES_DIR}/interactive_plots/{file}"
            
            # --- TO USE BULKFORMER, SET use_bulkformer=True HERE ---
            metrics_df = run_exploration_on_dataframe(
                data_df=df,
                labels_dict=labels_map,
                experiment_name=file,
                output_folder=output_dir,
                use_bulkformer=True,
                gene_length_dict=mock_gene_lengths,
                target_vocab=mock_target_vocab,
                ortholog_map=mock_ortholog_map
            )
            
            all_metrics[file] = metrics_df
            
        else:
            print(f"Error: Data file not found at {data_path}")
            
    if len(all_metrics) > 1:
        comparison_output_dir = f"{CLUSTER_EXPLORATION_FIGURES_DIR}/interactive_plots/Comparisons"
        plot_metrics_comparison(
            metrics_dict=all_metrics, 
            metadata_df=pd.DataFrame(labels_map),
            output_folder=comparison_output_dir
        )