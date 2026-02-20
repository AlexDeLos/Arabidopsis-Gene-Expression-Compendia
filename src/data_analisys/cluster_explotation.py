import pandas as pd
import os
import sys
import gc
from sklearn.metrics import silhouette_score, adjusted_rand_score
from sklearn.cluster import KMeans, MiniBatchKMeans

module_dir = './'
sys.path.append(module_dir)

from src.constants import *
from src.data_analisys.utils.cluster_exploration_utils_2 import (
    prepare_data_structure, align_labels_to_data, 
    run_pca, run_umap, run_tsne, plot_projection
)
# Assuming you have a function to load labels, e.g., load_labels_study

def run_exploration_on_dataframe(
    data_df: pd.DataFrame, 
    labels_dict: dict, 
    experiment_name: str,
    output_folder: str
):
    
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # 1. Transpose if needed (Samples x Genes)
    df_aligned = prepare_data_structure(data_df)
    
    categories = ['treatment', 'tissue', 'medium','study_id']
    
    results_summary = []

    for cat in categories:
        print(f"\n[Processing Category: {cat.upper()}]")
        
        try:
            # A. Align Data (Includes 'unspecified')
            X, text_labels, num_labels, n_classes = align_labels_to_data(df_aligned, labels_dict, cat)
            
            # Count 'unspecified'
            n_unspecified = text_labels.count('unspecified')
            print(f"  -> Samples: {X.shape[0]} | Unlabeled: {n_unspecified} | Classes: {n_classes}")

            if X.shape[0] < 5:
                print("  -> Skipping: Not enough samples.")
                continue

            # B. PCA (Pre-processing)
            X_pca, var_ratio = run_pca(X, n_components=50)
            
            # C. Metrics
            # Note: Unspecified labels will form their own 'class' in these metrics.
            # This is technically 'incorrect' for ARI (since unspecified isn't a ground truth),
            # but valid for checking if unspecified samples cluster together.
            print(f"  -> Calculating Metrics...")
            
            try:
                sil_score = silhouette_score(X_pca, num_labels, metric='euclidean', sample_size=5000)
            except ValueError:
                sil_score = -1 
            
            # ARI: K-Means vs Labels
            kmeans = MiniBatchKMeans(n_clusters=n_classes, batch_size=256, random_state=42).fit(X_pca)
            ari_score = adjusted_rand_score(num_labels, kmeans.labels_)
            
            print(f"     * Silhouette: {sil_score:.4f}")
            print(f"     * ARI: {ari_score:.4f}")
            
            results_summary.append({
                'Category': cat,
                'Silhouette': sil_score,
                'ARI': ari_score,
                'Num_Classes': n_classes,
                'Unspecified_Count': n_unspecified
            })

            # D. Visualizations 
            
            # UMAP
            umap_emb = run_umap(X_pca)
            plot_projection(
                umap_emb, text_labels, 
                title=f'UMAP - {cat.capitalize()} (ARI: {ari_score:.2f})', 
                output_path=f'{output_folder}/{experiment_name}_{cat}_UMAP.svg'
            )
            
            # t-SNE
            tsne_emb = run_tsne(X_pca)
            plot_projection(
                tsne_emb, text_labels, 
                title=f't-SNE - {cat.capitalize()} (Sil: {sil_score:.2f})', 
                output_path=f'{output_folder}/{experiment_name}_{cat}_TSNE.svg'
            )

            del X, X_pca, umap_emb, tsne_emb
            gc.collect()

        except Exception as e:
            print(f"  -> Error processing {cat}: {e}")
            import traceback
            traceback.print_exc()

    res_df = pd.DataFrame(results_summary)
    res_df.to_csv(f'{output_folder}/{experiment_name}_metrics.csv', index=False)
    print("\nProcessing Complete.")
    return res_df



from src.data_analisys.utils.cluster_exploration_utils import *
# --- Example Usage / Main Block ---
if __name__ == "__main__":
    for file in ['filter','study_corrected']:
        # 1. Load Labels
        print(f"Loading labels from: {LABELS_PATH}")
        labels = load_labels_study(LABELS_PATH)
        # labels = keys_upper(labels)

        # Prepare Label Dictionary
        labels_types = LABELS
        labels_df = make_df_from_labels(labels, labels_types)
        labels_map = labels_df.to_dict() # Structure: {'category': {'sample': 'value'}}
        
        del labels, labels_df # Clean up

        # 2. Load Data
        data_path = f'./new_storage/final_data/{file}.csv' # Overriding constant as per your snippet
        
        if os.path.exists(data_path):
            print(f"Loading expression data from: {data_path}")
            df = pd.read_csv(data_path, index_col=0)
            
            # --- PRE-PROCESSING: Clean Sample IDs ---
            print("  Cleaning sample IDs (stripping .1, .2 suffixes)...")
            df.columns = [c.split('.')[0] for c in df.columns]

            # --- NEW STEP: Backfill Missing Study IDs ---
            print("  Backfilling missing study_ids using get_study()...")
            
            # Ensure the category dictionary exists
            if 'study_id' not in labels_map:
                labels_map['study_id'] = {}
                
            count_filled = 0
            for sample in df.columns:
                # Check if this specific sample is missing a study_id label
                if sample not in labels_map['study_id']:
                    # Generate the label
                    study_val = get_study(sample)
                    # Assign it to the map
                    labels_map['study_id'][sample] = study_val
                    count_filled += 1
                    
            print(f"  -> Added study_id labels for {count_filled} samples.")

            # 3. Run Pipeline
            output_dir = f"{CLUSTER_EXPLORATION_FIGURES_DIR}/refactored_plots/{file}"
            
            print(f"Starting exploration pipeline for {EXPERIMENT_NAME}...")
            
            # Note: Ensure 'study_id' is added to the categories list inside 
            # run_exploration_on_dataframe if you want it plotted!
            results = run_exploration_on_dataframe(
                data_df=df,
                labels_dict=labels_map,
                experiment_name=file,
                output_folder=output_dir
            )
            
            print("\nPipeline finished successfully.")
            print(f"Outputs saved to: {output_dir}")
            
        else:
            print(f"Error: Data file not found at {data_path}")