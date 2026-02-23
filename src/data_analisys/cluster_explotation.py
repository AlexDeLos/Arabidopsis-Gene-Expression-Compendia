import pandas as pd
import os
import sys
import gc
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics import silhouette_score, adjusted_rand_score
from sklearn.cluster import KMeans, MiniBatchKMeans
from sklearn.neighbors import KNeighborsClassifier
from sklearn.preprocessing import MultiLabelBinarizer
from sklearn.model_selection import cross_val_score
from scipy.stats import chi2_contingency

module_dir = './'
sys.path.append(module_dir)

from src.constants import *
from src.data_analisys.utils.cluster_exploration_utils_2 import (
    prepare_data_structure, align_labels_to_data, 
    run_pca, run_umap, run_tsne, plot_projection
)
# Assuming you have a function to load labels, e.g., load_labels_study
import numpy as np
from sklearn.metrics import silhouette_score

def calculate_asw_batch_within_biology(X_pca, batch_labels, bio_labels):
    """
    Calculates the batch effect severity while controlling for biological confounding.
    Lower score = better batch mixing (closer to 0).
    Higher score = strong batch effect remaining.
    """
    scores = []
    bio_labels = np.array(bio_labels)
    batch_labels = np.array(batch_labels)
    
    # 1. Iterate through each unique biological label (e.g., 'leaf', 'root')
    for bio_class in np.unique(bio_labels):
        if bio_class == 'unspecified':
            continue # Skip unspecified as it's not a reliable ground truth
            
        mask = (bio_labels == bio_class)
        X_sub = X_pca[mask]
        batch_sub = batch_labels[mask]
        
        # 2. We can only evaluate batch mixing if this tissue exists in multiple studies
        if len(X_sub) > 2 and len(np.unique(batch_sub)) > 1:
            # 3. Calculate how strongly the batches cluster *within* this tissue
            score = silhouette_score(X_sub, batch_sub)
            scores.append(abs(score)) # Use absolute value so 0 is strictly the target
            
    if not scores:
        return np.nan
        
    # Return the average batch effect across all shared biological classes
    return np.mean(scores)


from sklearn.linear_model import LinearRegression

def variance_explained_by_label(X_pca, labels):
    """
    Approximates the percentage of PCA variance explained by a specific label.
    X_pca: The PCA transformed data (e.g., top 50 PCs)
    labels: The categorical labels (study_id or tissue)
    """
    # 1. Convert categorical labels to one-hot encoded dummy variables
    # This allows the linear model to understand the categories
    labels_encoded = pd.get_dummies(labels).values
    
    # 2. Fit a multivariate linear regression
    # Predicting the 50 PCs using only the label matrix
    model = LinearRegression()
    model.fit(labels_encoded, X_pca)
    
    # 3. Calculate the R-squared score (Variance Explained)
    # A score of 0.40 means this label explains 40% of the variation in the data
    variance_explained = model.score(labels_encoded, X_pca)
    
    return variance_explained

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
            # A. Align Data for the current category
            X, text_labels, num_labels, n_classes = align_labels_to_data(df_aligned, labels_dict, cat)
            
            # Extract batch (study_id) labels for these exact same samples to calculate ASW
            # (We use '_' to discard the X, num_labels, n_classes since we only need the text array)
            _, batch_text_labels, _, _ = align_labels_to_data(df_aligned, labels_dict, 'study_id')
            
            n_unspecified = text_labels.count('unspecified')
            print(f"  -> Samples: {X.shape[0]} | Unlabeled: {n_unspecified} | Classes: {n_classes}")

            if X.shape[0] < 5:
                print("  -> Skipping: Not enough samples.")
                continue

            # B. PCA (Pre-processing)
            X_pca, var_ratio = run_pca(X, n_components=50)
            
            # C. Metrics
            print(f"  -> Calculating Metrics...")
            
            # 1. Standard Silhouette
            try:
                sil_score = silhouette_score(X_pca, num_labels, metric='euclidean', sample_size=5000)
            except ValueError:
                sil_score = -1 
            
            # 2. ARI (Global clustering)
            kmeans = MiniBatchKMeans(n_clusters=n_classes, batch_size=256, random_state=42).fit(X_pca)
            ari_score = adjusted_rand_score(num_labels, kmeans.labels_)
            
            # 3. KNN Purity
            knn = KNeighborsClassifier(n_neighbors=5)
            knn_purity = cross_val_score(knn, X_pca, num_labels, cv=5).mean()
            
            # 4. PVCA (Variance Explained)
            var_explained = variance_explained_by_label(X_pca, text_labels)
            
            # 5. Batch ASW within Biology
            batch_asw = calculate_asw_batch_within_biology(X_pca, batch_text_labels, text_labels)
            
            print(f"     * Silhouette: {sil_score:.4f}")
            print(f"     * ARI: {ari_score:.4f}")
            print(f"     * KNN Purity: {knn_purity:.4f}")
            print(f"     * Variance Explained: {var_explained:.4f}")
            if cat != 'study_id':
                print(f"     * Batch ASW (within {cat}): {batch_asw:.4f}")
            
            results_summary.append({
                'Category': cat,
                'Silhouette': sil_score,
                'ARI': ari_score,
                'KNN_Purity': knn_purity,
                'Variance_Explained': var_explained,
                'Batch_ASW_within_Bio': batch_asw,
                'Num_Classes': n_classes,
                'Unspecified_Count': n_unspecified
            })
            
            # D. Visualizations 
            
            # UMAP
            umap_emb = run_umap(X_pca)
            plot_projection(
                umap_emb, text_labels, 
                title=f'UMAP - {cat.capitalize()} (KNN Purity: {knn_purity:.2f})', 
                output_path=f'{output_folder}/{experiment_name}_{cat}_UMAP.svg'
            )
            
            # t-SNE
            tsne_emb = run_tsne(X_pca)
            plot_projection(
                tsne_emb, text_labels, 
                title=f't-SNE - {cat.capitalize()} (Var Explained: {var_explained:.2f})', 
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



def calculate_cramers_v(series1: pd.Series, series2: pd.Series) -> float:
    """Calculates Cramer's V statistic for association between two categorical series."""
    # Create contingency table
    contingency_table = pd.crosstab(series1, series2)
    # If table is empty or too small, return 0 correlation
    if contingency_table.empty or contingency_table.shape[0] < 2 or contingency_table.shape[1] < 2:
        return 0.0
        
    # Chi-square test
    chi2, _, _, _ = chi2_contingency(contingency_table)
    
    # Cramer's V calculation
    n = contingency_table.sum().sum()
    min_dim = min(contingency_table.shape) - 1
    # Avoid division by zero if min_dim is 0
    if min_dim == 0 or n == 0:
        return 0.0
        
    v = np.sqrt(chi2 / (n * min_dim))
    return v

def calculate_multilabel_association(study_series: pd.Series, multilabel_series: pd.Series) -> float:
    """
    Calculates the mean Cramer's V association between a single categorical variable (study_id)
    and a multi-label variable (e.g., treatment lists).
    """
    # 1. Ensure everything is a list (convert strings/NaNs safely)
    clean_labels = multilabel_series.apply(
        lambda x: list(x) if isinstance(x, (list, tuple, set)) else ([str(x)] if pd.notna(x) else [])
    ).tolist()
    
    # 2. Binarize the multi-labels into a matrix (0s and 1s)
    mlb = MultiLabelBinarizer()
    binary_matrix = mlb.fit_transform(clean_labels)
    
    # If there are no labels, return 0
    if binary_matrix.shape[1] == 0:
        return 0.0
        
    v_scores = []
    
    # 3. Calculate association for each individual label vs the study_id
    for i, label in enumerate(mlb.classes_):
        binary_col = binary_matrix[:, i]
        
        # Skip labels that are present in EVERY sample or NO samples (zero variance)
        if len(np.unique(binary_col)) < 2:
            continue
            
        # Calculate Cramer's V between the study and this specific label's presence/absence
        v = calculate_cramers_v(study_series, binary_col)
        v_scores.append(v)
        
    # 4. Return the average association across all labels in this category
    if not v_scores:
        return 0.0
        
    return float(np.mean(v_scores))

# Updated Signature to accept metadata_df
def plot_metrics_comparison(metrics_dict: dict, 
                            metadata_df: pd.DataFrame, 
                            output_folder: str, 
                            experiment_name: str = "Normalization_Comparison"):
    """
    Plots a comparison of pre- and post-normalization metrics in a 2x3 grid.
    The 6th plot shows inherent confounding between study_id and biological variables.
    
    :param metrics_dict: A dictionary of DataFrames, e.g., {'Raw': filter_df, 'ComBat': combat_df}
    :param metadata_df: The original metadata DataFrame containing actual labels (tissue, treatment, study_id)
    :param output_folder: Where to save the plots.
    """
    print(f"\n[Generating Metric Comparison Plots...]")
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # --- 1. Prepare Metrics Data ---
    combined_data = []
    for stage_name, df in metrics_dict.items():
        df_copy = df.copy()
        df_copy['Stage'] = stage_name
        combined_data.append(df_copy)
    
    plot_df = pd.concat(combined_data, ignore_index=True)

    # --- 2. Prepare Confounding Correlation Data (New Step) ---
    bio_targets = LABELS
    confounding_data = []
    
    if 'study_id' in metadata_df.columns:
        for target in bio_targets:
            if target in metadata_df.columns:
                target_data = metadata_df[target]
                
                # Check if this column contains lists/sets
                has_lists = target_data.apply(lambda x: isinstance(x, (list, tuple, set))).any()
                
                if has_lists:
                    # Use the new multi-label metric
                    v_score = calculate_multilabel_association(metadata_df['study_id'], target_data)
                elif target_data.nunique() > 1:
                    # Use the standard metric for flat strings
                    v_score = calculate_cramers_v(metadata_df['study_id'], target_data)
                else:
                    v_score = 0.0
                    
                confounding_data.append({'Variable': target.capitalize(), 'Cramers_V': v_score})
            else:
                confounding_data.append({'Variable': target.capitalize() + " (Missing)", 'Cramers_V': 0})
                
    confounding_df = pd.DataFrame(confounding_data)

    # Set presentation-ready aesthetic
    sns.set_theme(style="whitegrid", context="talk")

    # 3. Setup the figure grid (2x3 plots)
    fig, axes = plt.subplots(2, 3, figsize=(24, 12))
    fig.suptitle('Batch Correction Evaluation & Confounding Check', fontsize=20, fontweight='bold', y=0.99)

    # --- Plot A: Variance Explained ---
    sns.barplot(data=plot_df, x='Category', y='Variance_Explained', hue='Stage', ax=axes[0, 0], palette='Set2')
    axes[0, 0].set_title('A. Variance Explained (Higher = More Influence)')
    axes[0, 0].set_ylabel('R² Score')
    if not plot_df['Variance_Explained'].dropna().empty:
        axes[0, 0].set_ylim(0, max(plot_df['Variance_Explained'].dropna()) * 1.2)

    # --- Plot B: KNN Purity ---
    sns.barplot(data=plot_df, x='Category', y='KNN_Purity', hue='Stage', ax=axes[0, 1], palette='Set2')
    axes[0, 1].set_title('B. KNN Purity (Higher = Better Local Grouping)')
    axes[0, 1].set_ylabel('Purity Score')
    axes[0, 1].set_ylim(0, 1.1)

    # --- Plot C: Batch ASW within Biology ---
    bio_only_df = plot_df[plot_df['Category'] != 'study_id']
    sns.barplot(data=bio_only_df, x='Category', y='Batch_ASW_within_Bio', hue='Stage', ax=axes[0, 2], palette='Set2')
    axes[0, 2].set_title('C. Batch ASW within Bio (Lower = Better Mixing)')
    axes[0, 2].set_ylabel('Silhouette Score of Batch')
    if not bio_only_df['Batch_ASW_within_Bio'].dropna().empty:
         # Ensure y-limit accommodates 0 if some scores are negative
        ymax = max(bio_only_df['Batch_ASW_within_Bio'].dropna()) * 1.2
        ymin = min(bio_only_df['Batch_ASW_within_Bio'].dropna().min(), 0) - 0.05
        axes[0, 2].set_ylim(ymin, ymax)

    # --- Plot D: Adjusted Rand Index (ARI) ---
    sns.barplot(data=plot_df, x='Category', y='ARI', hue='Stage', ax=axes[1, 0], palette='Set2')
    axes[1, 0].set_title('D. Adjusted Rand Index (ARI) (Lower = Better Mixed)')
    axes[1, 0].set_ylabel('ARI Score')

    # --- Plot E: Silhouette Score ---
    sns.barplot(data=plot_df, x='Category', y='Silhouette', hue='Stage', ax=axes[1, 1], palette='Set2')
    axes[1, 1].set_title('E. Silhouette Score (Higher = Tighter Clusters)')
    axes[1, 1].set_ylabel('Silhouette Score')

    # --- Plot F: Experimental Confounding (NEW) ---
    # Plotting Cramer's V correlation between study_id and other variables
    ax_conf = axes[1, 2]
    if not confounding_df.empty:
        # Use a distinct color, as this isn't comparing stages
        sns.barplot(data=confounding_df, x='Variable', y='Cramers_V', ax=ax_conf, color=sns.color_palette("deep")[0])
        ax_conf.set_title('F. Inherent Confounding: Study ID vs. Biology')
        ax_conf.set_ylabel("Association (Cramer's V)\n(1.0 = Perfect Confounding)")
        ax_conf.set_ylim(0, 1.1)
        
        # Add data labels on top of bars for clarity
        for p in ax_conf.patches:
             ax_conf.annotate(f'{p.get_height():.2f}', 
                              (p.get_x() + p.get_width() / 2., p.get_height()), 
                              ha = 'center', va = 'center', 
                              xytext = (0, 9), 
                              textcoords = 'offset points',
                              fontsize=12)
    else:
        ax_conf.text(0.5, 0.5, "Metadata missing or insufficient data\nfor confounding check.", 
                     ha='center', va='center')
        ax_conf.set_title('F. Inherent Confounding Check')


    # Formatting adjustments
    active_axes = axes.flatten()
    for i, ax in enumerate(active_axes):
        ax.set_xlabel('')
        # Handle legends: keep only on plot A, remove from B,C,D,E. Plot F doesn't have one.
        if ax.get_legend() is not None:
            if i == 0: 
                ax.legend(title='Pipeline Stage', loc='upper right', bbox_to_anchor=(1.0, 1.05))
            else:
                ax.get_legend().remove()

    plt.tight_layout(rect=[0, 0.03, 1, 0.96])
    
    # Save the figure
    output_path = os.path.join(output_folder, f"{experiment_name}_Summary_with_Confounding.svg")
    # Use try/except blocks for saving in case of tight layout issues on some systems
    try:
        plt.savefig(output_path, format='svg', bbox_inches='tight')
        plt.savefig(output_path.replace('.svg', '.png'), format='png', bbox_inches='tight', dpi=300)
    except Exception as e:
         print(f"[Warning] Could not save with tight bbox layout: {e}. Saving standard layout.")
         plt.savefig(output_path, format='svg')
         plt.savefig(output_path.replace('.svg', '.png'), format='png', dpi=300)

    print(f"  -> Saved comparison plots to {output_path.replace('.svg', '.png')}")
    plt.close()
        
from src.data_analisys.utils.cluster_exploration_utils import *
# --- Example Usage / Main Block ---
if __name__ == "__main__":
    all_metrics = {}
    print(f"Loading labels from: {LABELS_PATH}")
    labels = load_labels_study(LABELS_PATH)
    # labels = keys_upper(labels)

    # Prepare Label Dictionary
    labels_types = LABELS
    labels_df = make_df_from_labels(labels, labels_types)
    labels_map = labels_df.to_dict() # Structure: {'category': {'sample': 'value'}}
    
    del labels,labels_df
    for file in ['filter','imputed','study_corrected']:
        # 1. Load Labels

        # 2. Load Data
        data_path = f'{STORAGE_DIR}/final_data/{file}.csv' # Overriding constant as per your snippet
        
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
                    labels_map['study_id'][sample.upper()] = study_val
                    count_filled += 1
            df.columns = [c.upper() for c in df.columns]
            print(f"  -> Added study_id labels for {count_filled} samples.")

            # 3. Run Pipeline
            output_dir = f"{CLUSTER_EXPLORATION_FIGURES_DIR}/refactored_plots/{file}"
            
            print(f"Starting exploration pipeline for {EXPERIMENT_NAME}...")
            
            # Note: Ensure 'study_id' is added to the categories list inside 
            # run_exploration_on_dataframe if you want it plotted!
            metrics_df = run_exploration_on_dataframe(
                data_df=df,
                labels_dict=labels_map,
                experiment_name=file,
                output_folder=output_dir
            )
            
            # Store it for comparison plotting
            # We rename the keys to be more presentable ('filter' -> 'Raw Data')
            # stage_name = 'Raw Data' if file == 'filter' else 'ComBat Corrected'
            all_metrics[file] = metrics_df
            
        else:
            print(f"Error: Data file not found at {data_path}")
    if len(all_metrics) > 1:
        comparison_output_dir = f"{CLUSTER_EXPLORATION_FIGURES_DIR}/refactored_plots/Comparisons"
        plot_metrics_comparison(
            metrics_dict=all_metrics, 
            metadata_df=pd.DataFrame(labels_map),
            output_folder=comparison_output_dir
        )