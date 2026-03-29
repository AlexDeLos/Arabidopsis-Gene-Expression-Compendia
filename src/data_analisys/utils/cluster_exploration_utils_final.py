import os
import sys
import json
import glob
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from scipy.stats import chi2_contingency
from sklearn.preprocessing import MultiLabelBinarizer

import umap
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.metrics import silhouette_score
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import LabelEncoder


# =============================================================================
# 1. DATA AND LABEL PREPARATION (Updated for TULIP LLM Format)
# =============================================================================

def load_labels_study(labels_dir: str) -> dict:
    """
    Loads the TULIP LLM JSON label files.
    Reads all {GSE_ID}.json files from labels_dir.
    
    The new JSON format is:
    { "GSM_ID": { "axis": ["val1"], "axis2": [{"val": "Drought", "intensity": 2}] } }
    
    Returns a transposed dictionary ready for alignment, including sub-attributes: 
    { 
       'treatment': { GSM_ID: 'Drought' },
       'treatment_intensity': { GSM_ID: '2' }
    }
    """
    aggregated_data = {}
    
    # Support both a single file or a directory of files
    if os.path.isfile(labels_dir):
        files = [labels_dir]
    else:
        files = glob.glob(os.path.join(labels_dir, "*.json"))

    print(f"Loading labels from {len(files)} JSON files in {labels_dir}...")
    for file in files:
        with open(file, 'r') as f:
            try:
                data = json.load(f)
                aggregated_data.update(data)
            except json.JSONDecodeError:
                print(f"  ! Warning: Could not parse {file}")

    axis_map = {}
    for gsm_id, axes_dict in aggregated_data.items():
        if not isinstance(axes_dict, dict):
            continue
            
        for axis, values in axes_dict.items():
            if axis not in axis_map:
                axis_map[axis] = {}

            if isinstance(values, list):
                if len(values) == 0:
                    axis_map[axis][gsm_id.upper()] = "unspecified"
                
                elif isinstance(values[0], dict):
                    # Contains sub-attributes (e.g., 'val' and 'intensity')
                    vals = []
                    sub_attrs = {} 
                    
                    for v_dict in values:
                        # 1. Grab canonical value
                        canonical = str(v_dict.get('val', 'unspecified'))
                        vals.append(canonical)
                        
                        # 2. Grab any other sub-keys (like 'intensity')
                        for k, v in v_dict.items():
                            if k == 'val': continue
                            if k not in sub_attrs: sub_attrs[k] = []
                            sub_attrs[k].append(str(v))
                            
                    # Assign the main axis (e.g. 'treatment' = 'Chemical + Heat')
                    val_str = " + ".join(vals)
                    if val_str.lower() in ["none", "", "nan", "unknown"]: 
                        val_str = "unspecified"
                    axis_map[axis][gsm_id.upper()] = val_str
                    
                    # Assign the sub-attribute axes (e.g. 'treatment_intensity' = '2 + 1')
                    for sub_k, sub_list in sub_attrs.items():
                        sub_axis = f"{axis}_{sub_k}"
                        if sub_axis not in axis_map:
                            axis_map[sub_axis] = {}
                        axis_map[sub_axis][gsm_id.upper()] = " + ".join(sub_list)

                else:
                    # Standard flat list of strings
                    val_str = " + ".join([str(v) for v in values])
                    if val_str.lower() in ["none", "", "nan", "unknown"]: 
                        val_str = "unspecified"
                    axis_map[axis][gsm_id.upper()] = val_str
                    
            elif isinstance(values, str):
                val_str = values if values.lower() not in ["none", "", "nan", "unknown"] else "unspecified"
                axis_map[axis][gsm_id.upper()] = val_str
            else:
                axis_map[axis][gsm_id.upper()] = str(values)

    return axis_map

def make_df_from_labels(labels_dict: dict) -> pd.DataFrame:
    """
    Converts the parsed {axis: {sample: label}} dictionary into a pandas DataFrame.
    Samples become the index, and Label Axes become the columns.
    """
    df = pd.DataFrame(labels_dict)
    df.index.name = "Sample_ID"
    # Fill missing overlaps with 'unspecified'
    df = df.fillna("unspecified")
    return df

def get_study(sample_id: str) -> str:
    """
    Extracts the study ID (GSE) from a given sample ID (GSM) if possible.
    (This usually relies on an external mapping table/database. Add your specific
    SQL/CSV lookup logic here if needed, otherwise it returns a placeholder).
    """
    # Replace with specific mapping logic if you have a GSM -> GSE dictionary
    
    return "Unknown_Study"

def prepare_data_structure(df: pd.DataFrame) -> pd.DataFrame:
    """
    Ensures data is (Samples x Genes). 
    """
    if df.shape[0] > df.shape[1] and df.shape[0] > 10000:
        print(f"  -> Transposing dataframe from {df.shape} to (Samples, Genes)...")
        return df.T
    return df

def align_labels_to_data(df: pd.DataFrame, labels_dict: dict, label_category: str) -> list:
    """
    Aligns the dataframe samples with the labels dictionary.
    Includes logic to fallback to Upper Case keys and handle missing samples.
    """
    if label_category in labels_dict:
        sample_to_label_map = labels_dict[label_category]
    elif label_category.upper() in labels_dict:
        sample_to_label_map = labels_dict[label_category.upper()]
    else:
        print(f"    ! Warning: Category '{label_category}' not found. All samples will be 'unspecified'.")
        sample_to_label_map = {}

    cleaned_labels = []
    for s in df.index:
        s_upper = str(s).upper()
        # Look for upper case first, then exact match, fallback to 'unspecified'
        label = sample_to_label_map.get(s_upper, sample_to_label_map.get(s, 'unspecified'))
        cleaned_labels.append(label)
        
    return cleaned_labels


# =============================================================================
# 2. DIMENSIONALITY REDUCTION
# =============================================================================

def run_pca(df: pd.DataFrame, n_components=50):
    print(f"  Running PCA (n_components={n_components})...")
    pca = PCA(n_components=min(n_components, df.shape[0], df.shape[1]))
    embedding = pca.fit_transform(df)
    return embedding, pca

def run_umap(pca_embedding, n_components=2):
    print(f"  Running UMAP (n_components={n_components})...")
    reducer = umap.UMAP(n_components=n_components, random_state=42)
    embedding = reducer.fit_transform(pca_embedding)
    return embedding

def run_tsne(pca_embedding, n_components=2):
    print(f"  Running t-SNE (n_components={n_components})...")
    # Dynamically adjust perplexity based on sample size
    perplexity = min(30, max(5, pca_embedding.shape[0] // 3))
    tsne = TSNE(n_components=n_components, random_state=42, perplexity=perplexity)
    embedding = tsne.fit_transform(pca_embedding)
    return embedding


# =============================================================================
# 3. CLUSTER & METRIC EVALUATION
# =============================================================================

def calculate_asw_batch_within_biology(X_pca, batch_labels, bio_labels) -> float:
    """
    Calculates the Average Silhouette Width for Batch within Biological groups.
    Measures how well-mixed batches are within the same biological label.
    """
    scores = []
    bio_labels = np.array(bio_labels)
    batch_labels = np.array(batch_labels)
    
    for bio_class in np.unique(bio_labels):
        if bio_class in ['unspecified', 'unknown', 'None', 'nan']:
            continue 
            
        mask = (bio_labels == bio_class)
        X_sub = X_pca[mask]
        batch_sub = batch_labels[mask]
        
        # Need at least 2 batches and 2 samples to compute silhouette
        if len(X_sub) > 2 and len(np.unique(batch_sub)) > 1:
            try:
                score = silhouette_score(X_sub, batch_sub)
                scores.append(score)
            except ValueError:
                continue

    if len(scores) > 0:
        return np.mean(scores)
    return 0.0

def variance_explained_by_label(data, labels) -> float:
    """
    Calculates the variance in the dataset explained by the provided labels.
    Approximated using a Linear Regression R^2 score across principal components.
    """
    valid_mask = ~pd.Series(labels).isin(['unspecified', 'unknown', 'nan', 'None'])
    if valid_mask.sum() < 2:
        return 0.0
        
    data_sub = data[valid_mask]
    labels_sub = np.array(labels)[valid_mask]
    
    le = LabelEncoder()
    y_enc = le.fit_transform(labels_sub)
    
    if len(np.unique(y_enc)) < 2:
        return 0.0
        
    try:
        model = LinearRegression()
        model.fit(data_sub, y_enc)
        return model.score(data_sub, y_enc)
    except Exception:
        return 0.0

def calculate_cramers_v(series1: pd.Series, series2: pd.Series) -> float:
    contingency_table = pd.crosstab(series1, series2)
    if contingency_table.empty or contingency_table.shape[0] < 2 or contingency_table.shape[1] < 2:
        return 0.0
        
    chi2, _, _, _ = chi2_contingency(contingency_table)
    n = contingency_table.sum().sum()
    min_dim = min(contingency_table.shape) - 1
    
    if min_dim == 0 or n == 0:
        return 0.0
    return np.sqrt(chi2 / (n * min_dim))

def calculate_multilabel_association(study_series: pd.Series, multilabel_series: pd.Series) -> float:
    clean_labels = multilabel_series.apply(
        lambda x: list(x) if isinstance(x, (list, tuple, set)) else ([str(x)] if pd.notna(x) else [])
    ).tolist()
    
    mlb = MultiLabelBinarizer()
    binary_matrix = mlb.fit_transform(clean_labels)
    
    if binary_matrix.shape[1] == 0:
        return 0.0
        
    v_scores = []
    for i, label in enumerate(mlb.classes_):
        binary_col = binary_matrix[:, i]
        if len(np.unique(binary_col)) < 2:
            continue
        v = calculate_cramers_v(study_series, binary_col)
        v_scores.append(v)
        
    if not v_scores:
        return 0.0
    return float(np.mean(v_scores))

def plot_metrics_comparison(metrics_dict: dict, 
                            metadata_df: pd.DataFrame, 
                            output_folder: str,
                            bio_targets,
                            experiment_name: str = "Normalization_Comparison"):
    
    print(f"\n[Generating Metric Comparison Plots...]")
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    combined_data = []
    for stage_name, df in metrics_dict.items():
        df_copy = df.copy()
        df_copy['Stage'] = stage_name
        combined_data.append(df_copy)
    
    plot_df = pd.concat(combined_data, ignore_index=True)

    confounding_data = []
    
    if 'study_id' in metadata_df.columns:
        for target in bio_targets:
            if target in metadata_df.columns:
                target_data = metadata_df[target]
                has_lists = target_data.apply(lambda x: isinstance(x, (list, tuple, set))).any()
                
                if has_lists:
                    v_score = calculate_multilabel_association(metadata_df['study_id'], target_data)
                elif target_data.nunique() > 1:
                    v_score = calculate_cramers_v(metadata_df['study_id'], target_data)
                else:
                    v_score = 0.0
                    
                confounding_data.append({'Variable': target.capitalize(), 'Cramers_V': v_score})
            else:
                confounding_data.append({'Variable': target.capitalize() + " (Missing)", 'Cramers_V': 0})
                
    confounding_df = pd.DataFrame(confounding_data)

    sns.set_theme(style="whitegrid", context="talk")
    fig, axes = plt.subplots(2, 3, figsize=(24, 12))
    fig.suptitle('Batch Correction Evaluation & Confounding Check\n(Calculated exclusively on valid known labels)', fontsize=20, fontweight='bold', y=0.99)

    sns.barplot(data=plot_df, x='Category', y='Variance_Explained', hue='Stage', ax=axes[0, 0], palette='Set2')
    axes[0, 0].set_title('A. Variance Explained (Higher = More Influence)')
    axes[0, 0].set_ylabel('R² Score')

    sns.barplot(data=plot_df, x='Category', y='KNN_Purity', hue='Stage', ax=axes[0, 1], palette='Set2')
    axes[0, 1].set_title('B. KNN Purity (Higher = Better Local Grouping)')
    axes[0, 1].set_ylabel('Purity Score')
    axes[0, 1].set_ylim(0, 1.1)

    bio_only_df = plot_df[plot_df['Category'] != 'study_id']
    sns.barplot(data=bio_only_df, x='Category', y='Batch_ASW_within_Bio', hue='Stage', ax=axes[0, 2], palette='Set2')
    axes[0, 2].set_title('C. Batch ASW within Bio (Lower = +Study Mixing)')
    axes[0, 2].set_ylabel('Silhouette Score of Batch')

    sns.barplot(data=plot_df, x='Category', y='ARI', hue='Stage', ax=axes[1, 0], palette='Set2')
    axes[1, 0].set_title('D. Adjusted Rand Index (Align. w. clutsers)')
    axes[1, 0].set_ylabel('ARI Score')

    sns.barplot(data=plot_df, x='Category', y='Silhouette', hue='Stage', ax=axes[1, 1], palette='Set2')
    axes[1, 1].set_title('E. Silhouette Score (Higher = Tighter Clusters)')
    axes[1, 1].set_ylabel('Silhouette Score')

    ax_conf = axes[1, 2]
    if not confounding_df.empty:
        sns.barplot(data=confounding_df, x='Variable', y='Cramers_V', ax=ax_conf, color=sns.color_palette("deep")[0])
        ax_conf.set_title('F. Inherent Confounding: Study ID vs. Biology')
        ax_conf.set_ylabel("Association (Cramer's V)\n(1.0 = Perfect Confounding)")
        ax_conf.set_ylim(0, 1.1)
        
        for p in ax_conf.patches:
             ax_conf.annotate(f'{p.get_height():.2f}', 
                              (p.get_x() + p.get_width() / 2., p.get_height()), 
                              ha = 'center', va = 'center', 
                              xytext = (0, 9), textcoords = 'offset points', fontsize=12)
    else:
        ax_conf.text(0.5, 0.5, "Metadata missing or insufficient data\nfor confounding check.", ha='center', va='center')
        ax_conf.set_title('F. Inherent Confounding Check')

    active_axes = axes.flatten()
    for i, ax in enumerate(active_axes):
        ax.set_xlabel('')
        if ax.get_legend() is not None:
            if i == 0: 
                ax.legend(title='Pipeline Stage', loc='upper right', bbox_to_anchor=(1.0, 1.05))
            else:
                ax.get_legend().remove()

    plt.tight_layout(rect=[0, 0.03, 1, 0.96])
    
    output_path = os.path.join(output_folder, f"{experiment_name}_Summary_with_Confounding.svg")
    try:
        plt.savefig(output_path, format='svg', bbox_inches='tight')
        plt.savefig(output_path.replace('.svg', '.png'), format='png', bbox_inches='tight', dpi=300)
    except Exception as e:
         print(f"[Warning] Could not save with tight bbox layout: {e}. Saving standard layout.")
         plt.savefig(output_path, format='svg')
         plt.savefig(output_path.replace('.svg', '.png'), format='png', dpi=300)

    print(f"  -> Saved comparison plots to {output_path.replace('.svg', '.png')}")
    plt.close()

# def plot_metrics_comparison(metrics_dict: dict, metadata_df: pd.DataFrame, output_folder: str):
#     """
#     Aggregates and plots metrics across different stages/experiments.
#     Matches exact signature expected by cluster_explotation_new.py.
#     """
#     print(f"Saving comparison metrics to {output_folder}...")
#     os.makedirs(output_folder, exist_ok=True)
    
#     combined_rows = []
#     for exp_name, df in metrics_dict.items():
#         if isinstance(df, pd.DataFrame):
#             # Ensure the experiment name is tracked
#             df_copy = df.copy()
#             df_copy['Experiment'] = exp_name
#             combined_rows.append(df_copy)
            
#     if not combined_rows:
#         print("  ! No metrics DataFrames found to compare.")
#         return
        
#     combined_df = pd.concat(combined_rows, ignore_index=True)
#     combined_df.to_csv(f"{output_folder}/combined_metrics.csv", index=False)
    
#     # Generate Comparison Plot
#     if 'Metric' in combined_df.columns and 'Value' in combined_df.columns and 'Label_Axis' in combined_df.columns:
#         fig = px.bar(
#             combined_df, 
#             x='Label_Axis', 
#             y='Value', 
#             color='Experiment',
#             facet_col='Metric', 
#             barmode='group',
#             template='plotly_dark',
#             title="Cluster Metrics Comparison Across Pipeline Stages"
#         )
#         # Allow Y axes to adjust independently for different metric scales
#         fig.update_yaxes(matches=None) 
#         fig.write_html(f"{output_folder}/Metrics_Comparison_Barplot.html")
#         print("  -> Exported Metrics_Comparison_Barplot.html")