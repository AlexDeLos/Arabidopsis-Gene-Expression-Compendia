import os
import sys
import json
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from scipy.stats import chi2_contingency
from sklearn.preprocessing import MultiLabelBinarizer
from sklearn.model_selection import KFold

import umap
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.metrics import silhouette_score
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import LabelEncoder

import urllib.error
import time
import re
from Bio import Entrez
from src.constants import RNA_MA  # Pull in your new constant

Entrez.email = "alexdelossanto@tudelft.nl"

# --- GLOBAL CACHE --- 
# Prevents querying the exact same GSM across Tissue, Treatment, Medium, etc.
# --- PERSISTENT GLOBAL CACHE --- 
CACHE_FILE = "srr_mapping_cache.json"

def load_srr_cache() -> dict:
    """Loads the previously saved SRR mappings from disk."""
    if os.path.exists(CACHE_FILE):
        try:
            with open(CACHE_FILE, 'r') as f:
                return json.load(f)
        except Exception as e:
            print(f"    [!] Warning: Could not load SRR cache: {e}")
            return {}
    return {}

def save_srr_cache():
    """Saves the current SRR mappings to disk."""
    try:
        with open(CACHE_FILE, 'w') as f:
            json.dump(SRR_CACHE, f, indent=4)
    except Exception as e:
        print(f"    [!] Warning: Could not save SRR cache: {e}")

# Initialize the cache in memory when the module loads
SRR_CACHE = load_srr_cache()

def get_srr_ids(gsm_id: str, max_retries=5) -> list:
    gsm_id = gsm_id.strip().upper()
    
    # Return immediately if we already fetched this GSM ID earlier
    if gsm_id in SRR_CACHE:
        return SRR_CACHE[gsm_id]
        
    print(f"    [Cache Miss] Fetching SRR for {gsm_id} from NCBI...")
    
    for attempt in range(max_retries):
        try:
            # 1. Proactive Rate Limiting (ensures < 3 requests per second)
            time.sleep(0.4) 
            
            handle = Entrez.esearch(db="sra", term=gsm_id)
            record = Entrez.read(handle)
            handle.close()
            
            if not record['IdList']:
                SRR_CACHE[gsm_id] = []
                save_srr_cache() # Save the empty result so we don't try again next run
                return []
    
            handle = Entrez.esummary(db="sra", id=",".join(record['IdList']))
            summaries = Entrez.read(handle)
            handle.close()
            
            run_ids = []
            for summary in summaries:
                run_ids.extend(re.findall(r'acc="([A-Z0-9]+)"', summary.get('Runs', '')))
            
            res = list(set(run_ids))
            
            # Save to memory and immediately write to disk
            SRR_CACHE[gsm_id] = res 
            save_srr_cache()        
            
            return res
        
        except urllib.error.HTTPError as e:
            # 2. Reactive Backoff: If we still hit the limit, wait and retry
            if e.code == 429:
                wait_time = 2 ** attempt  # Waits 1s, 2s, 4s, 8s, 16s...
                print(f"    [!] HTTP 429 (Too Many Requests) for {gsm_id}. Waiting {wait_time}s before retry {attempt+1}/{max_retries}...")
                time.sleep(wait_time)
            else:
                print(f"    [!] HTTP Error {e.code} for {gsm_id}. Skipping.")
                SRR_CACHE[gsm_id] = []
                save_srr_cache()
                return []
        except Exception as e:
            print(f"    [!] Connection error for {gsm_id}: {e}. Retrying in 5s...")
            time.sleep(5)
            
    print(f"    [!] Failed to retrieve SRR for {gsm_id} after {max_retries} retries.")
    SRR_CACHE[gsm_id] = []
    save_srr_cache()
    return []
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
                
                # --- NEW LOGIC: Handle both Lists and Dictionaries ---
                if isinstance(data, list):
                    for item in data:
                        if isinstance(item, dict) and 'id' in item:
                            # Extract the 'id' to use as the key, and remove it from the values
                            item_id = item.pop('id')
                            aggregated_data[item_id] = item
                elif isinstance(data, dict):
                    aggregated_data.update(data)
                # -----------------------------------------------------
                
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
    Dynamically maps GSM -> SRR if RNA_MA is True.
    """
    if label_category in labels_dict:
        sample_to_label_map = labels_dict[label_category]
    elif label_category.upper() in labels_dict:
        sample_to_label_map = labels_dict[label_category.upper()]
    else:
        print(f"    ! Warning: Category '{label_category}' not found. All samples will be 'unspecified'.")
        sample_to_label_map = {}

    # -------------------------------------------------------------------------
    # --- RNA-SEQ LOGIC: Translate the GSM label keys to SRR keys dynamically
    # -------------------------------------------------------------------------
    if RNA_MA and sample_to_label_map:
        mapped_dict = {}
        for gsm_key, label_val in sample_to_label_map.items():
            gsm_str = str(gsm_key).upper()
            
            if gsm_str.startswith("GSM"):
                # Fetch SRR runs for this GSM (returns instantly if cached)
                srr_list = get_srr_ids(gsm_str)
                for srr in srr_list:
                    mapped_dict[srr] = label_val
            else:
                mapped_dict[gsm_str] = label_val
                
        # Overwrite map with our new SRR-keyed map
        sample_to_label_map = mapped_dict
    # -------------------------------------------------------------------------

    cleaned_labels = []
    for s in df.index:
        # Handles potential variations in index string (e.g., SRR123456_1)
        if RNA_MA:
            s_upper = str(s).upper().split('_')[1]
        else:
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
    
    tsne = TSNE(
        n_components=n_components, 
        random_state=42, 
        perplexity=perplexity,
        init='pca',      # More stable and faster convergence than 'random'
        verbose=1,       # Prints progress iterations to the console so you know it isn't frozen
        n_jobs=1         # Forces a single thread to prevent VS Code debugger deadlocks
    )
    print('done building TSNE')
    embedding = tsne.fit_transform(pca_embedding)
    print('done embeding TSNE')
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
                            bio_targets: list,
                            experiment_name: str = "Normalization_Comparison"):
    
    print(f"\n[Generating Metric Comparison Plots...]")
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    combined_data = []
    for stage_name, df in metrics_dict.items():
        df_copy = df.copy()
        
        # If the dataframe is in long format, pivot it back to wide format for Seaborn
        if 'Metric' in df_copy.columns and 'Value' in df_copy.columns:
            df_copy = df_copy.pivot(index='Label_Axis', columns='Metric', values='Value').reset_index()
            df_copy = df_copy.rename(columns={'Label_Axis': 'Category'})
            
        df_copy['Stage'] = stage_name
        combined_data.append(df_copy)
    
    plot_df = pd.concat(combined_data, ignore_index=True)

    # --- DEFINE BACKGROUND COLORS ---
    unique_categories = plot_df['Category'].unique()
    # Use a pastel palette so it doesn't clash with the main bar colors
    bg_palette = sns.color_palette("Pastel1", n_colors=len(unique_categories))
    bg_color_dict = dict(zip(unique_categories, bg_palette))

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
    
    # Moved the title up slightly to make room for the new legend
    fig.suptitle('Batch Correction Evaluation & Confounding Check\n(Calculated exclusively on valid known labels)', fontsize=20, fontweight='bold', y=0.98)

    # Set zorder=3 for barplots so they appear ON TOP of the background bands
    sns.barplot(data=plot_df, x='Category', y='Variance_Explained', hue='Stage', ax=axes[0, 0], palette='Set2', zorder=3)
    axes[0, 0].set_title('A. Variance Explained (Higher = More Influence)')
    axes[0, 0].set_ylabel('R² Score')

    sns.barplot(data=plot_df, x='Category', y='KNN_Purity', hue='Stage', ax=axes[0, 1], palette='Set2', zorder=3)
    axes[0, 1].set_title('B. KNN Purity (Higher = Better Local Grouping)')
    axes[0, 1].set_ylabel('Purity Score')
    axes[0, 1].set_ylim(0, 1.1)

    bio_only_df = plot_df[plot_df['Category'] != 'study_id']
    sns.barplot(data=bio_only_df, x='Category', y='Batch_ASW_within_Bio', hue='Stage', ax=axes[0, 2], palette='Set2', zorder=3)
    axes[0, 2].set_title('C. Batch ASW within Bio (Lower = +Study Mixing)')
    axes[0, 2].set_ylabel('Silhouette Score of Batch')

    sns.barplot(data=plot_df, x='Category', y='ARI', hue='Stage', ax=axes[1, 0], palette='Set2', zorder=3)
    axes[1, 0].set_title('D. Adjusted Rand Index (Align. w. clutsers)')
    axes[1, 0].set_ylabel('ARI Score')

    sns.barplot(data=plot_df, x='Category', y='Silhouette', hue='Stage', ax=axes[1, 1], palette='Set2', zorder=3)
    axes[1, 1].set_title('E. Silhouette Score (Higher = Tighter Clusters)')
    axes[1, 1].set_ylabel('Silhouette Score')

    ax_conf = axes[1, 2]
    if not confounding_df.empty:
        sns.barplot(data=confounding_df, x='Variable', y='Cramers_V', ax=ax_conf, color=sns.color_palette("deep")[0], zorder=3)
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

    # --- APPLY BACKGROUNDS & REMOVE X-LABELS ---
    active_axes = axes.flatten()
    for i, ax in enumerate(active_axes):
        ax.set_xlabel('')
        
        # Get the categories from the current axes to paint backgrounds
        ticks = [t.get_text() for t in ax.get_xticklabels()]
        
        for idx, tick_text in enumerate(ticks):
            # Normalization to match Plot F's capitalized 'Variable' names back to base categories
            clean_name = tick_text.replace(" (Missing)", "").lower()
            match_cat = next((c for c in unique_categories if c.lower() == clean_name), None)
            
            if match_cat:
                # Add a vertical shaded band behind the bars
                ax.axvspan(idx - 0.5, idx + 0.5, color=bg_color_dict[match_cat], alpha=0.4, zorder=0)

        # Clear out the text on the X-axis completely
        ax.set_xticks([]) 
        ax.grid(True, axis='y', zorder=1) # Ensure Y gridlines stay visible

        # Handle the existing 'Stage' legend
        if ax.get_legend() is not None:
            if i == 0: 
                ax.legend(title='Pipeline Stage', loc='upper right', bbox_to_anchor=(1.0, 1.05))
            else:
                ax.get_legend().remove()

    # --- ADD THE NEW SECONDARY CATEGORY LEGEND ---
    category_patches = [
        mpatches.Patch(color=color, alpha=0.4, label=cat.replace('_', ' ').capitalize()) 
        for cat, color in bg_color_dict.items()
    ]
    
    # Place it centrally at the very top, under the main title
    fig.legend(handles=category_patches, title='Label Type (Background Group)', 
               loc='upper center', ncol=len(unique_categories), 
               bbox_to_anchor=(0.5, 0.93), frameon=True, fontsize=12)

    # Adjusted top rect boundary (0.88 instead of 0.96) to make space for the new legend
    plt.tight_layout(rect=[0, 0.03, 1, 0.88])
    
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

module_dir = './'
sys.path.append(module_dir)

# --- NEW BULKFORMER IMPORTS ---
import torch
import umap
from tqdm import tqdm
from collections import OrderedDict
from torch_sparse import SparseTensor
from src.bulk.utils.BulkFormer import BulkFormer
from src.bulk.model.config import model_params
from src.constants import STORAGE_DIR
# ------------------------------
# ==========================================
# --- BULKFORMER INTEGRATION ---
# ==========================================
BULKFORMER_FILES = {
    'model_weights': f'{STORAGE_DIR}bulkformer/model/checkpoints_ath/BulkFormer_ath_best.pt',
    'graph_ei':      f'{STORAGE_DIR}bulkformer/graph/G_ath.pt',
    'graph_w':       f'{STORAGE_DIR}bulkformer/graph/G_ath_weight.pt',
    'gene_info':     f'{STORAGE_DIR}bulkformer/gene_metadata/arabidopsis_gene_info.csv',
}

def run_bulkformer(df_aligned: pd.DataFrame, batch_size=8):
    """
    Takes an aligned expression dataframe (Samples x Genes), runs it through 
    the BulkFormer model to extract 128-dim representations, and reduces 
    them to 2D via UMAP for visualization.
    """
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    print(f'  [BulkFormer] Device: {device}')

    # 1. Gene vocabulary mapping
    gene_info = pd.read_csv(BULKFORMER_FILES['gene_info'])
    id_col = 'tair_id' if 'tair_id' in gene_info.columns else gene_info.columns[0]
    all_genes = gene_info[id_col].drop_duplicates().tolist()

    # df_aligned is (Samples x Genes). Sync with known vocabulary
    genes_in_expr = set(df_aligned.columns)
    gene_list = [g for g in all_genes if g in genes_in_expr]
    GENE_LENGTH = len(gene_list)

    # Reorder columns and pad missing genes
    expr_df = df_aligned[[c for c in gene_list if c in df_aligned.columns]].copy()
    missing = [g for g in gene_list if g not in expr_df.columns]
    if missing:
        pad = pd.DataFrame(-10.0, index=expr_df.index, columns=missing)
        expr_df = pd.concat([expr_df, pad], axis=1)
    
    input_df = expr_df[gene_list]
    expr_arr = input_df.values.astype(np.float32)

    # 2. Load Graph
    graph_ei = torch.load(BULKFORMER_FILES['graph_ei'], map_location='cpu', weights_only=False)
    graph_w  = torch.load(BULKFORMER_FILES['graph_w'], map_location='cpu', weights_only=False)
    graph_cpu = SparseTensor(
        row=graph_ei[1], col=graph_ei[0],
        value=graph_w,
        sparse_sizes=(GENE_LENGTH, GENE_LENGTH)
    )

    # 3. Initialize Model
    model_params['graph'] = graph_cpu
    model_params['gene_emb'] = None
    model_params['gene_length'] = GENE_LENGTH
    model_params['dim'] = 128

    model = BulkFormer(**model_params).to(device)
    ckpt = torch.load(BULKFORMER_FILES['model_weights'], map_location='cpu', weights_only=False)
    sd = OrderedDict((k[7:] if k.startswith('module.') else k, v) for k, v in ckpt.items())
    
    model_sd = model.state_dict()
    to_load = {k: v for k, v in sd.items() if k in model_sd and model_sd[k].shape == v.shape}
    model.load_state_dict(to_load, strict=False)
    model.eval()
    
    # Move graph to device ONCE before inference
    model.graph = model.graph.to(device)

    # 4. Inference
    results = []
    with torch.no_grad():
        for i in tqdm(range(0, len(expr_arr), batch_size), desc='  [BulkFormer] Extracting embeddings'):
            batch = torch.tensor(expr_arr[i:i+batch_size], dtype=torch.float32).to(device)
            out = model(batch, mask_prob=0.0, output_expr=False)
            results.append(out.cpu())
            
    embeddings = torch.cat(results, dim=0)   # [Samples, Genes, Dim]
    sample_emb = embeddings.mean(dim=1).numpy() # [Samples, Dim]

    # 5. UMAP Reduction
    print('  [BulkFormer] Computing 2D UMAP projection...')
    reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, random_state=42)
    sample_umap = reducer.fit_transform(sample_emb)
    
    return sample_umap