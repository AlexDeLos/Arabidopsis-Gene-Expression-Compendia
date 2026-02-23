import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import LabelEncoder
import umap
from sklearn.manifold import TSNE
import matplotlib
matplotlib.use('Agg')

def prepare_data_structure(df: pd.DataFrame):
    """
    Ensures data is (Samples x Genes). 
    """
    if df.shape[0] > df.shape[1] and df.shape[0] > 10000:
        print(f"  -> Transposing dataframe from {df.shape} to (Samples, Genes)...")
        return df.T
    return df

def align_labels_to_data(df: pd.DataFrame, labels_dict: dict, label_category: str):
    """
    Aligns the dataframe samples with the labels dictionary.
    Includes logic to fallback to Upper Case keys and handle missing samples.
    """
    
    # 1. Access the specific sub-dictionary (Robust Logic)
    if label_category in labels_dict:
        sample_to_label_map = labels_dict[label_category]
    elif label_category.upper() in labels_dict:
        sample_to_label_map = labels_dict[label_category.upper()]
    else:
        print(f"    ! Warning: Category '{label_category}' not found. All samples will be 'unspecified'.")
        sample_to_label_map = {}

    # 2. Extract Labels for ALL samples in the DataFrame
    cleaned_labels = []
    
    for s in df.index:
        # Try to find the sample in the map
        if s in sample_to_label_map:
            raw_val = sample_to_label_map[s]
            
            # Handle List/Set types (common in TREATMENT)
            if isinstance(raw_val, (list, set, tuple)):
                if len(raw_val) == 0:
                    cleaned_labels.append('unspecified')
                elif len(raw_val) == 1:
                    cleaned_labels.append(str(list(raw_val)[0]))
                else:
                    # Join multiple tags: "Heat_Salt"
                    val_str = "_".join(sorted([str(x) for x in raw_val]))
                    cleaned_labels.append(val_str)
            else:
                cleaned_labels.append(str(raw_val))
        else:
            # S is in Data, but not in Labels -> 'unspecified'
            cleaned_labels.append('unspecified')

    # 3. No filtering needed anymore, we use the whole DF
    filtered_df = df
    
    # 4. Encode
    le = LabelEncoder()
    encoded_labels = le.fit_transform(cleaned_labels)
    
    return filtered_df, cleaned_labels, encoded_labels, len(le.classes_)

def run_pca(data, n_components=50):
    """
    Runs PCA to reduce dimensions before t-SNE/UMAP.
    This is crucial for performance and noise reduction.
    """
    print(f"  -> Running PCA (n={n_components})...")
    pca = PCA(n_components=n_components)
    pca_result = pca.fit_transform(data)
    return pca_result, pca.explained_variance_ratio_

def run_umap(data, n_neighbors=15, min_dist=0.1):
    print("  -> Running UMAP...")
    # Change init from default ('spectral') to 'pca' or 'random'
    reducer = umap.UMAP(
        n_neighbors=n_neighbors, 
        min_dist=min_dist, 
        init='pca',
        n_jobs=1
    )
    return reducer.fit_transform(data)

def run_tsne(data, perplexity=30):
    print("  -> Running t-SNE...")
    # n_jobs=-1 uses all processors
    tsne = TSNE(n_components=2, perplexity=perplexity, n_jobs=-1, random_state=42)
    return tsne.fit_transform(data)

def plot_projection(embedding, labels, title, output_path):
    """
    Generic plotting function for 2D embeddings (PCA, UMAP, t-SNE)
    with sample counts in the legend.
    """
    plt.figure(figsize=(10, 8))
    
    # Create DataFrame for Seaborn
    plot_df = pd.DataFrame(embedding, columns=['Dim1', 'Dim2'])
    plot_df['Raw_Label'] = labels
    
    # 1. Calculate the counts for each unique label
    label_counts = plot_df['Raw_Label'].value_counts().to_dict()
    
    # 2. Create a new formatted label string: "LabelName (n=15)"
    plot_df['Label_with_Count'] = plot_df['Raw_Label'].apply(
        lambda x: f"{x} (n={label_counts[x]})"
    )
    
    unique_label_count = len(label_counts)
    
    # Limit number of legend items if too many
    if unique_label_count > 40:
        sns.scatterplot(
            data=plot_df, x='Dim1', y='Dim2', hue='Label_with_Count', 
            alpha=0.6, s=15, legend=False, palette='tab20',
        )
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    else:
        # Use a dynamic palette to avoid errors if unique labels are between 11 and 20
        palette = 'tab20' if unique_label_count > 10 else 'tab10'
        
        sns.scatterplot(
            data=plot_df, x='Dim1', y='Dim2', hue='Label_with_Count', 
            alpha=0.6, s=15, palette=palette
        )
        # 3. Move the larger legend outside the plot box so it doesn't cover data
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
        
    plt.title(title)
    plt.xlabel('Dimension 1')
    plt.ylabel('Dimension 2')
    
    # Use tight_layout so the external legend isn't cut off when saving
    plt.tight_layout() 
    plt.savefig(output_path)
    plt.close()