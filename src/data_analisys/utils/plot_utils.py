import matplotlib.pyplot as plt
import torch
import networkx as nx
import os
import sys
import seaborn
from matplotlib.colors import LogNorm
import pandas as pd
import umap
import numpy as np
from sklearn.manifold import TSNE
from scipy.cluster import hierarchy
from typing import Optional
import re

def plot_losses(losses,output_dir,exp_name,iteration):
    plt.close()
    xs = [x for x in range(len(losses))]
    plt.plot(xs, losses)
    plt.savefig(f'{output_dir}{exp_name}_loss_{iteration}.svg')
    plt.close()

def plot_values(values,output_dir,exp_name):
    xs = [x for x in range(len(values))]
    plt.plot(xs, values)
    plt.savefig(f'{output_dir}/{exp_name}.svg')
    plt.close()

def plot_values_bar(values,output_dir,exp_name,title:str='',y_label='',x_label=''):
    xs = [x for x in range(len(values))]
    plt.bar(xs, values)
    plt.title(title)
    plt.ylabel(y_label)
    plt.xlabel(x_label)
    plt.savefig(f'{output_dir}/{exp_name}.svg')
    plt.close()

def plot_heat_map(df:pd.DataFrame,save_loc:str, name: str,cluster:bool=True,typ:str='png',title='',log_norm:bool = True, col: Optional[pd.DataFrame] = None,col_cluster:bool=False):
    # Create directories if they don't exist
    output_dir = os.path.join(save_loc, 'heat_map')
    os.makedirs(output_dir, exist_ok=True)
    o = sys.getrecursionlimit()
    sys.setrecursionlimit(10000)
    test = seaborn.clustermap(df,row_cluster=cluster, col_cluster=col_cluster,method='complete', norm=LogNorm() if log_norm else None, col_colors=col)
    # plt.title(title, fontsize=24)
    plt.savefig(f'{output_dir}/{name}.{typ}')
    plt.close()
    sys.setrecursionlimit(o)
    return test

def plot_predictions(final_predictions,final_targets,N,output_dir,exp_name,x_name='',y_name=''):
    xs_ = [x for x in range(len(final_predictions[0]))]
    for plo in range(len(final_predictions)):
        torch.save((final_predictions[plo]/N), output_dir+f'saves/{plo}_pred.pt')
        torch.save((final_targets[plo]/N), output_dir+f'saves/{plo}_targ.pt')
        plt.xlabel(x_name)
        plt.ylabel(y_name)
        plt.plot(xs_, (final_predictions[plo]/N),label='prediction')
        plt.plot(xs_, (final_targets[plo]/N), label = 'target')
        plt.legend(loc="upper left")
        plt.savefig(output_dir+exp_name+'_'+str(plo)+'.svg')
        plt.close()
        plt.xlabel(x_name)
        plt.ylabel(y_name)
        plt.plot(xs_, (final_predictions[plo]/N),label='prediction')
        plt.legend(loc="upper left")
        plt.savefig(output_dir+'pred_'+exp_name+'_'+str(plo)+'.svg')
        plt.close()

def plot_weights(plot_mat,output_dir,exp_name,x_name='',y_name=''):
    for i,w in enumerate(plot_mat):
        m = [x for x in range(len(w))]
        plt.xlabel(x_name)
        plt.ylabel(y_name)
        plt.plot(m, w.cpu(),label='weights')
        plt.legend(loc="upper left")
        plt.savefig(output_dir+exp_name+'_weights_'+str(i)+'.svg')
        plt.close()
    # plt.savefig(output_dir+exp_name+'_weights_overlap.svg')
    # plt.close()

def plot_matrix_graph(plot_mat_bool,output_dir,exp_name,remove_isolated:bool = False,x_name='',y_name=''):
    G = nx.from_numpy_array(plot_mat_bool.numpy())
    if remove_isolated:
        isolated_nodes = list(nx.isolates(G))
        G.remove_nodes_from(isolated_nodes)
    nx.draw(G,node_size=30, alpha = 0.5)
    plt.axis('equal')
    plt.xlabel(x_name)
    plt.ylabel(y_name)
    plt.savefig(output_dir+exp_name+'_graph.svg')
    plt.close()

def plot_matrix(mat,location,name,x_name='',y_name='',title=''):
    plt.imshow(mat, cmap='hot', interpolation='nearest')
    plt.xlabel(x_name)
    plt.ylabel(y_name)
    plt.title(title)
    plt.colorbar()
    plt.savefig(f'{location}{name}.svg')
    plt.close()

def plot_projection(embedding,
    colors: list,
    markers: list,
    title: Optional[str] = "t-SNE Projection",
    name: str = '',
    legend: bool = True,
    save_path: str = ''):

    os.makedirs(save_path, exist_ok=True)

    fig, ax = plt.subplots(figsize=(10, 10))
    # Define marker and color options
    available_markers = ['o', 's', '^', 'D', 'v', 'p', '*', 'X', 'h', '8', 'P', 'd', '>', '<', '1', '2', '3', '4']
    
    # Handle None cases for labels
    if colors is None:
        colors = ['default'] * len(tsne_results)
    if markers is None:
        markers = ['default'] * len(tsne_results)
        
    unique_colors = sorted(set(colors))
    unique_markers = sorted(set(markers))
    
    available_colors = plt.cm.rainbow(np.linspace(0, 1, len(unique_colors)))

    # Create mappings for colors and markers
    color_dict = {color: available_colors[i % len(available_colors)] 
                  for i, color in enumerate(unique_colors)}
    marker_dict = {marker: available_markers[i % len(available_markers)] 
                   for i, marker in enumerate(unique_markers)}
    
    # Create combined groups based on unique color-marker pairs
    unique_combinations = sorted(set(zip(colors, markers)))
    
    # Plot each unique combination
    for color_val, marker_val in unique_combinations:
        # Handle the default case
        if color_val == 'default' and marker_val == 'default':
            mask = np.ones(len(embedding), dtype=bool)
            plot_color = 'b'
            plot_marker = 'o'
            plot_label = None
        else:
            mask = (np.array(colors) == color_val) & (np.array(markers) == marker_val)
            plot_color = color_dict[color_val]
            plot_marker = marker_dict[marker_val]
            # Your original code only labels by color, we'll keep that logic
            plot_label = f'{color_val}' if legend else None

        ax.scatter(
            embedding[mask, 0],
            embedding[mask, 1],
            c=[plot_color] * sum(mask),
            marker=plot_marker,
            label=plot_label
        )
    
    ax.set_aspect('auto', 'datalim','C')
    if title:
        ax.set_title(title, fontsize=24)

    if legend:
        # 1. Get unique handles and labels
        # This de-duplicates the labels (e.g., 'Tissue A' appearing 5 times)
        handles, labels = ax.get_legend_handles_labels()
        unique = dict(zip(labels, handles)) # de-duplicate
        
        # 2. Automatically determine number of columns
        # Aim for max 20 rows per column (adjust as needed)
        num_items = len(unique)
        max_rows_per_col = 40  
        num_cols = (num_items + max_rows_per_col - 1) // max_rows_per_col
        num_cols = max(1, num_cols) # Ensure at least 1 column

        # 3. Place legend outside the plot
        ax.legend(
            unique.values(), 
            unique.keys(),
            loc='upper left',
            bbox_to_anchor=(1.02, 1.0), # (102% width, 100% height)
            borderaxespad=0.0,
            ncol=num_cols,
            fontsize='medium' # 'small' or 'medium' is good for papers
        )
    # plt.axis('square')
    # Construct save path and save figure
    if save_path:
        # We MUST use bbox_inches='tight' so the legend isn't cut off
        fig.savefig(f'{save_path}/{name}.svg', format='svg', bbox_inches='tight')
    
    plt.close(fig) # Close the figure object

def plot_tsne(
    df: pd.DataFrame,
    colors: list,
    markers: list,
    save_path: str,
    title: Optional[str] = "t-SNE Projection",
    name: str = '',
    legend: bool = True):
    """
    Plot t-SNE visualization with an improved legend placed outside the plot.

    Args:
        df: Input DataFrame (samples x features).
        colors: List of labels to be used for color mapping.
        markers: List of labels to be used for marker mapping.
        title: Plot title.
        name: Filename for saving.
        legend: Whether to show a legend.
        save_path: If provided, saves the plot to this path.
    """        
    # Compute t-SNE
    tsne = TSNE(n_components=2, random_state=42)
    tsne_results = tsne.fit_transform(df.values)
    plot_projection(tsne_results,
                    markers=markers,
                    colors=colors,
                    save_path= save_path,
                    title=title,
                    name=name,
                    legend=legend)


def plot_dendogram(embedding,linkage_method,number_of_clusters,figure_out_path,name=''):
        os.makedirs(f'{figure_out_path}/dendogram/', exist_ok=True)
        # 1. Compute linkage matrix
        Z = hierarchy.linkage(embedding, method=linkage_method)

        # 2. Plot dendrogram
        plt.figure(figsize=(12, 8))
        dendro = hierarchy.dendrogram(
            Z,
            truncate_mode='lastp',  # show only the last p merged clusters
            p=number_of_clusters,   # show only these many clusters
            show_leaf_counts=True,  # show number of samples in each cluster
            leaf_rotation=90.,      # rotate labels for better readability
            leaf_font_size=12.,     # font size for labels
            show_contracted=True,   # show contracted branches
            color_threshold=Z[-number_of_clusters+1, 2]  # color threshold for clusters
        )

        plt.title(f'{name} Dendrogram ({linkage_method} linkage)')
        plt.xlabel('Sample index or (cluster size)')
        plt.ylabel('Distance')
        plt.grid(False)
        plt.tight_layout()
        plt.axhline(y=Z[-number_of_clusters+1, 2], color='r', linestyle='--')
        plt.savefig(f'{figure_out_path}/dendogram/dendogram_{name}.svg')
        plt.close()

def evaluate_cluters(clusters,out_path):
    _,counts = np.unique(clusters,return_counts=True)
    plt.plot(counts)
    plt.savefig(f'{out_path}/test.svg')
    return

# new plots:
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import pandas as pd
import os
import re

# --- Configuration ---
name_map = {
    'robust': 'Robust',
    'standardized': 'Standardized',
    'robust+': 'Robust Norm\ntwo ways',
    'standardized+': 'Standardized\ntwo ways',
    '2_way_norm': 'Two-way\nNormalized',
    'study_corrected': 'Study\nCorrected',
    'imputed': 'No correction'
}

color_map = {
    'treatment': '#1f77b4',       # Blue
    'tissue': '#2ca02c',          # Green
    'study': '#d62728',           # Red
    'treatment_on_tissues': "#fbff00",
    'treatment_on_tissues_median': "#9900ff" # Handle potential casing differences
}

def plot_scores(scores: dict, color_map: dict, name_map: dict, title: str, file_name: str, output_dir: str):
    """
    Backbone plotting function using Seaborn.
    
    Args:
        scores (dict): Nested dict {data_name: {score_type: int/float, ...}, ...}
    """
    if not scores:
        print(f"Skipping plot '{title}' - empty data.")
        return

    # 1. Convert Dictionary to Tidy DataFrame for Seaborn
    data_list = []
    order_list = [] # To preserve insertion order of groups
    
    for group_key, metrics in scores.items():
        # Use name_map for the group label immediately, fallback to key if not found
        display_name = name_map.get(group_key, group_key)
        if display_name not in order_list:
            order_list.append(display_name)
            
        for metric_name, score in metrics.items():
            data_list.append({
                'Group': display_name,
                'Metric': metric_name,
                'Score': score
            })

    df = pd.DataFrame(data_list)

    # 2. Setup Plot
    plt.figure(figsize=(16, 8))
    
    # Create the barplot
    # 'hue' automatically creates the legend and grouped bars
    ax = sns.barplot(
        data=df, 
        x='Group', 
        y='Score', 
        hue='Metric', 
        palette=color_map,
        order=order_list,
        edgecolor='black', # Optional: adds definition to bars
        linewidth=0.5
    )

    # 3. Add Dotted Separators
    # In matplotlib/seaborn categorical plots, x-coords are 0, 1, 2...
    # We place lines at 0.5, 1.5, etc.
    for x in range(len(order_list) - 1):
        ax.axvline(x + 0.5, color='gray', linestyle='--', linewidth=1.5, alpha=0.7)

    # 4. Formatting
    ax.set_title(title, fontsize=16)
    ax.set_ylabel("Score", fontsize=14)
    ax.set_xlabel("") # Hide x-axis label (implied by ticks)
    ax.grid(axis='y', linestyle='--', alpha=0.5)
    
    # Add a zero line if data contains negatives (common in relative plots)
    ax.axhline(0, color='black', linewidth=1)

    # Adjust Legend title
    plt.legend(title='Metric', fontsize=12, title_fontsize=12)

    # Adjust layout to prevent label cutoff
    plt.tight_layout()

    # 5. Save
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(os.path.join(output_dir, file_name), bbox_inches='tight')
    plt.close()

def _parse_flat_dict(scores_dict: dict) -> dict:
    """
    Helper: Parses the flat 'regex-based' keys into the nested dictionary structure
    required by plot_scores.
    """
    parsed_data = {}
    
    # Preserve insertion order
    keys_order = []

    for label, score in scores_dict.items():
        # Logic adapted from original regex
        match = re.search(r'(.+?(?:val|study)) (\w+)$', label)
        if not match:
            match = re.search(r'(.+?) (study)$', label)
        
        if not match:
            continue

        matrix_val = match.group(1).strip()
        target = match.group(2).strip().lower() # treatment, tissue, study

        # Clean group key
        group_key = matrix_val.replace(' val', '').replace(' val study', '').replace(' study', '')

        if group_key not in parsed_data:
            parsed_data[group_key] = {}
            keys_order.append(group_key)
        
        parsed_data[group_key][target] = score
        
    # Return dict sorted by original appearance
    return {k: parsed_data[k] for k in keys_order}

def plot_summary_scores(scores_dict: dict, title: str, file_name: str, output_dir: str):
    """
    Wrapper: Parses data and calls the backbone plotting function.
    """
    # 1. Transform flat keys to nested dict
    nested_scores = _parse_flat_dict(scores_dict)
    
    # 2. Call backbone
    plot_scores(
        scores=nested_scores,
        color_map=color_map,
        name_map=name_map,
        title=title,
        file_name=file_name,
        output_dir=output_dir
    )

def plot_summary_scores_relative(scores_dict: dict, title: str, file_name: str, output_dir: str):
    """
    Wrapper: Parses data, calculates Relative Score (Value - Study), 
    removes 'study' column, and calls backbone.
    """
    # 1. Transform flat keys to nested dict
    nested_scores = _parse_flat_dict(scores_dict)
    
    processed_scores = {}

    # 2. Calculate Relative Scores
    for group, metrics in nested_scores.items():
        if 'study' not in metrics:
            # If no study score exists to subtract, skip or copy as is (depending on preference).
            # Here we skip to be safe, or you could handle imputation.
            continue
            
        study_val = metrics['study']
        
        new_metrics = {}
        for key, val in metrics.items():
            if key == 'study':
                continue # Remove study score
            new_metrics[key] = val - study_val
            
        processed_scores[group] = new_metrics

    # 3. Call backbone
    plot_scores(
        scores=processed_scores,
        color_map=color_map,
        name_map=name_map,
        title=title,
        file_name=file_name,
        output_dir=output_dir
    )