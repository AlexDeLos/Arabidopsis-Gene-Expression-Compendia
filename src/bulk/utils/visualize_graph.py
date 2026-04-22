import os
import torch
import pandas as pd
import networkx as nx
import sys
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
from pyvis.network import Network
sys.path.append(os.path.abspath("./"))
from src.constants import STORAGE_DIR

# --- Configuration ---
GRAPH_DIR = f'{STORAGE_DIR}graph_data'
EXPR_PATH = f'{STORAGE_DIR}final_data/filter.csv'
GENE_INFO = './src/bulk/metadata/arabidopsis_gene_info.csv'

# Choose a target gene to visualize its neighborhood
TARGET_GENE = 'AT1G01010' 
HOPS = 5# 1 = direct neighbors, 2 = neighbors of neighbors

def load_gene_mapping():
    """Recreates the exact gene list used in the build script to map indices to TAIR IDs."""
    print("Loading expression matrix to recreate gene mapping (this might take a moment)...")
    gene_info = pd.read_csv(GENE_INFO)
    gene_list = gene_info['tair_id'].tolist()
    
    # We only read 1 row of the expr data to save memory, just to get columns
    expr_cols = pd.read_csv(EXPR_PATH,index_col=0).index.tolist()
    
    # The columns in the build script were transposed, so original row indices became columns.
    # Adjust this if your EXPR_PATH transpose logic differs!
    final_genes = [c for c in gene_list if c in expr_cols]
    
    # Create dict mapping integer index to TAIR ID
    idx_to_gene = {i: gene for i, gene in enumerate(final_genes)}
    gene_to_idx = {gene: i for i, gene in enumerate(final_genes)}
    
    return idx_to_gene, gene_to_idx

def main():
    # 1. Load Data
    print("Loading graph tensors...")
    edge_index = torch.load(f'{GRAPH_DIR}/G_ath_MA.pt')
    edge_weight = torch.load(f'{GRAPH_DIR}/G_ath_weight_MA.pt')
    
    idx_to_gene, gene_to_idx = load_gene_mapping()
    
    if TARGET_GENE not in gene_to_idx:
        raise ValueError(f"Gene {TARGET_GENE} not found in the graph network.")

    # 2. Build the full NetworkX Directed Graph
    print("Building base NetworkX graph...")
    G_full = nx.DiGraph()
    
    # Convert tensors to numpy for fast iteration
    rows = edge_index[0].numpy()
    cols = edge_index[1].numpy()
    weights = edge_weight.numpy()
    
    # Add edges with attributes
    edges = [(int(r), int(c), {'weight': float(w)}) for r, c, w in zip(rows, cols, weights)]
    G_full.add_edges_from(edges)
    
    # 3. Extract Subgraph (Ego Graph)
    target_idx = gene_to_idx[TARGET_GENE]
    print(f"Extracting {HOPS}-hop neighborhood for {TARGET_GENE} (Index {target_idx})...")
    
    # nx.ego_graph gets the neighborhood; undirected=True allows in/out edges
    subgraph = nx.ego_graph(G_full, target_idx, radius=HOPS, undirected=True)
    
    # Relabel nodes from integer indices to TAIR IDs
    subgraph = nx.relabel_nodes(subgraph, idx_to_gene)
    
    print(f"Subgraph size: {subgraph.number_of_nodes()} nodes, {subgraph.number_of_edges()} edges.")

    # 4. Static Plot (Matplotlib)
    plt.figure(figsize=(10, 8))
    pos = nx.spring_layout(subgraph, seed=42) # Spring layout pulls connected nodes together
    
    # Color edges based on positive/negative correlation
    edge_colors = ['red' if subgraph[u][v]['weight'] < 0 else 'blue' for u, v in subgraph.edges()]
    edge_widths = [abs(subgraph[u][v]['weight']) * 2 for u, v in subgraph.edges()]
    
    nx.draw(subgraph, pos, with_labels=True, node_color='lightgreen', node_size=500, 
            font_size=8, edge_color=edge_colors, width=edge_widths, 
            alpha=0.8, arrowsize=10)
    
    plt.title(f"Co-expression Network Neighborhood ({HOPS}-hop) for {TARGET_GENE}")
    plt.savefig('subgraph_static.png', dpi=300, bbox_inches='tight')
    print("Saved static plot to 'subgraph_static.png'")

    # 5. Interactive Plot (PyVis)
    print("Generating interactive HTML...")
    net = Network(height='750px', width='100%', bgcolor='#222222', font_color='white', directed=True)
    
    # PyVis translates NetworkX graphs seamlessly
    net.from_nx(subgraph)
    
    # Customize edge appearance in PyVis
    for edge in net.edges:
        w = edge.get('weight', 0)
        edge['value'] = abs(w)  # Thickness
        edge['color'] = '#ff4d4d' if w < 0 else '#4d94ff'
        edge['title'] = f"PCC: {w:.3f}" # Hover tooltip
        
    # Add physics controls so you can play with the layout in browser
    net.show_buttons(filter_=['physics']) 
    net.save_graph('subgraph_interactive.html')
    print("Saved interactive plot to 'subgraph_interactive.html'")

if __name__ == "__main__":
    main()