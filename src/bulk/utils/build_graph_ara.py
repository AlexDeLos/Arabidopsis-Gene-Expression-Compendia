import os
import sys
import numpy as np
import pandas as pd
import torch
from tqdm import tqdm
from torch_geometric.utils import add_self_loops

sys.path.append(os.path.abspath("./"))
from src.constants import GRAPH_PATH, GRAPH_WEIGHT_PATH, GENE_INFO, EXPR_PATH

CHUNK_SIZE = 500
TOP_K      = 20
PCC_THRESH = 0.2

gene_info = pd.read_csv(GENE_INFO)
gene_list = gene_info["tair_id"].tolist()

expr = pd.read_csv(EXPR_PATH, index_col=0).T
expr = expr[[c for c in gene_list if c in expr.columns]]
print(f"Expression matrix from {EXPR_PATH}: {expr.shape}  (samples × genes)")

# Pre-standardise once
X = expr.values.astype(np.float32)
X = X - X.mean(axis=0, keepdims=True)
std = X.std(axis=0, keepdims=True)
std[std == 0] = 1.0
X = X / std

n_samples, G = X.shape
print(f"Computing chunked PCC for {G} genes...")

rows, cols, vals = [], [], []
for i_start in tqdm(range(0, G, CHUNK_SIZE)):
    i_end = min(i_start + CHUNK_SIZE, G)
    chunk = X[:, i_start:i_end]
    pcc_block = (chunk.T @ X) / n_samples  # (chunk_size, G)

    for local_i, global_i in enumerate(range(i_start, i_end)):
        row_pcc = pcc_block[local_i].copy()
        row_abs = np.abs(row_pcc)
        row_abs[global_i] = 0             # exclude self
        top_idx = np.argsort(row_abs)[-TOP_K:]

        for j in top_idx:
            if row_abs[j] >= PCC_THRESH:
                rows.append(global_i)
                cols.append(int(j))
                vals.append(float(row_pcc[j]))

edge_index = torch.tensor([rows, cols], dtype=torch.long)
edge_weight = torch.tensor(vals, dtype=torch.float32)

print(f"Raw graph: {edge_index.shape[1]} edges across {G} genes "
      f"(avg {edge_index.shape[1] / G:.1f} per gene)")

# ── Step 1: add self-loops with weight 1.0 ────────────────────────────────────

edge_index, edge_weight = add_self_loops(
    edge_index, edge_weight,
    fill_value=1.0,
    num_nodes=G
)
print(f"After self-loops: {edge_index.shape[1]} edges")

# ── Step 2: symmetric degree normalization D^{-1/2} A D^{-1/2} ───────────────
# This matches what GCNConv(normalize=True) would compute from binary weights,
# but applied to the PCC values so they can be used with normalize=False.
deg = torch.zeros(G)
deg.scatter_add_(0, edge_index[0], torch.ones(edge_index.shape[1]))

deg_inv_sqrt = deg.pow(-0.5)
deg_inv_sqrt[~torch.isfinite(deg_inv_sqrt)] = 0.0  # guard isolated nodes

edge_weight = deg_inv_sqrt[edge_index[0]] * edge_weight * deg_inv_sqrt[edge_index[1]]

# Safety clamp — PCC can exceed 1.0 slightly due to float32 accumulation
edge_weight = edge_weight.clamp(-1.0, 1.0)

# ── Sanity checks ─────────────────────────────────────────────────────────────
assert not torch.isnan(edge_weight).any(),  "NaNs in edge weights after normalization"
assert not torch.isinf(edge_weight).any(),  "Infs in edge weights after normalization"
assert (deg > 0).all(),                     "Isolated nodes found"

print(f"Edge weight stats after normalization: "
      f"min={edge_weight.min():.4f}  max={edge_weight.max():.4f}  "
      f"mean={edge_weight.mean():.4f}  std={edge_weight.std():.4f}")
print(f"Negative weights: {(edge_weight < 0).sum().item()} / {len(edge_weight)}")

# ── Save ──────────────────────────────────────────────────────────────────────
torch.save(edge_index,  GRAPH_PATH)
torch.save(edge_weight, GRAPH_WEIGHT_PATH)
print(f"Saved normalized graph with self-loops to:\n"
      f"  {GRAPH_PATH}\n  {GRAPH_WEIGHT_PATH}")