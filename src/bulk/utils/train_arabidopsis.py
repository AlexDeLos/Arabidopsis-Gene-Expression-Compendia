import os
import sys
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader, random_split
from torch_geometric.typing import SparseTensor

sys.path.append(os.path.abspath("./"))
from src.bulk.utils.BulkFormer import BulkFormer

from src.constants import GRAPH_PATH,GRAPH_WEIGHT_PATH,GENE_INFO,EXPR_PATH,WEIGHTS_PATH
# big
# # BulkFormer-127M
# model_params = {
#     'dim': 640,
#     "bins": 0,
#     "gb_repeat": 1,
#     "p_repeat": 8,
#     'bin_head': 12,
#     'full_head': 8,
#     'gene_length': 20010
# }

# BulkFormer-37M  Arabidopsis
# model_params = {
#     'dim':         128,
#     'bins':        0,
#     'gb_repeat':   1,
#     'p_repeat':    1,
#     'bin_head':    12,
#     'full_head':   8,
#     'gene_length': 21040   # from your graph build output
# }
# ── DEBUG SETTINGS ───────────────────────────────────────────────────────────
DEBUG = False  # Set to False for the full cluster run
DEBUG_GENES = 1000   # Number of genes to keep from 22600
DEBUG_SAMPLES = 2000  # Number of samples to keep from 13749
# ─────────────────────────────────────────────────────────────────────────────
# ── Config ────────────────────────────────────────────────────────────────────

LOAD_BEST   = True  # Set to True to load existing best weights

DIM         = 640
GB_REPEAT   = 1
P_REPEAT    = 1
FULL_HEAD   = 8
MASK_RATIO  = 0.15
BATCH_SIZE  = 32
LR          = 1e-5
EPOCHS      = 50
DEVICE      = 'cuda' if torch.cuda.is_available() else 'cpu'
print(f'Device: {DEVICE}')
configs = {
    "Graph": GRAPH_PATH,
    "Weights": GRAPH_WEIGHT_PATH,
    "Gene Info": GENE_INFO,
    "Expression": EXPR_PATH,
    "Pre-trained": WEIGHTS_PATH
}

print("\n--- Model Configuration ---")
for label, path in configs.items():
    print(f"{label:15}: {path}")
print("---------------------------\n")# ── Gene vocab ────────────────────────────────────────────────────────────────
gene_info = pd.read_csv(GENE_INFO)
all_genes = gene_info['tair_id'].drop_duplicates().tolist()

# ── Expression data ───────────────────────────────────────────────────────────
if DEBUG:
    # usecols handles the columns (samples). 
    # range(DEBUG_SAMPLES + 1) takes the first N columns + the gene ID index column.
    expr_df = pd.read_csv(
        EXPR_PATH, 
        index_col=0, 
        nrows=DEBUG_GENES, 
        usecols=range(DEBUG_SAMPLES + 1)
    ).T 
else:
    expr_df = pd.read_csv(EXPR_PATH, index_col=0).T
# Sync vocabulary
gene_list = [g for g in all_genes if g in expr_df.columns]

if DEBUG:
    # Subset genes and columns
    gene_list = gene_list[:DEBUG_GENES]
    expr_df = expr_df[gene_list]
    print(f"DEBUG: Using subset of {len(expr_df)} samples and {len(gene_list)} genes")
else:
    expr_df = expr_df[gene_list]
GENE_LENGTH = len(gene_list)
expr_df     = expr_df[gene_list]
# print(f'Vocabulary synced: {GENE_LENGTH} genes (matches graph node count)')

expr_np = expr_df.values.astype(np.float32)
print(f'Final matrix: {expr_np.shape}')

# ── Graph ─────────────────────────────────────────────────────────────────────
print('Loading and filtering graph...')
ei = torch.load(GRAPH_PATH,  weights_only=False)
ew = torch.load(GRAPH_WEIGHT_PATH, weights_only=False)

if DEBUG:
    # IMPORTANT: Use GENE_LENGTH (the actual count), not DEBUG_GENES
    mask = (ei[0] < GENE_LENGTH) & (ei[1] < GENE_LENGTH)
    ei = ei[:, mask]
    ew = ew[mask]
    print(f"DEBUG: Graph filtered to indices < {GENE_LENGTH}")
    print(f"DEBUG: Graph reduced to {ei.shape[1]} edges")

# Now sparse_sizes will match the max indices in ei
graph = SparseTensor(
    row=ei[1], 
    col=ei[0], 
    value=ew, 
    sparse_sizes=(GENE_LENGTH, GENE_LENGTH)
).to(DEVICE)
print(f'Graph: {ei.shape[1]} edges')

# ── Dataset ───────────────────────────────────────────────────────────────────
class ExprDataset(Dataset):
    def __init__(self, expr, mask_ratio=0.15):
        self.expr       = expr
        self.mask_ratio = mask_ratio

    def __len__(self):
        return len(self.expr)

    def __getitem__(self, idx):
        x    = self.expr[idx].copy()
        true = x.copy()
        obs  = np.where(x != 0)[0]
        k    = max(1, int(len(obs) * self.mask_ratio))
        chosen = np.random.choice(obs, size=k, replace=False)
        x[chosen] = -10.0
        mask = np.zeros(len(x), dtype=np.float32)
        mask[chosen] = 1.0
        return (torch.tensor(x,    dtype=torch.float32),
                torch.tensor(true, dtype=torch.float32),
                torch.tensor(mask, dtype=torch.float32))

n_val    = max(1, int(0.1 * len(expr_np)))
n_train  = len(expr_np) - n_val
train_ds, val_ds = random_split(
    ExprDataset(expr_np, MASK_RATIO), [n_train, n_val],
    generator=torch.Generator().manual_seed(42)
)
train_dl = DataLoader(train_ds, batch_size=BATCH_SIZE, shuffle=True,  num_workers=2, pin_memory=True)
val_dl   = DataLoader(val_ds,   batch_size=BATCH_SIZE, shuffle=False, num_workers=2, pin_memory=True)
print(f'Train: {n_train}  Val: {n_val}')

# ── Model ─────────────────────────────────────────────────────────────────────
model = BulkFormer(
    dim=DIM, graph=graph, gene_emb=None,
    gene_length=GENE_LENGTH,
    bin_head=12, full_head=FULL_HEAD,
    bins=0, gb_repeat=GB_REPEAT, p_repeat=P_REPEAT
).to(DEVICE)
if LOAD_BEST and os.path.exists(WEIGHTS_PATH):
    print(f"Loading weights from {WEIGHTS_PATH}...")
    try:
        # map_location ensures it loads correctly even if trained on a different GPU
        model.load_state_dict(torch.load(WEIGHTS_PATH, map_location=DEVICE))
        print("Successfully loaded pre-trained weights.")
    except Exception as e:
        print(f"Failed to load weights: {e}")
        print("Training from scratch instead.")
else:
    print(f'Model: {sum(p.numel() for p in model.parameters())/1e6:.1f}M params — training from scratch')

# ── Training loop ─────────────────────────────────────────────────────────────
optimizer = torch.optim.AdamW(model.parameters(), lr=LR, weight_decay=1e-4)
scheduler = torch.optim.lr_scheduler.OneCycleLR(
    optimizer, max_lr=LR,
    steps_per_epoch=len(train_dl), epochs=EPOCHS, pct_start=0.05
)

def run_epoch(loader, train=True):
    model.train() if train else model.eval()
    total_loss, n_batches = 0.0, 0
    with torch.set_grad_enabled(train):
        for x, true, mask in loader:
            x, true, mask = x.to(DEVICE), true.to(DEVICE), mask.to(DEVICE)
            pred = model(x, mask_prob=MASK_RATIO, output_expr=True)
            loss = ((pred - true) ** 2 * mask).sum() / (mask.sum() + 1e-8)
            if train:
                optimizer.zero_grad()
                loss.backward()
                nn.utils.clip_grad_norm_(model.parameters(), 1.0)
                optimizer.step()
                scheduler.step()
            total_loss += loss.item()
            # print(f'train={loss.item():.4f}', flush=True)
            n_batches  += 1
    return total_loss / n_batches

best_val = float('inf')
for epoch in range(1, EPOCHS + 1):
    train_loss = run_epoch(train_dl, train=True)
    val_loss   = run_epoch(val_dl,   train=False)
    print(f'Epoch {epoch:3d}  train={train_loss:.4f}  val={val_loss:.4f}  LR= {scheduler.get_last_lr()}', flush=True)

    # torch.save(model.state_dict(), f'{SAVE_DIR}/BulkFormer_ath_epoch{epoch:02d}.pt')
    if val_loss < best_val:
        best_val = val_loss
        torch.save(model.state_dict(), WEIGHTS_PATH)
        print(f'           → new best val: {best_val:.4f}')