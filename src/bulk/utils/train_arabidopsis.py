import os
import sys
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader, random_split
from torch_geometric.utils import add_self_loops

sys.path.append(os.path.abspath("./"))
from src.bulk.utils.BulkFormer import BulkFormer
from src.constants import GRAPH_PATH, GRAPH_WEIGHT_PATH, GENE_INFO, EXPR_PATH, WEIGHTS_PATH

# ── DEBUG SETTINGS ────────────────────────────────────────────────────────────
DEBUG          = False
DEBUG_GENES    = 1000
DEBUG_SAMPLES  = 2000
# ─────────────────────────────────────────────────────────────────────────────

# ── Config ────────────────────────────────────────────────────────────────────
LOAD_BEST  = False
DIM        = 640
GB_REPEAT  = 1
P_REPEAT   = 1
FULL_HEAD  = 8
MASK_RATIO = 0.15
BATCH_SIZE = 8
LR         = 1e-5
EPOCHS     = 50
DEVICE     = 'cuda' if torch.cuda.is_available() else 'cpu'
print(f'Device: {DEVICE}')

print("\n--- Model Configuration ---")
for label, path in [("Graph", GRAPH_PATH), ("Weights", GRAPH_WEIGHT_PATH),
                    ("Gene Info", GENE_INFO), ("Expression", EXPR_PATH),
                    ("Pre-trained", WEIGHTS_PATH)]:
    print(f"{label:15}: {path}")
print("---------------------------\n")

# ── Gene vocab ────────────────────────────────────────────────────────────────
gene_info = pd.read_csv(GENE_INFO)
all_genes = gene_info['tair_id'].drop_duplicates().tolist()

# ── Expression data ───────────────────────────────────────────────────────────
if DEBUG:
    expr_df = pd.read_csv(
        EXPR_PATH, index_col=0,
        nrows=DEBUG_GENES, usecols=range(DEBUG_SAMPLES + 1)
    ).T
else:
    expr_df = pd.read_csv(EXPR_PATH, index_col=0).T

gene_list = [g for g in all_genes if g in expr_df.columns]
if DEBUG:
    gene_list = gene_list[:DEBUG_GENES]
expr_df     = expr_df[gene_list]
GENE_LENGTH = len(gene_list)

print(f"Running with {GENE_LENGTH} genes.")

expr_np = expr_df.values.astype(np.float32)
print(f'Final matrix: {expr_np.shape}')
print(f"NaNs: {np.isnan(expr_np).any()}  Infs: {np.isinf(expr_np).any()}")
print(f"Max: {expr_np.max():.3f}  Min: {expr_np.min():.3f}  Mean: {expr_np.mean():.3f}")
if expr_np.max() > 500:
    print("WARNING: Values are very large — consider log1p normalisation.")

# ── Graph ─────────────────────────────────────────────────────────────────────
print('Loading graph...')
ei = torch.load(GRAPH_PATH,        weights_only=False)
ew = torch.load(GRAPH_WEIGHT_PATH, weights_only=False)

if DEBUG:
    mask = (ei[0] < GENE_LENGTH) & (ei[1] < GENE_LENGTH)
    ei, ew = ei[:, mask], ew[mask]
    

ei, ew = add_self_loops(
    ei, ew, 
    fill_value=1.0,
    num_nodes=GENE_LENGTH
)

graph = (ei.to(DEVICE), ew.to(DEVICE))
print(f"Unique edges: {ei.shape[1]}  "
      f"After dedup: {torch.unique(ei, dim=1).shape[1]}")

# Check edge weight distribution
print(f"Edge weight distribution: "
      f"min={ew.min():.4f} max={ew.max():.4f} "
      f"mean={ew.mean():.4f} std={ew.std():.4f} "
      f"negative: {(ew < 0).sum().item()} / {len(ew)}")

# Check if add_self_loops created duplicate self-loops
self_loop_mask = (ei[0] == ei[1])
print(f"Self-loops: {self_loop_mask.sum().item()}, "
      f"self-loop weights: min={ew[self_loop_mask].min():.4f} "
      f"max={ew[self_loop_mask].max():.4f}")
# Keep the diagnostic prints using ei directly:
deg = torch.zeros(GENE_LENGTH, device=DEVICE)
deg.scatter_add_(0, ei[0].to(DEVICE), torch.ones(ei.shape[1], device=DEVICE))
print(f"Degree stats — max: {deg.max():.0f}, mean: {deg.mean():.1f}, "
      f"p99: {deg.quantile(0.99):.0f}")
isolated = (deg == 0).sum().item()
if isolated > 0:
    print(f"⚠️ CRITICAL: Found {isolated} isolated genes!")
else:
    print("✅ All genes have at least one connection.")
print(f'Graph: {ei.shape[1]} edges')
print(f"Edge weight NaNs: {torch.isnan(ew).any()}  Infs: {torch.isinf(ew).any()}")

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

with torch.no_grad():
    w = model.gene_emb_onehot_layer.weight
    print(f"gene_emb_onehot_layer after init: shape={w.shape} "
          f"nan%={torch.isnan(w).float().mean():.3f} "
          f"min={w.min():.4f} max={w.max():.4f}")

if LOAD_BEST and os.path.exists(WEIGHTS_PATH):
    print(f"Attempting to load weights from {WEIGHTS_PATH}...")
    ckpt     = torch.load(WEIGHTS_PATH, map_location=DEVICE, weights_only=False)
    model_sd = model.state_dict()
    to_load  = {k: v for k, v in ckpt.items()
                if k in model_sd and model_sd[k].shape == v.shape}
    skipped  = [k for k in ckpt if k not in to_load]
    if skipped:
        print(f"  Shape mismatch — skipped {len(skipped)} layers: {skipped[:4]}...")
    model.load_state_dict(to_load, strict=False)
    print(f"  Loaded {len(to_load)}/{len(model_sd)} layers from checkpoint.")
else:
    print(f'Training from scratch — {sum(p.numel() for p in model.parameters())/1e6:.1f}M params')

def _safe_init(module):
    if isinstance(module, nn.Linear):
        nn.init.xavier_uniform_(module.weight)
        if module.bias is not None:
            nn.init.zeros_(module.bias)
    elif isinstance(module, nn.Embedding):
        nn.init.normal_(module.weight, mean=0.0, std=0.02)

layers_loaded = set(to_load.keys()) if (LOAD_BEST and os.path.exists(WEIGHTS_PATH)) else set()

for name, module in model.named_modules():
    params_loaded = any(f'{name}.{p}' in layers_loaded for p in ['weight', 'bias'])
    if not params_loaded:
        _safe_init(module)

print("Weight init complete.")

optimizer = torch.optim.AdamW(model.parameters(), lr=LR, weight_decay=1e-4)
scheduler = torch.optim.lr_scheduler.OneCycleLR(
    optimizer, max_lr=LR,
    steps_per_epoch=len(train_dl), epochs=EPOCHS, pct_start=0.05
)

def run_epoch(loader, train=True, grad_debug=False):
    model.train() if train else model.eval()
    total_loss, n_batches = 0.0, 0
    with torch.set_grad_enabled(train):
        for batch_idx, (x, true, mask) in enumerate(loader):
            x, true, mask = x.to(DEVICE), true.to(DEVICE), mask.to(DEVICE)
            pred = model(x, mask_prob=MASK_RATIO, output_expr=True)
            loss = ((pred - true) ** 2 * mask).sum() / (mask.sum() + 1e-8)

            if train:
                optimizer.zero_grad()
                loss.backward()

                if grad_debug:
                    # Print ALL parameters in forward order — OK and NaN alike
                    # First NaN in this list is the true source
                    print("\n--- Gradient report (forward order) ---")
                    for name, param in model.named_parameters():
                        if param.grad is None:
                            print(f"  NO_GRAD  {name}")
                        elif not torch.isfinite(param.grad).all():
                            nan_pct = (~torch.isfinite(param.grad)).float().mean().item()
                            abs_max = param.grad[torch.isfinite(param.grad)].abs().max().item() if torch.isfinite(param.grad).any() else float('nan')
                            print(f"  NaN      {name}  nan%={nan_pct:.3f}  finite_absmax={abs_max:.4f}")
                        else:
                            abs_max = param.grad.abs().max().item()
                            print(f"  OK       {name}  absmax={abs_max:.6f}")
                    print("--- End gradient report ---\n")
                    # return  # exit after one batch when debugging

                nn.utils.clip_grad_norm_(model.parameters(), 1.0)
                optimizer.step()
                scheduler.step()

            total_loss += loss.item()
            n_batches  += 1

    if n_batches == 0:
        return float('nan')
    return total_loss / n_batches

best_val = float('inf')
for epoch in range(1, EPOCHS + 1):
    # Only do grad debug on first epoch
    grad_debug = False
    train_loss = run_epoch(train_dl, train=True, grad_debug=grad_debug)
    if grad_debug:
        print("Grad debug complete — exiting.")
        break
    val_loss = run_epoch(val_dl, train=False)
    print(f'Epoch {epoch:3d}  train={train_loss:.4f}  val={val_loss:.4f}  LR= {scheduler.get_last_lr()}', flush=True)

    if torch.isfinite(torch.tensor(val_loss)) and val_loss < best_val:
        best_val = val_loss
        torch.save(model.state_dict(), WEIGHTS_PATH)
        print(f'           → new best val: {best_val:.4f}')
