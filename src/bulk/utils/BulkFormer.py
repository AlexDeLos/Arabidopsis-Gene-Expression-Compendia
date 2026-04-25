import sys

import torch
import torch.nn as nn

module_dir = "./"
sys.path.append(module_dir)
from src.bulk.utils.BulkFormer_block import BulkFormer_block  # noqa: E402
from src.bulk.utils.Rope import PositionalExprEmbedding  # noqa: E402

def check(name, t):
    if not torch.isfinite(t).all():
        finite = t[torch.isfinite(t)]
        nan_pct = (~torch.isfinite(t)).float().mean().item()
        min_val = finite.min().item() if finite.numel() > 0 else float('nan')
        max_val = finite.max().item() if finite.numel() > 0 else float('nan')
        raise RuntimeError(
            f"NaN/Inf in [{name}]  shape={tuple(t.shape)}  "
            f"nan%={nan_pct:.3f}  min={min_val:.3f}  max={max_val:.3f}"
        )


class BulkFormer(nn.Module):
    def __init__(self, 
                 dim, graph, gene_emb, gene_length,
                 bin_head=4, full_head=4, bins=10,
                 gb_repeat=3, p_repeat=1):
        super().__init__()
        self.dim = dim
        self.gene_length = gene_length
        self.graph = graph

        self.gene_emb_onehot_layer = nn.Embedding(gene_length, dim)
        nn.init.xavier_uniform_(self.gene_emb_onehot_layer.weight)

        self.gene_emb_proj = nn.Sequential(
            nn.Linear(dim, 4 * dim),
            nn.ReLU(),
            nn.Linear(4 * dim, dim)
        )

        self.expr_emb = PositionalExprEmbedding(dim)
        self.x_proj = nn.Sequential(
            nn.Linear(dim, 4 * dim),
            nn.ReLU(),
            nn.Linear(4 * dim, dim)
        )

        self.gb_formers = nn.ModuleList([
            BulkFormer_block(dim, gene_length, bin_head, full_head, bins, p_repeat)
            for _ in range(gb_repeat)
        ])

        self.layernorm = nn.LayerNorm(dim)

        # FIX: Linear(1, ...) instead of Linear(gene_length, ...).
        # The old Linear(34858, 4*dim) received gradients ~35x larger than in
        # debug mode (scales with gene_length), blowing up weights on the first
        # backward pass and corrupting every downstream layer including the GCN.
        self.global_expr_proj = nn.Sequential(
            nn.Linear(1, 4 * dim),
            nn.ReLU(),
            nn.Linear(4 * dim, dim)
        )

        self.head = nn.Sequential(
            nn.LayerNorm(dim + 3),
            nn.Linear(dim + 3, dim),
            nn.ReLU(),
            nn.Linear(dim, 1),
            nn.ReLU()
        )

    def forward(self, x, mask_prob=None, output_expr=False):
        b, g = x.shape
        x_input = x.clone()

        gene_emb_onehot = self.gene_emb_onehot_layer.weight
        check("gene_emb_onehot", gene_emb_onehot)

        gene_emb_proj = self.gene_emb_proj(gene_emb_onehot)
        check("gene_emb_proj", gene_emb_proj)

        # FIX: zero out sentinel values (-10.0) before computing the global
        # mean, so masked positions don't corrupt the sample-level summary.
        x_clean = x.clone()
        x_clean[x_clean == -10.0] = 0.0
        global_mean = x_clean.mean(dim=1)                               # (b,)
        global_proj_out = self.global_expr_proj(global_mean.unsqueeze(-1))  # (b, dim)
        check("global_expr_proj", global_proj_out)

        x = self.expr_emb(x) + gene_emb_proj + global_proj_out.unsqueeze(1).expand(-1, g, -1)
        check("after_sum", x)

        x = self.x_proj(x)
        check("x_proj", x)

        for i, layer in enumerate(self.gb_formers):
            x = layer(x, self.graph)
            check(f"gb_former_{i}", x)

        gene_emb = self.layernorm(x)
        check("layernorm", gene_emb)

        mask_scalar = torch.full((b, g, 1), mask_prob or 0.0, device=x.device)

        mask_token_val = -10.0
        mask = (x_input == mask_token_val).float()
        valid_mask = 1 - mask

        expr_mean = (x_input * valid_mask).sum(dim=1, keepdim=True) / (valid_mask.sum(dim=1, keepdim=True) + 1e-8)
        expr_mean = expr_mean.unsqueeze(-1).expand(-1, g, -1)

        nonzero_ratio = (x_input != 0).float().sum(dim=1, keepdim=True) / g
        nonzero_ratio = nonzero_ratio.unsqueeze(-1).expand(-1, g, -1)

        gene_emb_output = torch.cat([gene_emb, mask_scalar, expr_mean, nonzero_ratio], dim=-1)

        pred = self.head(gene_emb_output).squeeze(-1)

        pred_valid_mean = (pred * valid_mask).sum(dim=1, keepdim=True) / (valid_mask.sum(dim=1, keepdim=True) + 1e-8)
        observed_mean = (x_input * valid_mask).sum(dim=1, keepdim=True) / (valid_mask.sum(dim=1, keepdim=True) + 1e-8)
        pred_corrected = pred.clone()
        pred_corrected = pred_corrected - mask * (pred_valid_mean - observed_mean)

        if output_expr:
            return pred_corrected
        else:
            return gene_emb_output