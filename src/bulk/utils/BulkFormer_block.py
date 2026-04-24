import torch.nn as nn
import torch
from performer_pytorch import Performer
from torch_geometric.nn.conv import GCNConv


class BulkFormer_block(nn.Module):
    def __init__(self, dim, gene_length, bin_head=4, full_head=4, bins=10, p_repeat=1):
        super().__init__()
        self.dim = dim
        self.gene_length = gene_length
        self.bins = bins
        self.p_repeat = p_repeat
        self.bin_head = bin_head
        self.full_head = full_head

        # 图卷积层
        self.g = GCNConv(dim, dim, cached=True, add_self_loops=True)

        # 全局 performer
        self.f = nn.Sequential(*[Performer(dim=self.dim, heads=self.full_head, depth=1, dim_head=self.dim // self.full_head, attn_dropout=0.05, ff_dropout=0.1) for _ in range(self.p_repeat)])

        self.layernorm = nn.LayerNorm(self.dim)

    # def forward(self, x, graph):
    #     x = self.layernorm(x)

    #     gcn_out = self.g(x, graph)
    #     # gcn_out = torch.nan_to_num(gcn_out, nan=0.0, posinf=0.0, neginf=0.0)
    #     x = x + gcn_out
    #     x = self.f(x)
    #     return x
    def forward(self, x, graph):
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

        x = self.layernorm(x)
        check("after_layernorm", x)

        gcn_out = self.g(x, graph)
        check("after_gcnconv", gcn_out)
        x = x + gcn_out
        check("after_gcn_residual", x)

        x = self.f(x)
        check("after_performer", x)
        return x