import torch
import torch.nn as nn
from torch_geometric.nn.conv import GCNConv
from performer_pytorch import Performer


class BulkFormer_block(nn.Module):
    def __init__(self, dim, gene_length, bin_head=4, full_head=4, bins=10, p_repeat=1):
        super().__init__()
        self.dim = dim
        self.gene_length = gene_length
        self.bins = bins
        self.p_repeat = p_repeat
        self.bin_head = bin_head
        self.full_head = full_head

        # FIX: cached=False (was True) — caching is unsafe across batches.
        # Edge weights are not passed (None) because GCNConv's normalize=True
        # computes degree from weight sums when weights are provided, which
        # causes NaNs for nodes whose PCC weights nearly cancel. Binary
        # adjacency with normalize=True gives the correct D^{-1/2} A D^{-1/2}.
        self.g = GCNConv(dim, dim, cached=False, add_self_loops=False, normalize=True)

        self.f = nn.Sequential(*[
            Performer(dim=self.dim, heads=self.full_head, depth=1,
                      dim_head=self.dim // self.full_head,
                      attn_dropout=0.05, ff_dropout=0.1)
            for _ in range(self.p_repeat)
        ])

        self.layernorm = nn.LayerNorm(self.dim)

    def forward(self, x, graph_ei, graph_ew):
        # FIX: GCNConv expects (num_nodes, dim) but x is (b, g, dim).
        # Loop over the batch dimension explicitly — passing the 3D tensor
        # directly causes PyG to misinterpret the batch dim as nodes,
        # producing NaNs at full gene-count scale.
        b, g, d = x.shape

        x = self.layernorm(x)

        gcn_out = torch.stack(
            [self.g(x[i], graph_ei, graph_ew) for i in range(b)],
            dim=0
        )  # (b, g, dim)

        x = x + gcn_out

        x = self.f(x)

        return x