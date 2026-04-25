import torch.nn as nn
import torch
from performer_pytorch import Performer
from torch_geometric.nn.conv import GCNConv


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


class BulkFormer_block(nn.Module):
    def __init__(self, dim, gene_length, bin_head=4, full_head=4, bins=10, p_repeat=1):
        super().__init__()
        self.dim = dim
        self.gene_length = gene_length
        self.bins = bins
        self.p_repeat = p_repeat
        self.bin_head = bin_head
        self.full_head = full_head

        # FIX: cached=False — caching is unsafe when called with batched input
        # across different samples. normalize=True applies D^{-1/2} A D^{-1/2}
        # internally so we don't need manual edge weight renormalization.
        self.g = GCNConv(dim, dim, cached=False, add_self_loops=False, normalize=True)

        self.f = nn.Sequential(*[
            Performer(dim=self.dim, heads=self.full_head, depth=1,
                      dim_head=self.dim // self.full_head,
                      attn_dropout=0.05, ff_dropout=0.1)
            for _ in range(self.p_repeat)
        ])

        self.layernorm = nn.LayerNorm(self.dim)

    def forward(self, x, graph_ei, graph_ew):
        # x shape: (b, g, dim)
        b, g, d = x.shape

        x = self.layernorm(x)
        check("after_layernorm", x)
        gcn_out_list = []
        # TEMPORARY DIAGNOSTIC — remove after fixing
        print(f"[GCN input] min={x.min():.3f} max={x.max():.3f} "
            f"nan%={(~torch.isfinite(x)).float().mean():.4f}")

        # Run just node 0 from batch 0 to isolate
        with torch.no_grad():
            test_out = self.g(x[0], graph_ei, graph_ew)
            nan_nodes = (~torch.isfinite(test_out).all(dim=-1)).nonzero(as_tuple=True)[0]
            print(f"[GCN] NaN nodes in sample 0: {len(nan_nodes)}")
            if len(nan_nodes) > 0:
                print(f"  First 5: {nan_nodes[:5].tolist()}")
                # Check the actual weight matrix
                print(f"  GCN weight: min={self.g.lin.weight.min():.3f} "
                    f"max={self.g.lin.weight.max():.3f} "
                    f"nan%={(~torch.isfinite(self.g.lin.weight)).float().mean():.4f}")
                # Check what edge weights look like for a NaN node
                node = nan_nodes[0].item()
                edge_mask = (graph_ei[0] == node)
                print(f"  Node {node} edges: {edge_mask.sum().item()}, "
                    f"weights: {graph_ew[edge_mask][:5].tolist()}")
        raise RuntimeError("Diagnostic complete — stopping")
        for i in range(b):
            node_feat = x[i]                   # (g, dim)
            out_i = self.g(node_feat, graph_ei, graph_ew)   # (g, dim)
            gcn_out_list.append(out_i)
        gcn_out = torch.stack(gcn_out_list, dim=0)   # (b, g, dim)

        check("after_gcnconv", gcn_out)

        # Clamp before residual add to prevent any stray large values
        # from the GCN propagating forward
        gcn_out = torch.clamp(gcn_out, min=-50.0, max=50.0)

        x = x + gcn_out
        check("graph_and_sum", x)

        x = self.f(x)
        check("after_performer", x)

        return x