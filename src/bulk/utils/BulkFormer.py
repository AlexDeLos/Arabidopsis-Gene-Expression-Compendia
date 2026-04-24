import sys

import torch
import torch.nn as nn

module_dir = "./"
sys.path.append(module_dir)
from src.bulk.utils.BulkFormer_block import BulkFormer_block  # noqa: E402
from src.bulk.utils.Rope import PositionalExprEmbedding  # noqa: E402


class BulkFormer(nn.Module):
    def __init__(self, dim, graph, gene_emb, gene_length, bin_head=4, full_head=4, bins=10, gb_repeat=3, p_repeat=1):
        super().__init__()
        self.dim = dim
        self.gene_length = gene_length
        self.graph = graph

        # # === 基因 embedding ===
        self.gene_emb = nn.Parameter(gene_emb)

        # one-hot embedding 初始化
        self.gene_emb_onehot_layer = nn.Embedding(gene_length, dim)
        nn.init.xavier_uniform_(self.gene_emb_onehot_layer.weight)

        self.gene_emb_proj = nn.Sequential(nn.Linear(dim, 4 * dim), nn.ReLU(), nn.Linear(4 * dim, dim))

        self.expr_emb = PositionalExprEmbedding(dim)
        self.x_proj = nn.Sequential(nn.Linear(dim, 4 * dim), nn.ReLU(), nn.Linear(4 * dim, dim))

        # === 主体模块 ===
        self.gb_formers = nn.ModuleList([BulkFormer_block(dim, gene_length, bin_head, full_head, bins, p_repeat) for _ in range(gb_repeat)])

        self.layernorm = nn.LayerNorm(dim)

        # === 样本整体表达 embedding ===
        self.global_expr_proj = nn.Sequential(
            nn.LayerNorm(gene_length),
            nn.Linear(gene_length, 4 * dim), nn.ReLU(), nn.Linear(4 * dim, dim)
        )
        # === 输出头-逐基因预测 ===
        self.head = nn.Sequential(nn.LayerNorm(dim + 3), nn.Linear(dim + 3, dim), nn.ReLU(), nn.Linear(dim, 1), nn.ReLU())

    # def forward(self, x, mask_prob=None, output_expr=False):
    #     b, g = x.shape
    #     x_input = x.clone()

    #     gene_emb_onehot = self.gene_emb_onehot_layer.weight
    #     # === 基因表达嵌入 ===
    #     gene_emb_proj = self.gene_emb_proj(gene_emb_onehot)

    #     x = self.expr_emb(x) + gene_emb_proj + self.global_expr_proj(x).unsqueeze(1).expand(-1, g, -1)

    #     x = self.x_proj(x)

    #     # === 主干 ===
    #     for layer in self.gb_formers:
    #         x = layer(x, self.graph)

    #     # === gene token 输出 ===
    #     gene_emb = self.layernorm(x)

    #     # === 样本统计特征 ===
    #     # 掩码率
    #     mask_scalar = torch.full((b, g, 1), mask_prob or 0.0, device=x.device)
    #     # === 掩码信息 ===
    #     mask_token_val = -10.0
    #     mask = (x_input == mask_token_val).float()  # 1 表示mask基因
    #     valid_mask = 1 - mask  # 1 表示非mask基因

    #     # 平均表达值（排除掩码）  # noqa: RUF003
    #     expr_mean = (x_input * valid_mask).sum(dim=1, keepdim=True) / (valid_mask.sum(dim=1, keepdim=True) + 1e-8)
    #     expr_mean = expr_mean.unsqueeze(-1).expand(-1, g, -1)

    #     # 非零比例
    #     nonzero_ratio = (x_input != 0).float().sum(dim=1, keepdim=True) / g
    #     nonzero_ratio = nonzero_ratio.unsqueeze(-1).expand(-1, g, -1)

    #     # === 拼接特征（不含 CLS）=== [b, g, dim+3]  # noqa: RUF003
    #     gene_emb_output = torch.cat([gene_emb, mask_scalar, expr_mean, nonzero_ratio], dim=-1)

    #     # === 输出预测 ===
    #     pred = self.head(gene_emb_output).squeeze(-1)

    #     # === 均值校正 ===
    #     # 2. 非mask预测均值
    #     pred_valid_mean = (pred * valid_mask).sum(dim=1, keepdim=True) / (valid_mask.sum(dim=1, keepdim=True) + 1e-8)
    #     # 3. 观测均值（真实表达的非mask部分均值）  # noqa: RUF003
    #     observed_mean = (x_input * valid_mask).sum(dim=1, keepdim=True) / (valid_mask.sum(dim=1, keepdim=True) + 1e-8)
    #     # 4. 校正仅作用于mask位置
    #     pred_corrected = pred.clone()
    #     pred_corrected = pred_corrected - mask * (pred_valid_mean - observed_mean)

    #     if output_expr:
    #         return pred_corrected
    #     return gene_emb_output
    def forward(self, x, mask_prob=None, output_expr=False):
        b, g = x.shape
        x_input = x.clone()

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

        x_for_global = x_input.clone()
        x_for_global[x_for_global == -10.0] = 0.0

        gene_emb_onehot = self.gene_emb_onehot_layer.weight
        check("gene_emb_onehot", gene_emb_onehot)

        gene_emb_proj = self.gene_emb_proj(gene_emb_onehot)
        check("gene_emb_proj", gene_emb_proj)

        expr_emb_out = self.expr_emb(x)
        check("expr_emb", expr_emb_out)

        global_proj_out = self.global_expr_proj(x_for_global)
        check("global_expr_proj", global_proj_out)

        x = expr_emb_out + gene_emb_proj + global_proj_out.unsqueeze(1).expand(-1, g, -1)
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
        check("expr_mean", expr_mean)

        nonzero_ratio = (x_input != 0).float().sum(dim=1, keepdim=True) / g
        nonzero_ratio = nonzero_ratio.unsqueeze(-1).expand(-1, g, -1)

        gene_emb_output = torch.cat([gene_emb, mask_scalar, expr_mean, nonzero_ratio], dim=-1)
        check("gene_emb_output", gene_emb_output)

        pred = self.head(gene_emb_output).squeeze(-1)
        check("head_pred", pred)

        pred_valid_mean = (pred * valid_mask).sum(dim=1, keepdim=True) / (valid_mask.sum(dim=1, keepdim=True) + 1e-8)
        observed_mean = (x_input * valid_mask).sum(dim=1, keepdim=True) / (valid_mask.sum(dim=1, keepdim=True) + 1e-8)
        pred_corrected = pred.clone()
        pred_corrected = pred_corrected - mask * (pred_valid_mean - observed_mean)
        check("pred_corrected", pred_corrected)

        if output_expr:
            return pred_corrected
        return gene_emb_output