import torch
import torch.nn as nn
import pandas as pd
import numpy as np
from collections import OrderedDict
from tqdm import tqdm
from torch.utils.data import TensorDataset, DataLoader
from torch_sparse import SparseTensor
from src.models.BulkFormer import BulkFormer 

def normalize_to_tpm(count_df, gene_length_dict):
    """Normalized raw counts to log-transformed TPM as per official notebook."""
    gene_names = count_df.columns
    # Default to 1000bp if length is unknown
    gene_lengths_kb = np.array([gene_length_dict.get(gene, 1000) / 1000 for gene in gene_names])
    counts_matrix = count_df.values
    
    rate = counts_matrix / gene_lengths_kb
    sum_per_sample = rate.sum(axis=1)
    sum_per_sample[sum_per_sample == 0] = 1e-6 # Avoid division by zero
    
    tpm = (rate / sum_per_sample.reshape(-1, 1)) * 1e6
    return np.log1p(tpm)

def align_to_vocabulary(log_tpm_matrix, current_genes, target_vocab):
    """Aligns genes and pads missing ones with -10 as per official notebook."""
    df = pd.DataFrame(log_tpm_matrix, columns=current_genes)
    missing_genes = list(set(target_vocab) - set(current_genes))
    
    # Create padding dataframe filled with -10
    padding = pd.DataFrame(np.full((df.shape[0], len(missing_genes)), -10.0), 
                          columns=missing_genes, index=df.index)
    
    aligned_df = pd.concat([df, padding], axis=1)[target_vocab]
    mask_prob = len(missing_genes) / len(target_vocab)
    return aligned_df.values, mask_prob

def get_bulkformer_embeddings(
    count_df, 
    gene_length_dict, 
    target_vocab,
    model_path="src/models/weights/BulkFormer_147M.pt",
    graph_data_path="data/G_tcga.pt",
    graph_weights_path="data/G_tcga_weight.pt",
    gene_emb_path="data/esm2_feature_concat.pt"
):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    
    # 1. Load External Assets (Graph & ESM2 Embeddings)
    # The model uses a SparseTensor for PPI interactions
    g_info = torch.load(graph_data_path, map_location='cpu', weights_only=False)
    g_weights = torch.load(graph_weights_path, map_location='cpu', weights_only=False)
    graph = SparseTensor(row=g_info[1], col=g_info[0], value=g_weights).t().to(device)
    gene_emb = torch.load(gene_emb_path, map_location='cpu', weights_only=False)

    # 2. Instantiate Model with correct Hyperparameters (from PDF/Notebook)
    model = BulkFormer(
        dim=640, 
        graph=graph, 
        gene_emb=gene_emb, 
        gene_length=19393,
        gb_repeat=3, 
        p_repeat=4, # K=4 in Supplementary Table 5
        bin_head=4, 
        full_head=4
    ).to(device)

    # 3. Load Weights with Prefix Cleaning
    ckpt = torch.load(model_path, map_location=device, weights_only=False)
    state_dict = OrderedDict()
    for k, v in ckpt.items():
        name = k[7:] if k.startswith("module.") else k # Strip DDP prefix
        state_dict[name] = v
    model.load_state_dict(state_dict)
    model.eval()

    # 4. Data Preparation
    # Note: count_df should be (Samples x Genes)
    log_tpm = normalize_to_tpm(count_df, gene_length_dict)
    aligned_data, mask_prob = align_to_vocabulary(log_tpm, count_df.columns, target_vocab)
    
    # 5. Batched Inference
    dataset = TensorDataset(torch.tensor(aligned_data, dtype=torch.float32))
    loader = DataLoader(dataset, batch_size=8, shuffle=False)
    
    embeddings = []
    with torch.no_grad(), torch.autocast("cuda", enabled=(device.type == 'cuda')):
        for (X,) in tqdm(loader, desc="Extracting BulkFormer Features"):
            # Model returns per-gene embeddings: [Batch, Genes, Dim]
            gene_output = model(X.to(device), mask_prob=mask_prob, output_expr=False)
            
            # Aggregate to sample-level using 'mean' (matches common exploration use-case)
            sample_emb = torch.mean(gene_output, dim=1)
            embeddings.append(sample_emb.cpu().numpy())
            
    return pd.DataFrame(np.vstack(embeddings), index=count_df.index)