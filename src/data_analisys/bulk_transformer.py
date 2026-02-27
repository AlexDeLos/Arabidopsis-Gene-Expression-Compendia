import torch
import pandas as pd
import torch.nn as nn
from src.models.BulkFormer import BulkFormer 

def get_bulkformer_embeddings(count_df, model_path="weights/pretrain.pth"):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    
    # 1. Define Model Parameters
    # These must match the training configuration of the weights you are loading
    gene_length = 19393 
    dim = 640
    
    # 2. Prepare required positional arguments
    # Note: In the repo, 'graph' is typically a protein-protein interaction (PPI) 
    # graph object. If you don't have it, you may need to load it from the 
    # checkpoint or data folder provided in the GitHub.
    # We initialize a placeholder here, but ensure it matches the model's expected type.
    graph = None 
    gene_emb = torch.zeros(gene_length, dim) # Placeholder tensor

    # 3. Corrected Model Instantiation
    model = BulkFormer(
        dim=dim, 
        graph=graph, 
        gene_emb=gene_emb, 
        gene_length=gene_length
    ).to(device)

    # 4. Load State Dict
    # Using weights_only=True is recommended for security in newer PyTorch versions
    checkpoint = torch.load(model_path, map_location=device, weights_only=False)
    
    # If the checkpoint is a dictionary containing the state_dict
    if 'model_state_dict' in checkpoint:
        model.load_state_dict(checkpoint['model_state_dict'])
    else:
        model.load_state_dict(checkpoint)
        
    model.eval()

    # 5. Extract Embeddings
    data_tensor = torch.tensor(count_df.values).float().to(device)
    with torch.no_grad():
        # The model uses the forward pass or a specific 'encode' method.
        # Based on your snippet, you likely want the output of the transformer blocks.
        embeddings = model(data_tensor) 
        
    return pd.DataFrame(embeddings.cpu().numpy(), index=count_df.index)