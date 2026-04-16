import torch
import pandas as pd
import argparse
import os
import json

def load_labels_study(path):
    labels = {}
    for file in os.listdir(path):
        if file.endswith('.json'):
            file_path = os.path.join(path, file)
            study = file_path.split('/')[-1].split('.')[0]
            with open(file_path, 'r') as file:
                data = json.load(file)
            labels[study] = data
    return labels

def generate_labels_csv(json_dir, output_path):
    """
    Loads JSON labels and flattens them into a single DataFrame 
    saved as labels.csv with Sample IDs as the index.
    """
    print("Generating centralized labels.csv...")
    raw_labels = load_labels_study(json_dir)
    
    flattened_data = []
    
    for study_id, samples in raw_labels.items():
        for sample_id, metadata in samples.items():
            # Create a flat dictionary for the row
            row = {'sample_id': sample_id, 'study_id': study_id}
            
            # Handle lists in metadata (join them for limma compatibility)
            for key, value in metadata.items():
                if isinstance(value, list):
                    row[key] = '_'.join(map(str, sorted(value)))
                else:
                    row[key] = value
                    
            flattened_data.append(row)
            
    # Create DataFrame and set sample_id as the index
    labels_df = pd.DataFrame(flattened_data)
    
    # Optional: Fill missing values with 'unspecified'
    labels_df.fillna('unspecified', inplace=True)
    
    labels_df.set_index('sample_id', inplace=True)
    labels_df.to_csv(output_path+'/labels.csv')
    print(f"Saved labels for {len(labels_df)} samples to {output_path}")
    return labels_df


def train_loop(model, opt, loss_fn, dataloader:torch.Tensor, targated_genes:int,is_TF, lambda_reg:float =0.01,regularization_type:str = 'L1', previous_weights = None, kl_reg:float = 0.0, train = True):
    model.train()
    total_loss = 0
    count = 0
    # for i in targated_genes:
    # Iterate over targated genes
    dataloader_ = dataloader.T
    X:torch.Tensor = torch.cat((dataloader_[:,:targated_genes], dataloader_[:,(targated_genes+1):]), axis = 1)
    y:torch.Tensor = dataloader_[:,targated_genes]
    pred:tuple[torch.Tensor] = model(X)
    
    loss = loss_fn(y.reshape([len(y),1]), pred[0])

    # Apply L1 regularization
    if regularization_type == 'L1' and train:
        # TODO: apply the is_tF mask on this model
        l1_norm = sum(p.abs().sum() for p in model.parameters())
        # kl_loss = torch.nn.KLDivLoss(reduction="batchmean")
        # # output = kl_loss()
        # output = new_gene_loss(pred[1], previous_weights)
        # TODO: fix this returning Nan of the second iteration
         
        loss += lambda_reg * l1_norm
    
    # Apply L2 regularization
    elif regularization_type == 'L2' and train:
        l2_norm = sum(p.pow(2).sum() for p in model.parameters())
        # KL divergence
        kl_loss = torch.nn.KLDivLoss(reduction="batchmean")
        output = kl_loss(pred[1], previous_weights)
         
        loss += lambda_reg * l2_norm + output * kl_reg


    opt.zero_grad()
    loss.backward()
    opt.step()
    total_loss += loss.detach().item()
    count+=1
    return total_loss / (count),pred[1].detach(), y, pred[0].reshape([len(pred[0])])



def get_top_k_mask(tensor, k=10):
    """Returns a boolean mask where True indicates a value is in the top-k absolute values of its row."""
    abs_values = torch.abs(tensor)  # Compute absolute values
    # std = torch.std(tensor,dim=1)*2#! use 0 or 1 I belive 1 is rows

    # return abs_values < std
    # Get the k-th largest absolute value per row (smallest in top-k)
    top_k_thresholds = torch.topk(abs_values, k=k, dim=1)[0][:, -1]
    # Expand thresholds to the same shape as the original tensor for comparison
    expanded_thresholds = top_k_thresholds.unsqueeze(1).expand_as(abs_values)
    # True where the absolute value is >= the threshold (i.e., in top-k)
    mask = abs_values >= expanded_thresholds
    return mask

def mapping_names(x:str):
    split = x.split()
    while len(split[-1])<4:
        split[-1] = "0"+split[-1]
    return ' '.join(split)

def train_encoder(model, epochs,data):
    # Validation using MSE Loss function
    loss_function_ecoder = torch.nn.MSELoss()
    
    # Using an Adam Optimizer with lr = 0.1
    optimizer_ecoder = torch.optim.Adam(model.parameters(),
                                lr = 3e-6)
    losses = []
    for epoch in range(epochs):
        # for gene in data:
        reconstructed = model(data)
        loss = loss_function_ecoder(reconstructed, data)
        optimizer_ecoder.zero_grad()
        loss.backward()
        optimizer_ecoder.step()
        losses.append(loss)
    
    return model,losses

def load_data(path:str):
    data = pd.read_csv(path, index_col=0)
    try:
        del data['cluster']
    except:
        pass
    return data

# Print iterations progress
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 4, length = 100, fill = '█', printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    # Print New Line on Complete
    if iteration == total: 
        print()


def get_most_active_index(tensor,k:int):
    std_devs = torch.std(tensor, dim=1)
    # Step 2: Get the indices that would sort the standard deviations in ascending order
    sorted_indices = torch.argsort(std_devs)
    return sorted_indices[:k]

    # Step 3: Use the sorted indices to reorder the rows of the tensor 
    # sorted_tensor = tensor[sorted_indices]

    # return sorted_tensor[:k]



def get_most_active(tensor,k:int):
    std_devs = torch.std(tensor, dim=1)
    
    sorted_indices = torch.argsort(std_devs,descending=True)

    sorted_tensor = tensor[sorted_indices]

    return sorted_tensor[:k],sorted_indices[:k]

def insert_tensor(original_tensor, index, value):
    """Inserts value into original_tensor at index, returning a new tensor of size n+1."""
    # Convert value to a tensor if it's a scalar (e.g., float/int)
    if not isinstance(value, torch.Tensor):
        value = torch.tensor([value], dtype=original_tensor.dtype, device=original_tensor.device)
    
    # Split the original tensor into left and right parts around the insertion point
    left = original_tensor[:index]
    right = original_tensor[index:]
    
    # Concatenate left + value + right
    new_tensor = torch.cat((left, value, right))
    return new_tensor

def parse(Full:bool, selected_genes: int = 1000,folder_name:str = 'default',version:str='no_version',lambda_reg:float = 0.1,number_runs:int=1):
    parser = argparse.ArgumentParser(
                    prog='GRN extraction',
                    description='Trains models to extract GRN from the ML weights',
                    epilog='Text at the bottom of help')

    parser.add_argument('-N','--number_runs', type=int, default=number_runs)
    parser.add_argument('-e','--epochs', type=int, default=50)
    parser.add_argument('-lr','--learning_rate', type=float, default=1e-3)
    parser.add_argument('-p','--percent', type=float, default=0.9)
    parser.add_argument('-f','--full', type=bool, default=Full)
    parser.add_argument('-p_w','--plot_weights', type=bool, default=not Full)
    parser.add_argument('-p_p','--plot_pred', type=bool, default=not Full)
    parser.add_argument('-c','--clustering', type=bool, default=False)
    parser.add_argument('-l_r','--lambda_reg', type=float, default=lambda_reg)
    parser.add_argument('-r','--regularization', type=str, default='L1')
    parser.add_argument('-n_s_g','--N_selected_genes', type=int, default=selected_genes)
    parser.add_argument('-v','--version', type=str, default=version)
    parser.add_argument('-fol','--folder_name', type=str, default=folder_name)


    return parser.parse_args()