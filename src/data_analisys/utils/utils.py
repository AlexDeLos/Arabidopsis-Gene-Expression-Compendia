import pandas as pd
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
