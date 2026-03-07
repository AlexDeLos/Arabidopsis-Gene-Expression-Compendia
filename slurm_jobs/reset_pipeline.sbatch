#!/bin/bash

#SBATCH --partition=short
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --job-name=clean_cache_reset_tracker
#SBATCH --output=/home/nfs/alexdelossanto/Dataset_fusion_Microarray/logs_slurm/stdout-%x-%j.txt

echo "=================================================="
echo "Starting Cleanup and Tracker Reset..."
echo "=================================================="

# 1. Clear Nextflow global work cache
echo "-> Clearing Nextflow work directory..."
rm -rf /tudelft.net/staff-umbrella/GeneExpressionStorage/work/*

# 2. Clear corrupted batch processing folders and hidden local caches
echo "-> Clearing local batch folders..."
rm -rf /tudelft.net/staff-umbrella/GeneExpressionStorage/rnaseq_data/batch_processing/batch_GSE*

# 3. Load Conda Environment
echo "-> Loading Python environment..."
module purge
module use /opt/insy/modulefiles
module load miniconda/3.9
eval "$(conda shell.bash hook)"
conda activate download_env

# 4. Safely Reset Tracker Files
echo "-> Resetting file tracker status (4 -> 1)..."

python3 - << 'EOF'
import os
import json
import pandas as pd

tracker_dir = "/tudelft.net/staff-umbrella/GeneExpressionStorage/rnaseq_data/file_tracker/"
updated_files = 0

for file in os.listdir(tracker_dir):
    filepath = os.path.join(tracker_dir, file)
    
    # Handle CSV trackers
    if file.endswith(".csv"):
        try:
            df = pd.read_csv(filepath)
            # Safely replace exact integer 4 or string '4' with 1 across the dataframe
            # This ignores substrings, so 'GSE40216' stays completely safe.
            df = df.replace({4: 1, "4": "1"})
            df.to_csv(filepath, index=False)
            updated_files += 1
        except Exception as e:
            print(f"  Failed to process CSV {file}: {e}")
            
    # Handle JSON trackers
    elif file.endswith(".json"):
        try:
            with open(filepath, 'r') as f:
                data = json.load(f)
                
            # If it's a flat dictionary (e.g. {"GSE123": 4})
            if isinstance(data, dict):
                for k, v in data.items():
                    if v == 4 or v == "4":
                        data[k] = 1
                        
            with open(filepath, 'w') as f:
                json.dump(data, f, indent=4)
            updated_files += 1
        except Exception as e:
            print(f"  Failed to process JSON {file}: {e}")

print(f"-> Successfully scanned and updated {updated_files} tracker file(s).")
EOF

echo "=================================================="
echo "Cleanup complete! Ready for pipeline resubmission."
echo "=================================================="