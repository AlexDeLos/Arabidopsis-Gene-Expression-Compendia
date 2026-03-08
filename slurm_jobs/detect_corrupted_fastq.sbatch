#!/bin/bash

#SBATCH --partition=general
#SBATCH --qos=medium
#SBATCH --time=32:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --job-name=check_fastqs
#SBATCH --output=/home/nfs/alexdelossanto/Dataset_fusion_Microarray/logs_slurm/stdout-%x-%j.txt

module purge
module use /opt/insy/modulefiles
module load miniconda/3.9
eval "$(conda shell.bash hook)"
conda activate download_env

echo "=================================================="
echo "Starting STRICT Python Integrity Scan..."
echo "=================================================="

python3 - << 'EOF'
import os
import gzip
from glob import glob
from multiprocessing import Pool

def check_gzip_integrity(file_path):
    try:
        # We must open and read the file to force Python to evaluate every single multi-part stream marker
        with gzip.open(file_path, 'rb') as f:
            # Read in 10MB chunks to avoid filling up the RAM, while verifying the entire file
            while f.read(1024 * 1024 * 10): 
                pass
        return None # File is 100% healthy
    except Exception as e:
        return file_path # File crashed Python, it is silently corrupted!

if __name__ == '__main__':
    directory = '/tudelft.net/staff-umbrella/GeneExpressionStorage/rnaseq_data/fastq_storage'
    print(f"Locating files in {directory}...")
    
    # Recursively find all fastq.gz files
    files = glob(f"{directory}/**/*.fastq.gz", recursive=True)
    print(f"Scanning {len(files)} files using 8 strict Python cores. This may take a while...")
    
    with Pool(8) as p:
        results = p.map(check_gzip_integrity, files)
        
    corrupted = [r for r in results if r is not None]
    
    if corrupted:
        print("\n==================================================")
        for bad_file in corrupted:
            print(f"🚨 CORRUPTED (Python strict-mode): {bad_file} -> DELETING")
            os.remove(bad_file)
        print(f"==================================================")
        print(f"Deleted {len(corrupted)} silently corrupted files.")
    else:
        print("\nAll files passed the strict Python integrity check!")
EOF

echo "=================================================="
echo "Scan Complete!"
echo "=================================================="