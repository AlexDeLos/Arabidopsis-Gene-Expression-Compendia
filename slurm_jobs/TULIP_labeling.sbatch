#!/bin/sh

#SBATCH --partition=general
#SBATCH --qos=medium
#SBATCH --time=20:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --job-name=tulip_label_gen
#SBATCH --output=/home/nfs/alexdelossanto/Dataset_fusion_Microarray/logs_slurm/stdout-%x-%j.txt

# Print job information
echo "=================================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Job Name: $SLURM_JOB_NAME"
echo "Node: $(hostname)"
echo "Date: $(date)"
echo "=================================================="
echo ""

# 1. Environment Variables
export BIG_STORAGE="/tudelft.net/staff-umbrella/GeneExpressionStorage"

# Load variables from the .env file (automatically exports them)
if [ -f .env ]; then
    echo "Loading .env file..."
    set -a
    source .env
    set +a
else
    echo "Warning: .env file not found!"
fi

echo "Starting TULIP label generation..."

# 2. Run the Script using Apptainer
# Make sure the file is saved as 'tulip_label_generation.py' in the target directory
apptainer run /tudelft.net/staff-umbrella/GeneExpressionStorage/dataset_fusion.sif python3 -u src/meta_data_processing/tulip_label_generation.py

echo "Job finished successfully."
echo "Date: $(date)"