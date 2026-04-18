#!/bin/sh
#SBATCH --account=testusers
#SBATCH --partition=all
#SBATCH --time=50:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --job-name=tulip_label_gen
#SBATCH --output=/home/alexdelossanto/linuxhome/Dataset_fusion_Microarray/logs_slurm/stdout-%x-%A_%a.txt

# Print job information
echo "=================================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Job Name: $SLURM_JOB_NAME"
echo "Node: $(hostname)"
echo "Date: $(date)"
echo "=================================================="
echo ""

echo "Starting TULIP label generation..."

# Define your directories
REPO_DIR="/home/alexdelossanto/linuxhome/Dataset_fusion_Microarray"
BIG_STORAGE="/tudelft.net/staff-umbrella/GeneExpressionStorage"
SIF="$BIG_STORAGE/dataset_fusion.sif"

echo "running Microarray labels"
apptainer exec \
    --bind "$BIG_STORAGE":"$BIG_STORAGE" \
    --bind "$REPO_DIR":/repo \
    --pwd /repo \
    "$SIF" \
    python -u src/meta_data_processing/tulip_label_generation.py \
    --in-folder "$BIG_STORAGE/processed_microarray_data/"

# echo "running RNA labels"
# apptainer exec \
#     --bind "$BIG_STORAGE":"$BIG_STORAGE" \
#     --bind "$REPO_DIR":/repo \
#     --pwd /repo \
#     "$SIF" \
#     python -u src/meta_data_processing/tulip_label_generation.py \
#     --in-folder "$BIG_STORAGE/rnaseq_data/metadata/"



echo "Job finished successfully."
echo "Date: $(date)"