#!/bin/sh
#SBATCH --account=testusers
#SBATCH --partition=all
#SBATCH --time=5:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --job-name=biological_pipeline
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

echo "running biological pipeline"
apptainer exec \
    --bind "$BIG_STORAGE":"$BIG_STORAGE" \
    --bind "$REPO_DIR":/repo \
    --pwd /repo \
    "$SIF" \
    python -u src/data_analisys/biological_analisys/diff_and_GSEA_pipeline.py --rna


echo "Job finished successfully."
echo "Date: $(date)"