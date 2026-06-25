#!/bin/sh
#SBATCH --job-name=rna_processing_RNA
#SBATCH --account=ewi-insy-prb
#SBATCH --partition=general
#SBATCH --qos=medium
#SBATCH --time=24:00:00
#SBATCH --mem=150G
#SBATCH --job-name=biological_pipeline
#SBATCH --output=./logs_slurm/stdout-%x-%A_%a.txt

# Print job information
echo "=================================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Job Name: $SLURM_JOB_NAME"
echo "Node: $(hostname)"
echo "Date: $(date)"
echo "=================================================="
echo ""

# Define your directories
REPO_DIR="/home/nfs/alexdelossanto/Dataset_fusion_Microarray"
BIG_STORAGE="/tudelft.net/staff-umbrella/GeneExpressionStorage"
SIF=/home/nfs/alexdelossanto/fused.sif

echo "running biological pipeline on RNA"

apptainer exec \
    --bind "$BIG_STORAGE":"$BIG_STORAGE" \
    --bind "$REPO_DIR":/repo \
    --pwd /repo \
    "$SIF" \
    python -u src/data_analisys/biological_analisys/diff_and_GSEA_pipeline.py --rna


echo "Job finished successfully."
echo "Date: $(date)"