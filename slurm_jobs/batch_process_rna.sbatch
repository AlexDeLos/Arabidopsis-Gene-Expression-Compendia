#!/bin/sh

#SBATCH --partition=general
#SBATCH --qos=medium
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --job-name=processing_RNA_seq
#SBATCH --output=/home/nfs/alexdelossanto/Dataset_fusion_Microarray/logs_slurm/stdout-%x-%A_%a.txt
#SBATCH --array=0-19

# Print job information
echo "=================================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Node: $(hostname)"
echo "Date: $(date)"
echo "=================================================="

module use /opt/insy/modulefiles
module load miniconda/3.9

conda activate download_env

# Ensure clean environment
unset JAVA_HOME
unset NXF_JAVA_HOME

# Define your storage paths
export BIG_STORAGE="/tudelft.net/staff-umbrella/GeneExpressionStorage"
export NXF_WORK="$BIG_STORAGE/work"
mkdir -p "$NXF_WORK"

# Singularity Cache Setup (Critical for Nextflow)
export NXF_SINGULARITY_CACHEDIR="$BIG_STORAGE/singularity_images"
export SINGULARITY_TMPDIR="$TMPDIR"
export SINGULARITY_CACHEDIR="$BIG_STORAGE/singularity_images"
mkdir -p "$NXF_SINGULARITY_CACHEDIR"

# Clean previous work for this specific job if needed (Be careful with shared work dirs)
# rm -rf ./work 

# Define Batch Size
BATCH_SIZE=5

echo "Work Directory is set to: $NXF_WORK"
echo "Processing Batch Index: $SLURM_ARRAY_TASK_ID with Size: $BATCH_SIZE"

# --- RUN THE SCRIPT ---
# We pass the array index and the batch size as arguments
python3 -u src/data_importing/new_download.py \
    --array_index $SLURM_ARRAY_TASK_ID \
    --batch_size $BATCH_SIZE