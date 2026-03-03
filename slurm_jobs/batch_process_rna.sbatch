#!/bin/bash

#SBATCH --partition=general
#SBATCH --qos=long
#SBATCH --time=168:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --mail-type=END
#SBATCH --job-name=batch_RNA_seq_processing
#SBATCH --output=/home/nfs/alexdelossanto/Dataset_fusion_Microarray/logs_slurm/stdout-%x-%A_%a.txt
#SBATCH --array=0-10
#SBATCH --export=NONE

# 2. Properly initialize Conda for non-interactive scripts
eval "$(conda shell.bash hook)"
conda activate download_env

# Print job information
echo "=================================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Node: $(hostname)"
echo "Date: $(date)"
echo "Environment activated: $CONDA_DEFAULT_ENV"
echo "Python being used: $(which python3)"
echo "=================================================="

# Ensure clean Java environment
# unset JAVA_HOME
# unset NXF_JAVA_HOME

# # Define your base storage paths
# export BIG_STORAGE="/tudelft.net/staff-umbrella/GeneExpressionStorage"
# export BIG_STORAGE_temp='/tmp/alexdelossanto/'
# mkdir -p $BIG_STORAGE_temp
# export NXF_WORK="$BIG_STORAGE/work"
# mkdir -p "$NXF_WORK"

# # ==========================================
# # --- THE FINAL SINGULARITY CACHE FIX ---
# # ==========================================
# # Define your new home-based cache directory
# # "/home/nfs/alexdelossanto/Dataset_fusion_Microarray/singularity_cache"
# export MY_LOCAL_CACHE='/tmp/alexdelossanto/'
# mkdir -p "$MY_LOCAL_CACHE"

# # 1. Tell Nextflow where to store the finalized images (Shared Network Drive, no chmod issues)
# export NXF_SINGULARITY_CACHEDIR="$MY_LOCAL_CACHE"

# # 2. Tell Singularity/Apptainer where to build the layers (Shared Network Drive)
# export SINGULARITY_CACHEDIR="$MY_LOCAL_CACHE"
# export APPTAINER_CACHEDIR="$MY_LOCAL_CACHE"

# # 3. Keep temporary Apptainer extractions strictly in /tmp to prevent network permission crashes
# export APPTAINER_TMPDIR="/tmp/apptainer_tmp_$USER"
# export SINGULARITY_TMPDIR="/tmp/apptainer_tmp_$USER"
# mkdir -p "$APPTAINER_TMPDIR"

# # Prevent host Conda paths from bleeding into the container
# export SINGULARITYENV_PREPEND_PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
# export APPTAINERENV_PREPEND_PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"

# # Define Batch Size
# BATCH_SIZE=10
# echo "Work Directory is set to: $NXF_WORK"
# echo "Processing Batch Index: $SLURM_ARRAY_TASK_ID with Size: $BATCH_SIZE"

# # --- RUN THE PIPELINE ---
# echo "Running RNA SEQ download and processing..."
# python3 -u src/data_importing/new_download.py \
#     -o "$BIG_STORAGE/" \
#     --array_index $SLURM_ARRAY_TASK_ID \
#     --batch_size $BATCH_SIZE \
#     --rna

# echo "Job finished"
# echo "Date: $(date)"