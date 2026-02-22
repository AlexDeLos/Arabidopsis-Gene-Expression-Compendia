#!/bin/bash

#SBATCH --partition=general
#SBATCH --qos=medium
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --job-name=batch_RNA_seq_processing
#SBATCH --output=/home/nfs/alexdelossanto/Dataset_fusion_Microarray/logs_slurm/stdout-%x-%A_%a.txt
#SBATCH --array=0-0
#SBATCH --export=NONE   # <-- CRITICAL: Prevents your login environment from bleeding into the job

# Print job information
echo "=================================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Node: $(hostname)"
echo "Date: $(date)"
echo "=================================================="

# 1. Clean environment and load modules
module purge
module use /opt/insy/modulefiles
module load miniconda/3.9

# 2. Properly initialize Conda for non-interactive scripts
eval "$(conda shell.bash hook)"
conda activate download_env

# Ensure clean Java environment
unset JAVA_HOME
unset NXF_JAVA_HOME

# Define your base storage paths
export BIG_STORAGE="/tudelft.net/staff-umbrella/GeneExpressionStorage"
export NXF_WORK="$BIG_STORAGE/work"
mkdir -p "$NXF_WORK"

# --- APPTAINER/SINGULARITY NFS CHMOD FIX ---
# 1. Define permanent storage on the network drive
PERM_CACHE="$BIG_STORAGE/singularity_images"
mkdir -p "$PERM_CACHE"

# 2. Set Nextflow and Apptainer to exclusively use the local node's /tmp drive
export NXF_SINGULARITY_CACHEDIR="/tmp/nxf_cache_$USER"
export APPTAINER_CACHEDIR="/tmp/apptainer_cache_$USER"
export APPTAINER_TMPDIR="/tmp/apptainer_tmp_$USER"
export SINGULARITY_CACHEDIR="$APPTAINER_CACHEDIR"
export SINGULARITY_TMPDIR="$APPTAINER_TMPDIR"
mkdir -p "$NXF_SINGULARITY_CACHEDIR" "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR"

# 3. Pre-load any previously downloaded images from NFS to local /tmp
echo "Pre-loading cached Singularity images to local node..."
cp -u "$PERM_CACHE"/*.img "$PERM_CACHE"/*.sif "$NXF_SINGULARITY_CACHEDIR/" 2>/dev/null || true

# Prevent host Conda paths from bleeding into the container
export SINGULARITYENV_PREPEND_PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
export APPTAINERENV_PREPEND_PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"

# Define Batch Size
BATCH_SIZE=5
echo "Environment activated: $CONDA_DEFAULT_ENV"
echo "Python being used: $(which python3)"
echo "Work Directory is set to: $NXF_WORK"
echo "Processing Batch Index: $SLURM_ARRAY_TASK_ID with Size: $BATCH_SIZE"

# --- RUN THE PIPELINE ---
echo "Running RNA SEQ download and processing..."
python3 -u src/data_importing/new_download.py \
    -o "$BIG_STORAGE/" \
    --array_index $SLURM_ARRAY_TASK_ID \
    --batch_size $BATCH_SIZE \
    --rna

# --- CLEANUP & SAVE ---
# 4. Save any newly pulled images back to the NFS for future jobs
echo "Saving Singularity images back to permanent storage..."
cp -u "$NXF_SINGULARITY_CACHEDIR"/*.img "$NXF_SINGULARITY_CACHEDIR"/*.sif "$PERM_CACHE/" 2>/dev/null || true

echo "Job finished successfully."