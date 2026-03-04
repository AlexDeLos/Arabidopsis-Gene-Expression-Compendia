#!/bin/bash

#SBATCH --partition=general
#SBATCH --qos=long
#SBATCH --time=168:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=6G
#SBATCH --mail-type=END
#SBATCH --job-name=batch_RNA_seq_processing
#SBATCH --output=/home/nfs/alexdelossanto/Dataset_fusion_Microarray/logs_slurm/stdout-%x-%A_%a.txt
#SBATCH --array=0-10
#SBATCH --export=NONE

module purge
module use /opt/insy/modulefiles
module load miniconda/3.9

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
unset JAVA_HOME
unset NXF_JAVA_HOME

# Define your base storage paths
export BIG_STORAGE="/tudelft.net/staff-umbrella/GeneExpressionStorage"
# export BIG_STORAGE_temp='/tmp/alexdelossanto/'
# mkdir -p $BIG_STORAGE_temp
export NXF_WORK="$BIG_STORAGE/work"
mkdir -p "$NXF_WORK"

# ==========================================
# --- THE FINAL SINGULARITY CACHE FIX ---
# ==========================================
# Define your new home-based cache directory
# "/home/nfs/alexdelossanto/Dataset_fusion_Microarray/singularity_cache"
# 1. FINAL IMAGES -> Shared Umbrella Storage
# Nextflow will put the compiled .img files here so all nodes can use them.
export MY_SHARED_IMAGES='/tudelft.net/staff-umbrella/GeneExpressionStorage/singularity_cache'
mkdir -p "$MY_SHARED_IMAGES"
export NXF_SINGULARITY_CACHEDIR="$MY_SHARED_IMAGES"

# 2. TEMPORARY LAYERS -> Local Compute Node Storage
# Singularity will do its downloads here, bypassing the umbrella storage's 
# strict chmod/ACL rules entirely.
export LOCAL_TEMP_CACHE="/tmp/apptainer_cache_$USER"
mkdir -p "$LOCAL_TEMP_CACHE"

export SINGULARITY_CACHEDIR="$LOCAL_TEMP_CACHE"
export APPTAINER_CACHEDIR="$LOCAL_TEMP_CACHE"
export SINGULARITY_TMPDIR="$LOCAL_TEMP_CACHE"
export APPTAINER_TMPDIR="$LOCAL_TEMP_CACHE"

# Prevent host Conda paths from bleeding into the container
export SINGULARITYENV_PREPEND_PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
export APPTAINERENV_PREPEND_PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
# Define Batch Size
BATCH_SIZE=10
echo "Work Directory is set to: $NXF_WORK"
echo "Processing Batch Index: $SLURM_ARRAY_TASK_ID with Size: $BATCH_SIZE"

# --- RUN THE PIPELINE ---
echo "Running RNA SEQ download and processing..."
python3 -u src/data_importing/new_download.py \
    -o "$BIG_STORAGE/" \
    --array_index $SLURM_ARRAY_TASK_ID \
    --batch_size $BATCH_SIZE \
    --rna

echo "Job finished"
echo "Date: $(date)"