#!/bin/bash

#SBATCH --time=36:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --mail-type=END
#SBATCH --job-name=batch_RNA_seq_processing
#SBATCH --output=/home/alexdelossanto/linuxhome/Dataset_fusion_Microarray/logs_slurm/stdout-%x-%A_%a.txt
#SBATCH --array=6
#SBATCH --export=NONE

# Print job information
echo "=================================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Node: $(hostname)"
echo "Date: $(date)"
echo "SIF: ${SIF}"
echo "=================================================="

# Ensure clean Java environment
unset JAVA_HOME
unset NXF_JAVA_HOME

SIF=/tudelft.net/staff-umbrella/GeneExpressionStorage/dataset_fusion.sif

# Define your base storage paths
export STORAGE="/tudelft.net/staff-umbrella/GeneExpressionStorage"

#working:
# export NXF_WORK="$STORAGE/work"
# mkdir -p "$NXF_WORK"
#temp (testing):
# Create a unique local scratch directory for THIS specific array job
export LOCAL_SCRATCH="/tmp/nxf_scratch_${USER}_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$LOCAL_SCRATCH"

# Point Nextflow's work and temp directories to the local node's disk
export NXF_WORK="$LOCAL_SCRATCH/work"
export NXF_TEMP="$LOCAL_SCRATCH/temp"
mkdir -p "$NXF_WORK" "$NXF_TEMP"

echo "Local Nextflow Work Directory is set to: $NXF_WORK"

# ==========================================
# --- THE FINAL SINGULARITY CACHE FIX ---
# ==========================================
# Define your new home-based cache directory
# "/home/nfs/alexdelossanto/Dataset_fusion_Microarray/singularity_cache"
# 1. FINAL IMAGES -> Shared Umbrella Storage
# Nextflow will put the compiled .img files here so all nodes can use them.
export MY_SHARED_IMAGES='/home/nfs/alexdelossanto/singularity_cache'
mkdir -p "$MY_SHARED_IMAGES"
export NXF_SINGULARITY_CACHEDIR="$MY_SHARED_IMAGES"

# Increase timeout so the 10 waiting jobs don't crash while the 1st job downloads
export NXF_SINGULARITY_PULLTIMEOUT="120m"

# 2. TEMPORARY LAYERS -> Local Compute Node Storage (/tmp)
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
srun apptainer exec \
    --bind "${STORAGE}:${STORAGE}" \
    --bind "${REPO}:${REPO}" \
    "${SIF}" \
    python3 -u src/data_importing/new_download.py \
    -o "$STORAGE/" \
    --array_index $SLURM_ARRAY_TASK_ID \
    --batch_size $BATCH_SIZE \
    --rna

echo "Job finished"
echo "Date: $(date)"