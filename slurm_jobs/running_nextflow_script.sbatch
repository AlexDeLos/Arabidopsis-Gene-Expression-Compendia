#!/bin/bash

#SBATCH --partition=general
#SBATCH --qos=medium
#SBATCH --time=12:00:00         # 12 hours is plenty for index building
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8       # 8 cores for building the STAR/Salmon indices
#SBATCH --mem=32G               # High RAM is required to compile STAR indices
#SBATCH --job-name=build_rna_indices
#SBATCH --output=/home/nfs/alexdelossanto/Dataset_fusion_Microarray/logs_slurm/stdout-%x-%j.txt
#SBATCH --export=NONE

# Print job information
echo "=================================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Date: $(date)"
echo "Task: Compiling Nextflow Genome Indices"
echo "=================================================="

# 1. Clean environment and load modules
module purge
module use /opt/insy/modulefiles
module load miniconda/3.9

# 2. Properly initialize Conda
eval "$(conda shell.bash hook)"
conda activate download_env

# Ensure clean Java environment
unset JAVA_HOME
unset NXF_JAVA_HOME

# Define your base storage paths
export BIG_STORAGE="/tudelft.net/staff-umbrella/GeneExpressionStorage"
export NXF_WORK="$BIG_STORAGE/work"
mkdir -p "$NXF_WORK"

# ==========================================
# --- THE FINAL SINGULARITY CACHE FIX ---
# ==========================================
export MY_LOCAL_CACHE="/home/nfs/alexdelossanto/Dataset_fusion_Microarray/singularity_cache"
mkdir -p "$MY_LOCAL_CACHE"

export NXF_SINGULARITY_CACHEDIR="$MY_LOCAL_CACHE"
export SINGULARITY_CACHEDIR="$MY_LOCAL_CACHE"
export APPTAINER_CACHEDIR="$MY_LOCAL_CACHE"

export APPTAINER_TMPDIR="/tmp/apptainer_tmp_$USER"
export SINGULARITY_TMPDIR="/tmp/apptainer_tmp_$USER"
mkdir -p "$APPTAINER_TMPDIR"

export SINGULARITYENV_PREPEND_PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
export APPTAINERENV_PREPEND_PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"

# Run the python script to trigger the build
python src/build_genome_indices.py

echo "Job Finished at $(date)"