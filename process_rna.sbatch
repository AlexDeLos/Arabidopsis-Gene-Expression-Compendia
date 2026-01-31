#!/bin/sh

#SBATCH --partition=general
#SBATCH --qos=medium
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=90G
# source ~/.bashrc

module use /opt/insy/modulefiles
module load miniconda/3.9

# Ensure clean environment
unset JAVA_HOME
unset NXF_JAVA_HOME

conda activate download_env

# --- STORAGE CONFIGURATION ---

# 1. Define where the FINAL .sif images will live (Network Storage - Persistent)
export NXF_SINGULARITY_CACHEDIR="/tudelft.net/staff-umbrella/GeneExpressionStorage/singularity_cache"
mkdir -p "$NXF_SINGULARITY_CACHEDIR"

# 2. Define where the TEMPORARY layers/blobs are downloaded (Local Node - Disposable)
# We force this to $TMPDIR so that 'chmod' operations happen on the local disk, avoiding the error.
export SINGULARITY_CACHEDIR="$TMPDIR"
export APPTAINER_CACHEDIR="$TMPDIR"
export SINGULARITY_TMPDIR="$TMPDIR"
export APPTAINER_TMPDIR="$TMPDIR"

# --- RUN ---
python3 -u src/data_importing/new_download.py -o '/tudelft.net/staff-umbrella/GeneExpressionStorage/'

echo "done with main.py"