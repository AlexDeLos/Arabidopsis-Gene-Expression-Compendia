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

# --- 1. CLEANUP (Critical) ---
# Delete the locked/corrupt cache from your previous failed runs
rm -rf /tudelft.net/staff-umbrella/GeneExpressionStorage/singularity_cache/cache

# --- 2. STORAGE CONFIGURATION (The Fix) ---

# A. Where the FINAL .sif images will be saved (Network Storage - Persistent)
export NXF_SINGULARITY_CACHEDIR="/tudelft.net/staff-umbrella/GeneExpressionStorage/singularity_cache"
mkdir -p "$NXF_SINGULARITY_CACHEDIR"

# B. Where the raw layers are downloaded and built (Local Node - Disposable)
# We point ALL "temp" and "cache" variables to the local node's scratch space.
# This ensures 'chmod' happens locally where it is allowed.
export SINGULARITY_CACHEDIR="$TMPDIR"
export APPTAINER_CACHEDIR="$TMPDIR"
export SINGULARITY_TMPDIR="$TMPDIR"
export APPTAINER_TMPDIR="$TMPDIR"

# --- RUN ---
python3 -u src/data_importing/new_download.py -o '/tudelft.net/staff-umbrella/GeneExpressionStorage/'

echo "done with main.py"