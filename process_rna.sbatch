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
export BIG_STORAGE="/tudelft.net/staff-umbrella/GeneExpressionStorage/singularity_cache"
mkdir -p "$BIG_STORAGE"

# 1. Store FINAL images on the Big Storage (To save space)
export SINGULARITY_CACHEDIR="$BIG_STORAGE"
export APPTAINER_CACHEDIR="$BIG_STORAGE"
export NXF_SINGULARITY_CACHEDIR="$BIG_STORAGE"

# 2. Do the BUILD/EXTRACT work on the Local Node (To fix permission errors)
# --- CRITICAL CHANGE HERE: Use $TMPDIR, not $BIG_STORAGE ---
export SINGULARITY_TMPDIR="$TMPDIR"
export APPTAINER_TMPDIR="$TMPDIR"

# --- RUN ---
python3 -u src/data_importing/new_download.py -o '/tudelft.net/staff-umbrella/GeneExpressionStorage/'

echo "done with main.py"