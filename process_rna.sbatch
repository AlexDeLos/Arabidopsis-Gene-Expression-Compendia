#!/bin/sh

#SBATCH --partition=general
#SBATCH --qos=medium
#SBATCH --mail-type=END
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=90G
# source ~/.bashrc

module use /opt/insy/modulefiles
module load miniconda/3.9

# Ensure the environment is clean before activating
unset JAVA_HOME
unset NXF_JAVA_HOME

conda activate download_env

# --- FIX: CACHE REDIRECTION ---
# 1. Set the path to your large storage folder
export BIG_STORAGE="/tudelft.net/staff-umbrella/GeneExpressionStorage/singularity_cache"
mkdir -p "$BIG_STORAGE"

# 2. Tell Singularity/Apptainer to store cached images there
export SINGULARITY_CACHEDIR="$BIG_STORAGE"
export APPTAINER_CACHEDIR="$BIG_STORAGE"

# 3. Tell Singularity to use it for temporary download files
export SINGULARITY_TMPDIR="$BIG_STORAGE"
export APPTAINER_TMPDIR="$BIG_STORAGE"

# 4. Tell Nextflow to look for images there
export NXF_SINGULARITY_CACHEDIR="$BIG_STORAGE"

# --- RUN ---
# Added '-u' to see logs immediately in the output file
python3 -u src/data_importing/new_download.py -o '/tudelft.net/staff-umbrella/GeneExpressionStorage/'

echo "done with main.py"