#!/bin/sh

#SBATCH --partition=general
#SBATCH --qos=medium
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=32G
# source ~/.bashrc

module use /opt/insy/modulefiles
module load miniconda/3.9

# Ensure clean environment
unset JAVA_HOME
unset NXF_JAVA_HOME

conda activate download_env

export BIG_STORAGE="/tudelft.net/staff-umbrella/GeneExpressionStorage"

export NXF_WORK="$BIG_STORAGE/work"
mkdir -p "$NXF_WORK"

export NXF_SINGULARITY_CACHEDIR="$TMPDIR"
export SINGULARITY_TMPDIR="$TMPDIR"
export SINGULARITY_CACHEDIR="$TMPDIR"
export APPTAINER_TMPDIR="$TMPDIR"
export APPTAINER_CACHEDIR="$TMPDIR"

rm -rf ./work

echo "Work Directory is set to: $NXF_WORK"
python3 -u src/data_importing/new_download.py -o '/tudelft.net/staff-umbrella/GeneExpressionStorage/'

echo "done with main.py"