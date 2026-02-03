#!/bin/sh

#SBATCH --partition=general
#SBATCH --qos=medium
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=75G
#SBATCH --array=0-49%5
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
python3 -u src/data_importing/new_download.py -o '/tudelft.net/staff-umbrella/GeneExpressionStorage/' --array_index $SLURM_ARRAY_TASK_ID

echo "done with main.py"