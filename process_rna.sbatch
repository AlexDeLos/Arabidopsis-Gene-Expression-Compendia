#!/bin/sh

#SBATCH --partition=general
#SBATCH --qos=medium
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=90G
# source ~/.bashrc

module use /opt/insy/modulefiles
module load miniconda/3.9

# Ensure clean environment
unset JAVA_HOME
unset NXF_JAVA_HOME

conda activate download_env

rm -rf /tudelft.net/staff-umbrella/GeneExpressionStorage/singularity_cache/cache

export SINGULARITY_TMPDIR="$TMPDIR"
export APPTAINER_TMPDIR="$TMPDIR"

export SINGULARITY_CACHEDIR="$TMPDIR"
export APPTAINER_CACHEDIR="$TMPDIR"

export NXF_SINGULARITY_CACHEDIR="$TMPDIR"

python3 -u src/data_importing/new_download.py -o '/tudelft.net/staff-umbrella/GeneExpressionStorage/'

echo "done with main.py"