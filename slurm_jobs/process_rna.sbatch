#!/bin/sh

#SBATCH --partition=general
#SBATCH --qos=medium
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G
# SBATCH --array=0-0%1
#SBATCH --job-name=processing_RNA_seq
#SBATCH --output=/home/nfs/alexdelossanto/Dataset_fusion_Microarray/logs_slurm/stdout-%x-%j.txt
# source ~/.bashrc

# Print job information
echo "=================================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Job Name: $SLURM_JOB_NAME"
echo "Node: $(hostname)"
echo "Date: $(date)"
echo "=================================================="
echo ""

module use /opt/insy/modulefiles
module load miniconda/3.9

conda activate download_env

# Ensure clean environment
unset JAVA_HOME
unset NXF_JAVA_HOME


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