#!/bin/sh

#SBATCH --partition=general
#SBATCH --qos=medium
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G
#SBATCH --job-name=2_microarray_analysis
#SBATCH --output=/home/nfs/alexdelossanto/Dataset_fusion_Microarray/logs_slurm/stdout-%x-%j.txt

# Print job information
echo "=================================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Job Name: $SLURM_JOB_NAME"
echo "Node: $(hostname)"
echo "Date: $(date)"
echo "=================================================="
echo ""

# # 1. Load Conda Module
# module use /opt/insy/modulefiles
# module load miniconda/3.9

# 3. Environment Variables
export BIG_STORAGE="/tudelft.net/staff-umbrella/GeneExpressionStorage"

# Optional: Force R to use the environment's R if rpy2 gets confused
export R_HOME=$CONDA_PREFIX/lib/R

echo "Environment activated: $CONDA_PREFIX"
echo "Python being used: $(which python3)"
# 4. Run the Script
echo "Starting analysis..."
python3 -u src/data_importing/data_norm_and_analisys.py

echo "Job finished successfully."