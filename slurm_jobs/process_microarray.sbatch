#!/bin/sh

#SBATCH --partition=general
#SBATCH --qos=medium
#SBATCH --time=36:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --job-name=microarray_analysis
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

# # 2. Activate the PACKED environment directly
# # We use 'source' instead of 'conda activate' to ensure it uses the 
# # exact folder where we patched the R paths.
# source /home/nfs/alexdelossanto/.conda/envs/pre_processing_packed/bin/activate

# 3. Environment Variables
export BIG_STORAGE="/tudelft.net/staff-umbrella/GeneExpressionStorage"

# Optional: Force R to use the environment's R if rpy2 gets confused
export R_HOME=$CONDA_PREFIX/lib/R

echo "Environment activated: $CONDA_PREFIX"
echo "Python being used: $(which python3)"

echo "Gethering data..."
python3 -u src/data_importing/new_download.py -o '/tudelft.net/staff-umbrella/GeneExpressionStorage/' --ma
# 4. Run the Script
echo "Starting analysis..."
python3 -u src/data_importing/data_norm_and_analisys.py

echo "Job finished successfully."