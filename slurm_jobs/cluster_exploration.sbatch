#!/bin/sh

#SBATCH --partition=general
#SBATCH --qos=medium
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G
#SBATCH --job-name=cluster_exploration
#SBATCH --output=/home/nfs/alexdelossanto/Dataset_fusion_Microarray/logs_slurm/stdout-%x-%j.txt

# Print job information
echo "=================================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Job Name: $SLURM_JOB_NAME"
echo "Node: $(hostname)"
echo "Date: $(date)"
echo "=================================================="
echo ""


echo "Environment activated: $CONDA_PREFIX"
echo "Python being used: $(which python3)"

echo "Analizing data..."
python3 -u src/data_analisys/cluster_explotation_new.py

echo "Job finished successfully."