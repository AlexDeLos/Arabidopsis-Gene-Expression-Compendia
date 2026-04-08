#!/bin/sh

#SBATCH --partition=general
#SBATCH --qos=medium
#SBATCH --time=36:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=6G
#SBATCH --job-name=tulip_label_gen
#SBATCH --output=/home/nfs/alexdelossanto/Dataset_fusion_Microarray/logs_slurm/stdout-%x-%j.txt

# Print job information
echo "=================================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Job Name: $SLURM_JOB_NAME"
echo "Node: $(hostname)"
echo "Date: $(date)"
echo "=================================================="
echo ""

echo "Starting TULIP label generation..."

# 2. Run the Script using Apptainer
# echo "running MA labels"
# srun apptainer exec --bind /tudelft.net/staff-umbrella/GeneExpressionStorage:/tudelft.net/staff-umbrella/GeneExpressionStorage /tudelft.net/staff-umbrella/GeneExpressionStorage/dataset_fusion.sif python -u src/meta_data_processing/tulip_label_generation.py --in-folder "/tudelft.net/staff-umbrella/GeneExpressionStorage/processed_microarray_data/"
echo "running RNA labels"
apptainer exec --bind /tudelft.net/staff-umbrella/GeneExpressionStorage:/tudelft.net/staff-umbrella/GeneExpressionStorage /tudelft.net/staff-umbrella/GeneExpressionStorage/dataset_fusion.sif python -u src/meta_data_processing/tulip_label_generation.py --in-folder "/tudelft.net/staff-umbrella/GeneExpressionStorage/rnaseq_data/metadata/"


echo "Job finished successfully."
echo "Date: $(date)"