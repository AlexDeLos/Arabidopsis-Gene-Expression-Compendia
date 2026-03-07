#!/bin/bash

#SBATCH --partition=short
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --job-name=clean_cache_reset_tracker
#SBATCH --output=/home/nfs/alexdelossanto/Dataset_fusion_Microarray/logs_slurm/stdout-%x-%j.txt

echo "=================================================="
echo "Starting Cleanup and Tracker Reset..."
echo "=================================================="

# 1. Clear Nextflow global work cache
echo "-> Clearing Nextflow work directory..."
rm -rf /tudelft.net/staff-umbrella/GeneExpressionStorage/work/*

# 2. Clear corrupted batch processing folders and hidden local caches
echo "-> Clearing local batch folders..."
rm -rf /tudelft.net/staff-umbrella/GeneExpressionStorage/rnaseq_data/batch_processing/batch_GSE*


echo "-> Clearing local singularity cache..."
rm -rf /home/nfs/alexdelossanto/singularity_cache

# 3. Reset Tracker Files
echo "-> Resetting file tracker status (4 -> 1)..."

TRACKER_DIR="/tudelft.net/staff-umbrella/GeneExpressionStorage/rnaseq_data/file_tracker/"
cd "$TRACKER_DIR" || { echo "Directory not found!"; exit 1; }

# Since the file ONLY contains a digit, a raw sed substitution is perfectly safe!
grep -l "4" * | xargs --no-run-if-empty sed -i 's/4/1/g'

echo "=================================================="
echo "Cleanup complete! Ready for pipeline resubmission."
echo "=================================================="