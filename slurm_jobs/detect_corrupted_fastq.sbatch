#!/bin/bash

#SBATCH --partition=short
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --job-name=check_fastqs
#SBATCH --output=/home/nfs/alexdelossanto/Dataset_fusion_Microarray/logs_slurm/stdout-%x-%j.txt

echo "=================================================="
echo "Starting FastQ Integrity Scan..."
echo "Directory: /tudelft.net/staff-umbrella/GeneExpressionStorage/rnaseq_data/fastq_storage"
echo "=================================================="

# Run the integrity test across all files using 8 cores.
# Healthy files will remain completely silent.
# Errors from corrupted files will be captured directly in the SLURM log.
find /tudelft.net/staff-umbrella/GeneExpressionStorage/rnaseq_data/fastq_storage -type f -name "*.fastq.gz" -print0 | xargs -0 -P 8 gzip -t

echo "=================================================="
echo "Scan Complete!"
echo "=================================================="