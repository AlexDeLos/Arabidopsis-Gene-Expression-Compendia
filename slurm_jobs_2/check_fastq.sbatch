#!/bin/bash
#SBATCH --job-name=fastq_check
#SBATCH --account=testusers
#SBATCH --partition=all
#SBATCH --qos=short
#SBATCH --time=0:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=/home/alexdelossanto/linuxhome/Dataset_fusion_Microarray/logs_slurm/fastq_check_%A_%a.txt

GZ_PATH=$1
RESULT_FILE=$2

echo "Checking: ${GZ_PATH}"

# Stage 1: gzip integrity test
if ! gzip -t -q "${GZ_PATH}" 2>/dev/null; then
    echo "CORRUPT"
    echo "CORRUPT" > "${RESULT_FILE}"
    exit 1
fi

# Stage 2: content check
LINE_COUNT=$(zcat "${GZ_PATH}" | head -4 | wc -l)
if [ "${LINE_COUNT}" -lt 4 ]; then
    echo "TRUNCATED"
    echo "CORRUPT" > "${RESULT_FILE}"
    exit 1
fi

echo "OK"
echo "OK" > "${RESULT_FILE}"
exit 0