#!/bin/bash
#SBATCH --job-name=fastq_download
#SBATCH --time=04:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=2
#SBATCH --account=ewi-insy-prb
#SBATCH --partition=ewi-insy-prb,prb,ewi-insy,insy

# Properly initialize Conda for non-interactive scripts
eval "$(conda shell.bash hook)"
conda activate download_env
# Arguments passed by the Python script
SRR_LIST=$1
OUT_DIR=$2

# Extract the specific SRR for this array task ID
SRR=$(awk "NR==${SLURM_ARRAY_TASK_ID}" "$SRR_LIST")

if [ -z "$SRR" ]; then
    echo "Error: Could not extract SRR for task $SLURM_ARRAY_TASK_ID"
    exit 1
fi

echo "Starting download for $SRR into $OUT_DIR"

MAX_RETRIES=3
ATTEMPT=0
SUCCESS=false

# Load conda environment if fastq-dump is installed via conda
# source ~/.bashrc
# conda activate download_env

while [ $ATTEMPT -lt $MAX_RETRIES ]; do
    echo "Attempt $((ATTEMPT+1)) of $MAX_RETRIES..."
    
    fastq-dump --gzip --split-files --outdir "$OUT_DIR" "$SRR"
    
    # Check if any .gz files were successfully created
    FILE_COUNT=$(ls "$OUT_DIR"/"$SRR"*.gz 2>/dev/null | wc -l)
    
    if [ "$FILE_COUNT" -gt 0 ]; then
        echo "Successfully downloaded $SRR"
        SUCCESS=true
        break
    else
        echo "Download failed or incomplete for $SRR. Cleaning up partial files..."
        rm -f "$OUT_DIR"/"$SRR"*
        sleep 10
        ATTEMPT=$((ATTEMPT+1))
    fi
done

if [ "$SUCCESS" = false ]; then
    echo "FAILED to download $SRR after $MAX_RETRIES attempts."
    exit 1
fi