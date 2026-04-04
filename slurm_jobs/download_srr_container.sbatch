#!/bin/bash
#SBATCH --job-name=new_fastq_download
#SBATCH --account=testusers
#SBATCH --partition=all
#SBATCH --time=2:30:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=4

BIG_STORAGE="/tudelft.net/staff-umbrella/GeneExpressionStorage"
SIF="$BIG_STORAGE/dataset_fusion.sif"

# Arguments passed by the Python script
SRR_LIST=$1
OUT_DIR=$2

# Extract the specific SRR for this array task ID
SRR=$(awk "NR==${SLURM_ARRAY_TASK_ID}" "$SRR_LIST")

echo "=================================================="
echo "Job ID:        $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Node:          $(hostname)"
echo "Date:          $(date)"
echo "SIF:           $SIF"
echo "SRR:           $SRR"
echo "OUT_DIR:       $OUT_DIR"
echo "=================================================="

if [ -z "$SRR" ]; then
    echo "Error: Could not extract SRR for task $SLURM_ARRAY_TASK_ID"
    exit 1
fi

MAX_RETRIES=3
ATTEMPT=0
SUCCESS=false

while [ $ATTEMPT -lt $MAX_RETRIES ]; do
    echo "Attempt $((ATTEMPT+1)) of $MAX_RETRIES..."

    apptainer exec \
        --no-pid \
        --bind /tudelft.net/staff-umbrella:/tudelft.net/staff-umbrella \
        "$SIF" \
        fastq-dump --gzip --split-files --outdir "$OUT_DIR" "$SRR"

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