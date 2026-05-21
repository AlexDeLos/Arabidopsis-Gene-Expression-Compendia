#!/bin/bash
#SBATCH --job-name=new_fastq_download
#SBATCH --account=ewi-insy-prb
#SBATCH --partition=general
#SBATCH --qos=short
#SBATCH --time=4:00:00
#SBATCH --mem=5G
#SBATCH --cpus-per-task=1

# Fix the /repo chdir bug by moving to the safe storage drive before running
cd /tudelft.net/staff-umbrella/GeneExpressionStorage || exit 1
# ══════════════════════════════════════════════════
#  START LOGGING
# ══════════════════════════════════════════════════
BIG_STORAGE="/tudelft.net/staff-umbrella/GeneExpressionStorage"
SIF="$BIG_STORAGE/dataset_fusion.sif"

SRR_LIST=$1
OUT_DIR=$2

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

# ══════════════════════════════════════════════════
#  PURGE INHERITED PARENT ENVIRONMENT VARIABLES
#  Prevents child jobs from inheriting the parent's
#  RAM disk paths or bind mounts (e.g. /repo)
# ══════════════════════════════════════════════════
unset APPTAINER_BIND
unset SINGULARITY_BIND
unset APPTAINER_TMPDIR
unset SINGULARITY_TMPDIR
unset APPTAINER_CACHEDIR
unset SINGULARITY_CACHEDIR

# ══════════════════════════════════════════════════
#  FASTQ-DUMP TEMP DIR ON SHARED STORAGE
#  Avoids /tmp exhaustion on compute nodes for
#  large SRR files
# ══════════════════════════════════════════════════
FASTQ_TMP="$OUT_DIR/tmp_${SRR}"
mkdir -p "$FASTQ_TMP"

while [ $ATTEMPT -lt $MAX_RETRIES ]; do
    echo "Attempt $((ATTEMPT+1)) of $MAX_RETRIES..."

    apptainer exec \
        --no-pid \
        --cleanenv \
        --bind /tudelft.net/staff-umbrella:/tudelft.net/staff-umbrella \
        --env TMPDIR="$FASTQ_TMP" \
        --env HOME="$FASTQ_TMP" \
        "$SIF" \
        fastq-dump --gzip --split-files --outdir "$OUT_DIR" "$SRR"

    FILE_COUNT=$(ls "$OUT_DIR"/"$SRR"*.gz 2>/dev/null | wc -l)

    if [ "$FILE_COUNT" -gt 0 ]; then
        echo "Successfully downloaded $SRR"
        SUCCESS=true
        rm -rf "$FASTQ_TMP"
        break
    else
        echo "Download failed or incomplete for $SRR. Cleaning up partial files..."
        rm -f "$OUT_DIR"/"$SRR"*
        rm -rf "$FASTQ_TMP"
        sleep 10
        ATTEMPT=$((ATTEMPT+1))
    fi
done

# ══════════════════════════════════════════════════
#  END LOGGING
# ══════════════════════════════════════════════════
if [ "$SUCCESS" = false ]; then
    echo "FAILED to download $SRR after $MAX_RETRIES attempts."
    echo "=================================================="
    echo "Job ID:        $SLURM_JOB_ID"
    echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
    echo "Node:          $(hostname)"
    echo "Date:          $(date)"
    echo "=================================================="
    exit 1
fi

echo "=================================================="
echo "Job ID:        $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Node:          $(hostname)"
echo "Date:          $(date)"
echo "=================================================="
