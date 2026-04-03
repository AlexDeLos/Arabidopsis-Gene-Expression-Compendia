#!/bin/bash

#SBATCH --time=36:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=6G
#SBATCH --mail-type=END
#SBATCH --job-name=batch_RNA_seq_processing
#SBATCH --output=/home/alexdelossanto/linuxhome/Dataset_fusion_Microarray/logs_slurm/stdout-%x-%A_%a.txt
#SBATCH --array=6
#SBATCH --export=NONE

# Print job information
echo "=================================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Node: $(hostname)"
echo "Date: $(date)"
echo "SIF: ${SIF}"
echo "=================================================="

# Ensure clean Java environment
unset JAVA_HOME
unset NXF_JAVA_HOME

BIG_STORAGE="/tudelft.net/staff-umbrella/GeneExpressionStorage"
SIF="$BIG_STORAGE/dataset_fusion.sif"
NF_BIN="$BIG_STORAGE/bin"

export NXF_WORK="$BIG_STORAGE/work/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
export NXF_SINGULARITY_CACHEDIR="/home/nfs/alexdelossanto/singularity_cache"
export APPTAINER_TMPDIR="/tmp/apptainer_tmp_${USER}"
export SINGULARITY_TMPDIR="$APPTAINER_TMPDIR"
mkdir -p "$NXF_WORK"

# ── Run Python INSIDE the container ──
apptainer exec \
    --bind /tudelft.net/staff-umbrella:/tudelft.net/staff-umbrella \
    --bind /home/nfs/alexdelossanto:/home/nfs/alexdelossanto \
    --bind /usr/bin/sbatch:/usr/bin/sbatch \          # host sbatch
    --bind /usr/bin/squeue:/usr/bin/squeue \          # host squeue (used by NF)
    --bind /usr/bin/scancel:/usr/bin/scancel \
    --bind "$NF_BIN/nextflow":/usr/local/bin/nextflow \ # host nextflow
    --env NXF_WORK="$NXF_WORK" \
    --env NXF_SINGULARITY_CACHEDIR="$NXF_SINGULARITY_CACHEDIR" \
    --env APPTAINER_TMPDIR="$APPTAINER_TMPDIR" \
    --env SINGULARITY_TMPDIR="$SINGULARITY_TMPDIR" \
    --env SLURM_JOB_ID="$SLURM_JOB_ID" \
    --env SLURM_ARRAY_TASK_ID="$SLURM_ARRAY_TASK_ID" \
    --env PATH="/usr/local/bin:/usr/bin:/bin:$NF_BIN" \
    "$SIF" \
    python3 -u src/data_importing/new_download.py \
        -o "$BIG_STORAGE/" \
        --array_index $SLURM_ARRAY_TASK_ID \
        --batch_size 10 --rna

echo "Job finished"
echo "Date: $(date)"