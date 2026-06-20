#!/bin/bash
#SBATCH --account=ewi-insy-prb
#SBATCH --partition=general
#SBATCH --time=10:00:00
#SBATCH --qos=medium
#SBATCH --mem=200G
#SBATCH --cpus-per-task=1
#SBATCH --output=./logs_slurm/stdout-%x-%A_%a.txt

REPO_DIR="/home/nfs/alexdelossanto/Dataset_fusion_Microarray"
BIG_STORAGE="/tudelft.net/staff-umbrella/GeneExpressionStorage"
SIF="/home/nfs/alexdelossanto/fused.sif"
# ──────────────────────────────────────────────────────────────────────────────

# Make sure the log directory exists before the job starts writing to it
mkdir -p "$(dirname "$SLURM_STDOUT")" 2>/dev/null || true

echo "========================================"
echo "Job ID      : ${SLURM_JOB_ID}"
echo "Node        : $(hostname)"
echo "Started     : $(date)"
echo "SIF         : ${SIF}"
echo "========================================"

# Verify the .sif file is reachable from this compute node before attempting exec.
# If this fails the path is wrong — check with: find /tudelft.net -name "*.sif" 2>/dev/null
if [ ! -f "${SIF}" ]; then
    echo "ERROR: Container image not found at ${SIF}"
    echo "Available .sif files under ${STORAGE}:"
    find "${STORAGE}" -name "*.sif" 2>/dev/null
    exit 1
fi

# apptainer exec --nv \
#     --bind "$BIG_STORAGE":"$BIG_STORAGE" \
#     --bind "$REPO_DIR":/repo \
#     --pwd /repo \
#     "$SIF" \
#     python src/data_analisys/cluster_explotation_new.py

echo "+++++++++++++++++ NOW RUNNING RNA DATA +++++++++++++++++"
apptainer exec --nv \
    --bind "$BIG_STORAGE":"$BIG_STORAGE" \
    --bind "$REPO_DIR":/repo \
    --pwd /repo \
    "$SIF" \
    python src/data_analisys/cluster_explotation_new.py --rna

EXIT_CODE=$?
echo "========================================"
echo "Finished    : $(date)"
echo "Exit code   : ${EXIT_CODE}"
echo "========================================"
exit ${EXIT_CODE}
