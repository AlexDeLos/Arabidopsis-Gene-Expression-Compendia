#!/bin/bash
#SBATCH --account=ewi-insy-prb
#SBATCH --partition=general
#SBATCH --time=00:30:00
#SBATCH --qos=short
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=./logs_slurm/stdout-%x-%A_%a.txt

REPO_DIR="/home/nfs/alexdelossanto/Dataset_fusion_Microarray"
BIG_STORAGE="/tudelft.net/staff-umbrella/GeneExpressionStorage"
SIF="/home/nfs/alexdelossanto/fused.sif"

# Script + args
SCRIPT="src/data_analisys/biological_analisys/scripts/comapre_gsea.py"
CONFIG="${1:-All_tissues_full_mixed_min_group_0}"
FDR="${2:-0.25}"
OUT_CSV="${3:-./gsea_norm_comparison.csv}"

# ──────────────────────────────────────────────────────────────────────────────
# Make sure the log directory exists before the job starts writing to it
mkdir -p ./logs_slurm 2>/dev/null || true

echo "========================================"
echo "Job ID      : ${SLURM_JOB_ID}"
echo "Node        : $(hostname)"
echo "Started     : $(date)"
echo "SIF         : ${SIF}"
echo "Config      : ${CONFIG}"
echo "FDR cutoff  : ${FDR}"
echo "Output CSV  : ${OUT_CSV}"
echo "========================================"

# Verify the .sif file is reachable from this compute node before attempting exec.
# If this fails the path is wrong — check with: find /tudelft.net -name "*.sif" 2>/dev/null
if [ ! -f "${SIF}" ]; then
    echo "ERROR: Container image not found at ${SIF}"
    echo "Available .sif files under ${BIG_STORAGE}:"
    find "${BIG_STORAGE}" -name "*.sif" 2>/dev/null
    exit 1
fi

# Sanity-check the script and storage mount before burning a job slot
if [ ! -f "${REPO_DIR}/${SCRIPT}" ]; then
    echo "ERROR: script not found at ${REPO_DIR}/${SCRIPT}"
    exit 1
fi
if [ ! -d "${BIG_STORAGE}" ]; then
    echo "ERROR: storage mount not found at ${BIG_STORAGE}"
    exit 1
fi

echo "+++++++++++++++++ NOW RUNNING GSEA NORMALIZATION OVERLAP COMPARISON +++++++++++++++++"

apptainer exec \
    --bind "$BIG_STORAGE":"$BIG_STORAGE" \
    --bind "$REPO_DIR":/repo \
    --pwd /repo \
    "$SIF" \
    python /repo/"${SCRIPT}" \
        --config "${CONFIG}" \
        --fdr "${FDR}" \
        --out-csv "${OUT_CSV}"

EXIT_CODE=$?

echo "========================================"
echo "Finished    : $(date)"
echo "Exit code   : ${EXIT_CODE}"
echo "========================================"

exit ${EXIT_CODE}