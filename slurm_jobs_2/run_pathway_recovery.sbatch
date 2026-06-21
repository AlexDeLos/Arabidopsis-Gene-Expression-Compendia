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

SCRIPT="src/data_analisys/biological_analisys/scripts/pathway_recovery_analysis.py"
CONFIG="${1:-All_tissues_full_mixed_min_group_0}"
FDR="${2:-0.25}"
OUT_DIR="${3:-./pathway_recovery_plots}"
OUT_CSV="${4:-pathway_recovery_table.csv}"

mkdir -p ./logs_slurm 2>/dev/null || true

echo "========================================"
echo "Job ID      : ${SLURM_JOB_ID}"
echo "Node        : $(hostname)"
echo "Started     : $(date)"
echo "SIF         : ${SIF}"
echo "Config      : ${CONFIG}"
echo "FDR cutoff  : ${FDR}"
echo "Out dir     : ${OUT_DIR}"
echo "Out CSV     : ${OUT_CSV}"
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

echo "+++++++++++++++++ NOW RUNNING PATHWAY RECOVERY ANALYSIS +++++++++++++++++"

apptainer exec \
    --bind "$BIG_STORAGE":"$BIG_STORAGE" \
    --bind "$REPO_DIR":/repo \
    --pwd /repo \
    "$SIF" \
    python /repo/"${SCRIPT}" \
        --config "${CONFIG}" \
        --fdr "${FDR}" \
        --out-dir "${OUT_DIR}" \
        --out-csv "${OUT_CSV}"

EXIT_CODE=$?

echo "========================================"
echo "Finished    : $(date)"
echo "Exit code   : ${EXIT_CODE}"
echo "Report dir  : ${OUT_DIR}"
echo "========================================"

exit ${EXIT_CODE}
