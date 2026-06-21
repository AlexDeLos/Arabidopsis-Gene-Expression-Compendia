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

SCRIPT="src/data_analisys/biological_analisys/scripts/gsea_normalization_impact.py"
CONFIG_A="${1:-All_tissues_full_mixed_min_group_0}"
CONFIG_B="${2:-leaf_full_mixed_min_group_0}"
AXIS_B_LABEL="${3:-leaf}"
OUT_DIR="${4:-./gsea_norm_impact_plots}"

mkdir -p ./logs_slurm 2>/dev/null || true

echo "========================================"
echo "Job ID      : ${SLURM_JOB_ID}"
echo "Node        : $(hostname)"
echo "Started     : $(date)"
echo "SIF         : ${SIF}"
echo "Config A    : ${CONFIG_A}"
echo "Config B    : ${CONFIG_B} (${AXIS_B_LABEL})"
echo "Out dir     : ${OUT_DIR}"
echo "========================================"

if [ ! -f "${SIF}" ]; then
    echo "ERROR: Container image not found at ${SIF}"
    find "${BIG_STORAGE}" -name "*.sif" 2>/dev/null
    exit 1
fi
if [ ! -f "${REPO_DIR}/${SCRIPT}" ]; then
    echo "ERROR: script not found at ${REPO_DIR}/${SCRIPT}"
    exit 1
fi
if [ ! -d "${BIG_STORAGE}" ]; then
    echo "ERROR: storage mount not found at ${BIG_STORAGE}"
    exit 1
fi

echo "+++++++++++++++++ NOW RUNNING GSEA NORMALIZATION/TISSUE IMPACT PLOTS +++++++++++++++++"

apptainer exec \
    --bind "$BIG_STORAGE":"$BIG_STORAGE" \
    --bind "$REPO_DIR":/repo \
    --pwd /repo \
    "$SIF" \
    python /repo/"${SCRIPT}" \
        --config-a "${CONFIG_A}" \
        --config-b "${CONFIG_B}" \
        --axis-b-label "${AXIS_B_LABEL}" \
        --out-dir "${OUT_DIR}"

EXIT_CODE=$?

echo "========================================"
echo "Finished    : $(date)"
echo "Exit code   : ${EXIT_CODE}"
echo "========================================"

exit ${EXIT_CODE}
