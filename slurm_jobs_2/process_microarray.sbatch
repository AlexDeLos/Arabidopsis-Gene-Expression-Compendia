#!/bin/bash
#SBATCH --job-name=microarray_processing
#SBATCH --account=ewi-insy-prb
#SBATCH --partition=insy
#SBATCH --qos=short
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=140G
#SBATCH --output=./logs_slurm/stdout-%x-%A_%a.txt

# ── User configuration ─────────────────────────────────────────────────────────
SIF=/tudelft.net/staff-umbrella/GeneExpressionStorage/fused.sif
STORAGE=/tudelft.net/staff-umbrella/GeneExpressionStorage
REPO="/home/nfs/alexdelossanto/Dataset_fusion_Microarray"

# ADD THESE TWO LINES TO FIX THE UMAP/NUMBA CACHE CRASH:
mkdir -p "${STORAGE}/numba_cache"
export APPTAINERENV_NUMBA_CACHE_DIR="${STORAGE}/numba_cache"
# ──────────────────────────────────────────────────────────────────────────────
# Make sure the log directory exists before the job starts writing to it
mkdir -p "$(dirname "$SLURM_STDOUT")" 2>/dev/null || true
mkdir -p "${STORAGE}/logs"

echo "========================================"
echo "Job ID      : ${SLURM_JOB_ID}"
echo "Node        : $(hostname)"
echo "Started     : $(date)"
echo "SIF         : ${SIF}"
echo "========================================"

if [ ! -f "${SIF}" ]; then
    echo "ERROR: Container image not found at ${SIF}"
    exit 1
fi

echo "Done with data combination. generating new data matrices"

apptainer exec \
    --bind "${STORAGE}:${STORAGE}" \
    --bind "${REPO}:/repo" \
    --pwd "/repo" \
    "${SIF}" \
    python -u src/data_importing/data_norm_and_analisys.py

EXIT_CODE=$?
echo "========================================"
echo "Finished    : $(date)"
echo "Exit code   : ${EXIT_CODE}"
echo "========================================"
exit ${EXIT_CODE}