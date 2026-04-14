#!/bin/bash
#SBATCH --job-name=rna_processing
#SBATCH --account=testusers
#SBATCH --partition=all
#SBATCH --time=4:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=2
#SBATCH --output=/home/alexdelossanto/linuxhome/Dataset_fusion_Microarray/logs_slurm/stdout-%x-%A_%a.txt

# ── User configuration ─────────────────────────────────────────────────────────
# Set this to the actual location of your .sif file on the shared storage.
# It MUST be on /tudelft.net/... — paths under /home/ are not visible on compute nodes.
SIF=/tudelft.net/staff-umbrella/GeneExpressionStorage/dataset_fusion.sif

# Storage root (bound into the container so all data paths resolve correctly)
STORAGE=/tudelft.net/staff-umbrella/GeneExpressionStorage
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

# Verify the .sif file is reachable from this compute node before attempting exec.
# If this fails the path is wrong — check with: find /tudelft.net -name "*.sif" 2>/dev/null
if [ ! -f "${SIF}" ]; then
    echo "ERROR: Container image not found at ${SIF}"
    echo "Available .sif files under ${STORAGE}:"
    find "${STORAGE}" -name "*.sif" 2>/dev/null
    exit 1
fi
srun apptainer exec \
    --bind "${STORAGE}:${STORAGE}" \
    --bind "${REPO}:${REPO}" \
    "${SIF}" \
    python src/data_importing/scripts/combine_rna.py

echo "Done with data combination. generating new data matrices"

srun apptainer exec \
    --bind "${STORAGE}:${STORAGE}" \
    --bind "${REPO}:${REPO}" \
    "${SIF}" \
    python src/data_importing/rnaseq_norm_and_analysis.py

EXIT_CODE=$?
echo "========================================"
echo "Finished    : $(date)"
echo "Exit code   : ${EXIT_CODE}"
echo "========================================"
exit ${EXIT_CODE}
