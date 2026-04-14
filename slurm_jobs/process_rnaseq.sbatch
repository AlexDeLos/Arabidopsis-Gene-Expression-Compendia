#!/bin/bash
#SBATCH --job-name=rna_processing
#SBATCH --account=testusers
#SBATCH --partition=all
#SBATCH --time=4:00:00
#SBATCH --nodes=1                 # <--- ADDED: Explicitly state 1 node
#SBATCH --ntasks=1                # <--- ADDED: Explicitly state 1 task
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --output=/home/alexdelossanto/linuxhome/Dataset_fusion_Microarray/logs_slurm/stdout-%x-%A_%a.txt

# ── User configuration ─────────────────────────────────────────────────────────
SIF=/tudelft.net/staff-umbrella/GeneExpressionStorage/dataset_fusion.sif
STORAGE=/tudelft.net/staff-umbrella/GeneExpressionStorage

# Explicitly define your repository directory!
REPO="/home/alexdelossanto/linuxhome/Dataset_fusion_Microarray"
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

srun --ntasks=1 --cpus-per-task=4 --exact apptainer exec \
    --bind "${STORAGE}:${STORAGE}" \
    --bind "${REPO}:/repo" \
    --pwd "/repo" \
    "${SIF}" \
    python src/data_importing/scripts/combine_rna.py

echo "Done with data combination. generating new data matrices"

srun --ntasks=1 --cpus-per-task=4 --exact apptainer exec \
    --bind "${STORAGE}:${STORAGE}" \
    --bind "${REPO}:/repo" \
    --pwd "/repo" \
    "${SIF}" \
    python src/data_importing/rnaseq_norm_and_analysis.py

EXIT_CODE=$?
echo "========================================"
echo "Finished    : $(date)"
echo "Exit code   : ${EXIT_CODE}"
echo "========================================"
exit ${EXIT_CODE}