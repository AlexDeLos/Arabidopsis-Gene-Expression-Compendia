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

# Inputs / outputs for the status report
STUDY_LIST="./study_ids/RNA_seq_ids.txt"
TRACKER_DIR="${BIG_STORAGE}/rnaseq_data/file_tracker/"
SAMPLE_MAP="${BIG_STORAGE}/final_data/rnaseq_processed/Salmon_RNAseq_Combined_TPM_sample_map.csv"
OUTPUT_DIR="./rnaseq_status_report/"

# ──────────────────────────────────────────────────────────────────────────────
# Make sure the log directory exists before the job starts writing to it
mkdir -p ./logs_slurm 2>/dev/null || true

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
  echo "Available .sif files under ${BIG_STORAGE}:"
  find "${BIG_STORAGE}" -name "*.sif" 2>/dev/null
  exit 1
fi

# Sanity-check the tracker dir and study list before burning a job slot
if [ ! -d "${TRACKER_DIR}" ]; then
  echo "ERROR: tracker dir not found at ${TRACKER_DIR}"
  exit 1
fi
if [ ! -f "${STUDY_LIST}" ]; then
  echo "ERROR: study list not found at ${STUDY_LIST}"
  exit 1
fi

mkdir -p "${OUTPUT_DIR}"

echo "+++++++++++++++++ NOW RUNNING RNA-SEQ STATUS REPORT +++++++++++++++++"
apptainer exec \
  --bind "$BIG_STORAGE":"$BIG_STORAGE" \
  --bind "$REPO_DIR":/repo \
  --pwd /repo \
  "$SIF" \
  python /repo/src/scripts/rnaseq_processing_status_report.py \
    --study-list "${STUDY_LIST}" \
    --tracker-dir "${TRACKER_DIR}" \
    --sample-map "${SAMPLE_MAP}" \
    --output-dir "${OUTPUT_DIR}"

EXIT_CODE=$?

echo "========================================"
echo "Finished    : $(date)"
echo "Exit code   : ${EXIT_CODE}"
echo "Report dir  : ${OUTPUT_DIR}"
echo "========================================"
exit ${EXIT_CODE}