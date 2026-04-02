#!/bin/bash

#SBATCH --partition=general
#SBATCH --qos=medium
#SBATCH --time=36:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --mail-type=END
#SBATCH --job-name=batch_RNA_seq_processing
#SBATCH --output=/home/nfs/alexdelossanto/Dataset_fusion_Microarray/logs_slurm/stdout-%x-%A_%a.txt
#SBATCH --array=0-9%1
#SBATCH --export=NONE

module purge
module use /opt/insy/modulefiles
module load miniconda/3.9

# 2. Properly initialize Conda for non-interactive scripts
eval "$(conda shell.bash hook)"
conda activate download_env

# Print job information
echo "=================================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Node: $(hostname)"
echo "Date: $(date)"
echo "Environment activated: $CONDA_DEFAULT_ENV"
echo "Python being used: $(which python3)"
echo "=================================================="

# Ensure clean Java environment
unset JAVA_HOME
unset NXF_JAVA_HOME

# Define your base storage paths
export BIG_STORAGE="/tudelft.net/staff-umbrella/GeneExpressionStorage"

# =============================================================================
# WORK DIR STRATEGY
# =============================================================================
# NXF_WORK *must* be on shared storage (umbrella) — NOT /tmp.
#
# Reason: Nextflow with SLURM executor submits each process as a separate
# SLURM job. Those jobs can land on ANY compute node. /tmp is node-local:
# a job on node B cannot read files written to /tmp on node A.
# Using /tmp for NXF_WORK causes every sub-job to fail immediately with
# "terminated by external system" and exit status "-" (SLURM kills it
# because the work directory doesn't exist on the assigned node).
#
# Per-job subdirectory prevents parallel array tasks from colliding in work/.
# cleanup = true in .nextflow.config deletes successful task files automatically.
# The Python run_and_delete logic deletes FASTQs after each batch succeeds.
# These two together keep storage from growing uncontrollably.
export NXF_WORK="$BIG_STORAGE/work/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$NXF_WORK"

echo "Nextflow Work Directory: $NXF_WORK"

# =============================================================================
# SINGULARITY / APPTAINER CACHE
# =============================================================================
# Singularity image cache (.img files) -> shared home dir so all nodes reuse
# the same pulled images. Shared across jobs — intentionally NOT per-job.
export NXF_SINGULARITY_CACHEDIR="/home/nfs/alexdelossanto/singularity_cache"
mkdir -p "$NXF_SINGULARITY_CACHEDIR"
export NXF_SINGULARITY_PULLTIMEOUT="120m"

# Apptainer TMPDIR for SIF building -> /tmp is fine here because SIF builds
# happen on the login/submit node, not inside the sub-jobs. The nodev issue
# only matters on compute nodes; if you still hit it, change to:
#   /home/nfs/alexdelossanto/apptainer_tmp
export APPTAINER_TMPDIR="/tmp/apptainer_tmp_${USER}"
export SINGULARITY_TMPDIR="$APPTAINER_TMPDIR"
mkdir -p "$APPTAINER_TMPDIR"

# Layer/blob cache for Singularity pulls -> /tmp is fine (node-local, ephemeral)
export APPTAINER_CACHEDIR="/tmp/apptainer_cache_${USER}"
export SINGULARITY_CACHEDIR="$APPTAINER_CACHEDIR"
mkdir -p "$APPTAINER_CACHEDIR"

# Pass TMPDIR and CACHEDIR into sub-jobs so they inherit the correct values
export SINGULARITYENV_PREPEND_PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
export APPTAINERENV_PREPEND_PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"


BATCH_SIZE=10
echo "Work Directory: $NXF_WORK"
echo "Processing Batch Index: $SLURM_ARRAY_TASK_ID with Size: $BATCH_SIZE"

# --- RUN THE PIPELINE ---
echo "Running RNA SEQ download and processing..."
python3 -u src/data_importing/new_download.py \
    -o "$BIG_STORAGE/" \
    --array_index $SLURM_ARRAY_TASK_ID \
    --batch_size $BATCH_SIZE \
    --rna
PIPELINE_EXIT=$?

# =============================================================================
# POST-RUN CLEANUP
# Remove the per-job work directory after the job finishes.
# cleanup = true in .nextflow.config already deletes files for successful tasks
# during the run; this catches anything left over (failed task dirs, metadata).
# Comment this out if you need to debug a failed run.
# =============================================================================
if [ -d "$NXF_WORK" ]; then
    echo "Cleaning up work directory: $NXF_WORK"
    rm -rf "$NXF_WORK"
fi

echo "Job finished with exit code $PIPELINE_EXIT"
echo "Date: $(date)"
exit $PIPELINE_EXIT