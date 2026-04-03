#!/bin/bash
#SBATCH --time=36:00:00
#SBATCH --mem=6G
#SBATCH --cpus-per-task=1
#SBATCH --array=0-5%10
#SBATCH --output=/home/alexdelossanto/linuxhome/Dataset_fusion_Microarray/logs_slurm/stdout-%x-%A_%a.txt

# ══════════════════════════════════════════════════
#  PATHS
# ══════════════════════════════════════════════════
BIG_STORAGE="/tudelft.net/staff-umbrella/GeneExpressionStorage"
SIF="$BIG_STORAGE/dataset_fusion.sif"

# ══════════════════════════════════════════════════
#  START LOGGING
# ══════════════════════════════════════════════════
echo "=================================================="
echo "Job ID:          $SLURM_JOB_ID"
echo "Array Task ID:   $SLURM_ARRAY_TASK_ID"
echo "Node:            $(hostname)"
echo "Date:            $(date)"
echo "SIF:             ${SIF}"
echo "=================================================="

# ══════════════════════════════════════════════════
#  NEXTFLOW + JAVA (installed in BIG_STORAGE, no conda/module needed)
# ══════════════════════════════════════════════════
export JAVA_HOME="$BIG_STORAGE/java"
export NXF_HOME="$BIG_STORAGE/.nextflow"

# ══════════════════════════════════════════════════
#  NEXTFLOW WORK DIR — must be on shared storage, never /tmp
#  Each array task gets its own subdir to avoid collisions
# ══════════════════════════════════════════════════
export NXF_WORK="$BIG_STORAGE/work/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$NXF_WORK"

# ══════════════════════════════════════════════════
#  SINGULARITY / APPTAINER CACHE
#  Shared across all jobs — never cleaned automatically
# ══════════════════════════════════════════════════
export NXF_SINGULARITY_CACHEDIR="$BIG_STORAGE/singularity_cache"
export APPTAINER_TMPDIR="/tmp/apptainer_tmp_${USER}"
export SINGULARITY_TMPDIR="$APPTAINER_TMPDIR"
mkdir -p "$APPTAINER_TMPDIR"

# ══════════════════════════════════════════════════
#  RUN PYTHON ORCHESTRATION INSIDE THE CONTAINER
#
#  Bind summary (all discovered empirically on DAIC):
#    SLURM binaries    — sbatch/squeue/scancel need to be host versions
#    libslurmfull.so   — SLURM shared lib, not present in container
#    /usr/lib64/slurm  — auth_munge.so + cred_munge.so plugins
#    libmunge.so*      — munge auth library
#    /var/run/munge    — munge daemon socket (live auth)
#    /etc/slurm        — SLURM cluster config files
#    /etc/passwd+group — needed so SLURM can validate SlurmUser
#    nextflow binary   — host nextflow, bound into container PATH
#    java dir          — JDK installed in BIG_STORAGE
#    staff-umbrella    — all data, configs, outputs, SIF itself
# ══════════════════════════════════════════════════
apptainer exec \
    --bind /usr/bin/sbatch:/usr/bin/sbatch \
    --bind /usr/bin/squeue:/usr/bin/squeue \
    --bind /usr/bin/scancel:/usr/bin/scancel \
    --bind /usr/lib64/slurm:/usr/lib64/slurm \
    --bind /usr/lib64/libmunge.so.2:/usr/lib64/libmunge.so.2 \
    --bind /usr/lib64/libmunge.so.2.0.0:/usr/lib64/libmunge.so.2.0.0 \
    --bind /var/run/munge:/var/run/munge \
    --bind /etc/slurm:/etc/slurm \
    --bind /etc/passwd:/etc/passwd \
    --bind /etc/group:/etc/group \
    --bind "$BIG_STORAGE/bin/nextflow":/usr/local/bin/nextflow \
    --bind "$BIG_STORAGE/java":"$BIG_STORAGE/java" \
    --bind /tudelft.net/staff-umbrella:/tudelft.net/staff-umbrella \
    --env JAVA_HOME="$BIG_STORAGE/java" \
    --env NXF_HOME="$NXF_HOME" \
    --env NXF_WORK="$NXF_WORK" \
    --env NXF_SINGULARITY_CACHEDIR="$NXF_SINGULARITY_CACHEDIR" \
    --env APPTAINER_TMPDIR="$APPTAINER_TMPDIR" \
    --env SINGULARITY_TMPDIR="$SINGULARITY_TMPDIR" \
    --env SLURM_JOB_ID="$SLURM_JOB_ID" \
    --env SLURM_ARRAY_TASK_ID="$SLURM_ARRAY_TASK_ID" \
    --env PATH="$BIG_STORAGE/java/bin:/usr/local/bin:/usr/bin:/bin" \
    # --pwd "/home/alexdelossanto/linuxhome/Dataset_fusion_Microarray" \
    "$SIF" \
    python3 -u src/data_importing/new_download.py \
        -o "$BIG_STORAGE/" \
        --array_index "$SLURM_ARRAY_TASK_ID" \
        --batch_size 10 --rna

# ══════════════════════════════════════════════════
#  EPILOG — clean up work dir after job completes
#  (Nextflow cleanup=true removes per-task files during run;
#   this rm removes the top-level work subdir afterwards)
# ══════════════════════════════════════════════════
rm -rf "$NXF_WORK"

# ══════════════════════════════════════════════════
#  END LOGGING
# ══════════════════════════════════════════════════
echo "=================================================="
echo "Job ID:          $SLURM_JOB_ID"
echo "Array Task ID:   $SLURM_ARRAY_TASK_ID"
echo "Node:            $(hostname)"
echo "Date:            $(date)"
echo "SIF:             ${SIF}"
echo "=================================================="
