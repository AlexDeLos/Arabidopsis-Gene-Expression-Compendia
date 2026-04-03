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
#  SINGULARITY / APPTAINER CACHING (SEPARATED)
# ══════════════════════════════════════════════════
# 1. Nextflow's Image Cache: Must be on shared network storage ($BIG_STORAGE) 
#    so all SLURM nodes can find the final pulled images.
export NXF_SINGULARITY_CACHEDIR="$BIG_STORAGE/singularity_cache"
mkdir -p "$NXF_SINGULARITY_CACHEDIR"

# 2. Apptainer's Internal Cache & Tmp: Must be on local compute node (/tmp) 
#    because Apptainer strictly enforces chmod 0700 during layer downloads, 
#    which the NFS network drive blocks.
export APPTAINER_TMPDIR="/tmp/apptainer_tmp_${USER}_${SLURM_JOB_ID}"
export SINGULARITY_TMPDIR="$APPTAINER_TMPDIR"
export APPTAINER_CACHEDIR="/tmp/apptainer_cache_${USER}_${SLURM_JOB_ID}"
export SINGULARITY_CACHEDIR="$APPTAINER_CACHEDIR"
mkdir -p "$APPTAINER_TMPDIR" "$APPTAINER_CACHEDIR"

# ══════════════════════════════════════════════════
#  RUN PYTHON ORCHESTRATION INSIDE THE CONTAINER
# ══════════════════════════════════════════════════
REPO_DIR="$(cd "$(dirname "$0")/.." && pwd)"

# ══════════════════════════════════════════════════
#  CREATE HYBRID PASSWD/GROUP FILES FOR CONTAINER
# ══════════════════════════════════════════════════
MY_PASSWD="$BIG_STORAGE/tmp_passwd_$SLURM_JOB_ID"
MY_GROUP="$BIG_STORAGE/tmp_group_$SLURM_JOB_ID"

# Copy the system file, then append your specific network user info!
cat /etc/passwd > "$MY_PASSWD"
getent passwd $(whoami) >> "$MY_PASSWD"

cat /etc/group > "$MY_GROUP"
getent group $(id -g -n) >> "$MY_GROUP"

# ══════════════════════════════════════════════════
#  RUN THE PIPELINE
# ══════════════════════════════════════════════════
apptainer exec \
    --bind /usr/bin/apptainer:/usr/bin/apptainer \
    --bind /usr/bin/singularity:/usr/bin/singularity \
    --bind /etc/apptainer:/etc/apptainer \
    --bind /lib64/libsubid.so.3:/usr/lib64/libsubid.so.3 \
    --bind /lib64/libcrypt.so.2:/usr/lib64/libcrypt.so.2 \
    --bind /usr/bin/sbatch:/usr/bin/sbatch \
    --bind /usr/bin/squeue:/usr/bin/squeue \
    --bind /usr/bin/scancel:/usr/bin/scancel \
    --bind /usr/bin/scontrol:/usr/bin/scontrol \
    --bind /usr/lib64/slurm:/usr/lib64/slurm \
    --bind /usr/lib64/libslurm.so.43:/usr/lib64/libslurm.so.43 \
    --bind /usr/lib64/libmunge.so.2:/usr/lib64/libmunge.so.2 \
    --bind /usr/lib64/libmunge.so.2.0.0:/usr/lib64/libmunge.so.2.0.0 \
    --bind /usr/lib64/libauksapi.so.0:/usr/lib64/libauksapi.so.0 \
    --bind /var/run/munge:/var/run/munge \
    --bind /etc/slurm:/etc/slurm \
    --bind /etc/auks:/etc/auks \
    --bind "$MY_PASSWD":/etc/passwd \
    --bind "$MY_GROUP":/etc/group \
    --bind "$BIG_STORAGE/bin/nextflow":/usr/local/bin/nextflow \
    --bind "$BIG_STORAGE/java":"$BIG_STORAGE/java" \
    --bind /tudelft.net/staff-umbrella:/tudelft.net/staff-umbrella \
    --bind /home/alexdelossanto/linuxhome/Dataset_fusion_Microarray:/repo \
    --env LD_LIBRARY_PATH="/usr/lib64:/usr/lib64/slurm" \
    --env JAVA_HOME="$BIG_STORAGE/java" \
    --env NXF_HOME="$NXF_HOME" \
    --env NXF_WORK="$NXF_WORK" \
    --env NXF_SINGULARITY_CACHEDIR="$NXF_SINGULARITY_CACHEDIR" \
    --env APPTAINER_CACHEDIR="$APPTAINER_CACHEDIR" \
    --env SINGULARITY_CACHEDIR="$SINGULARITY_CACHEDIR" \
    --env APPTAINER_TMPDIR="$APPTAINER_TMPDIR" \
    --env SINGULARITY_TMPDIR="$SINGULARITY_TMPDIR" \
    --env APPTAINER_TMPDIR="$APPTAINER_TMPDIR" \
    --env SINGULARITY_TMPDIR="$SINGULARITY_TMPDIR" \
    --env SLURM_JOB_ID="$SLURM_JOB_ID" \
    --env SLURM_ARRAY_TASK_ID="$SLURM_ARRAY_TASK_ID" \
    --env PATH="$BIG_STORAGE/java/bin:/usr/local/bin:/usr/bin:/bin" \
    --pwd /repo \
    "$SIF" \
    python3 -u /repo/src/data_importing/new_download.py \
        -o "$BIG_STORAGE/" \
        --array_index "$SLURM_ARRAY_TASK_ID" \
        --batch_size 1 --rna --container

# Clean up the temporary files after the job finishes
rm -f "$MY_PASSWD" "$MY_GROUP"
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