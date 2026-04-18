#!/bin/bash
#SBATCH --account=testusers
#SBATCH --partition=all
#SBATCH --time=36:00:00
#SBATCH --mem=16G  
#SBATCH --cpus-per-task=1
#SBATCH --array=16-99%5
#SBATCH --output=/home/alexdelossanto/linuxhome/Dataset_fusion_Microarray/logs_slurm/stdout-%x-%A_%a.txt

# ══════════════════════════════════════════════════
#  PATHS
# ══════════════════════════════════════════════════
BIG_STORAGE="/tudelft.net/staff-umbrella/GeneExpressionStorage"
SIF="$BIG_STORAGE/dataset_fusion.sif"

echo "=================================================="
echo "Job ID:          $SLURM_JOB_ID"
echo "Array Task ID:   $SLURM_ARRAY_TASK_ID"
echo "Node:            $(hostname)"
echo "Date: $(date)"
echo "SIF:             ${SIF}"
echo "=================================================="

# ══════════════════════════════════════════════════
#  NEXTFLOW + JAVA CONFIGURATION
# ══════════════════════════════════════════════════
export JAVA_HOME="$BIG_STORAGE/java"
export NXF_HOME="$BIG_STORAGE/.nextflow"
export NXF_WORK="$BIG_STORAGE/work/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$NXF_WORK"

# ══════════════════════════════════════════════════
#  SINGULARITY CACHES (RAM DISK + BIG STORAGE)
# ══════════════════════════════════════════════════
export NXF_SINGULARITY_CACHEDIR="$BIG_STORAGE/singularity_cache"
mkdir -p "$NXF_SINGULARITY_CACHEDIR"

export APPTAINER_TMPDIR="/dev/shm/apptainer_tmp_${USER}_${SLURM_JOB_ID}"
export SINGULARITY_TMPDIR="$APPTAINER_TMPDIR"
export APPTAINER_CACHEDIR="/dev/shm/apptainer_cache_${USER}_${SLURM_JOB_ID}"
export SINGULARITY_CACHEDIR="$APPTAINER_CACHEDIR"
mkdir -p "$APPTAINER_TMPDIR" "$APPTAINER_CACHEDIR"

# ══════════════════════════════════════════════════
#  CREATE HYBRID PASSWD/GROUP FILES FOR CONTAINER
# ══════════════════════════════════════════════════
MY_PASSWD="$BIG_STORAGE/tmp_passwd_$SLURM_JOB_ID"
MY_GROUP="$BIG_STORAGE/tmp_group_$SLURM_JOB_ID"

cat /etc/passwd > "$MY_PASSWD"
getent passwd $(whoami) >> "$MY_PASSWD"

cat /etc/group > "$MY_GROUP"
for g in $(id -G -n); do
    getent group "$g" >> "$MY_GROUP" || true
done
sort -u "$MY_GROUP" -o "$MY_GROUP"

# ══════════════════════════════════════════════════
#  CREATE WRAPPER TO BYPASS NFS CHMOD BUG
# ══════════════════════════════════════════════════
FAKE_BIN="$BIG_STORAGE/fake_bin_$SLURM_JOB_ID"
mkdir -p "$FAKE_BIN"

cat << 'EOF' > "$FAKE_BIN/singularity"
#!/bin/bash
if [[ "$1" == "pull" ]]; then
    # Nextflow runs `singularity pull` directly on the NFS drive, causing chmod errors.
    # We intercept it, run the pull in the High-Speed RAM disk, and move the result back.
    TEMP_DIR=$(mktemp -d -p /dev/shm)
    ORIG_DIR="$PWD"
    cd "$TEMP_DIR"
    
    /usr/bin/apptainer "$@"
    EXIT_CODE=$?
    
    if [[ $EXIT_CODE -eq 0 ]]; then
        # Copy ALL downloaded files (including .pulling temporary files)
        find . -maxdepth 1 -type f -exec cp {} "$ORIG_DIR/" \;
    fi
    
    cd /tmp
    rm -rf "$TEMP_DIR"
    exit $EXIT_CODE
else
    # For all other commands, pass through directly
    exec /usr/bin/apptainer "$@"
fi
EOF
chmod +x "$FAKE_BIN/singularity"

# ══════════════════════════════════════════════════
#  RUN THE PIPELINE INSIDE THE CONTAINER
# ══════════════════════════════════════════════════
REPO_DIR="$(cd "$(dirname "$0")/.." && pwd)"

apptainer exec \
    --bind /usr/bin/apptainer:/usr/bin/apptainer \
    --bind "$FAKE_BIN/singularity":/usr/bin/singularity \
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
    --env SLURM_JOB_ID="$SLURM_JOB_ID" \
    --env SLURM_ARRAY_TASK_ID="$SLURM_ARRAY_TASK_ID" \
    --env PATH="$BIG_STORAGE/java/bin:/usr/local/bin:/usr/bin:/bin" \
    --env SLURM_ACCOUNT="testusers" \
    --pwd /repo \
    "$SIF" \
    python3 -u /repo/src/data_importing/new_download.py \
        -o "$BIG_STORAGE/" \
        --array_index "$SLURM_ARRAY_TASK_ID" \
        --batch_size 10 --rna --container

# ══════════════════════════════════════════════════
#  CLEANUP
# ══════════════════════════════════════════════════
echo "Cleaning up all temporary directories..."
rm -f "$MY_PASSWD" "$MY_GROUP"
rm -rf "$NXF_WORK"
rm -rf "$FAKE_BIN"
rm -rf "$APPTAINER_TMPDIR"
rm -rf "$APPTAINER_CACHEDIR"

echo "=================================================="
echo "Job Completed"
echo "=================================================="
echo "Job ID:          $SLURM_JOB_ID"
echo "Array Task ID:   $SLURM_ARRAY_TASK_ID"
echo "Node:            $(hostname)"
echo "Date: $(date)"
echo "SIF:             ${SIF}"
echo "=================================================="