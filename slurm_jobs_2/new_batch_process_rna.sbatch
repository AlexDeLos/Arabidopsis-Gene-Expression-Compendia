#!/bin/bash
#SBATCH --account=ewi-insy-prb
#SBATCH --partition=general
#SBATCH --time=10:00:00
#SBATCH --qos=medium
#SBATCH --mem=4G
#SBATCH --tmp=800000M
#SBATCH --cpus-per-task=1
#SBATCH --array=2%1
#SBATCH --output=./logs_slurm/stdout-%x-%A_%a.txt

# ══════════════════════════════════════════════════
#  PATHS
# ══════════════════════════════════════════════════
BIG_STORAGE="/tudelft.net/staff-umbrella/GeneExpressionStorage"
HOME="/home/nfs/alexdelossanto/"
TMP="/tmp/${USER}/${SLURM_JOBID}"
export TMP
mkdir -p "$TMP"
SIF="/home/nfs/alexdelossanto/Dataset_fusion_Microarray/dataset_fusion.sif"

# ══════════════════════════════════════════════════
#  /tmp KILL-SWITCH THRESHOLD
#  If the local scratch dir ($TMP) grows past this many
#  GB, the background monitor below scancels the job.
# ══════════════════════════════════════════════════

echo "=================================================="
echo "Job ID:          $SLURM_JOB_ID"
echo "Array Task ID:   $SLURM_ARRAY_TASK_ID"
echo "Node:            $(hostname)"
echo "Date:            $(date)"
echo "SIF:             ${SIF}"
echo "=================================================="

# ══════════════════════════════════════════════════
#  NEXTFLOW + JAVA CONFIGURATION
# ══════════════════════════════════════════════════
export JAVA_HOME="$HOME/java"
export NXF_HOME="$HOME/.nextflow"
export NXF_WORK="$TMP/work_alexdelossatos_nextflow/${SLURM_ARRAY_TASK_ID}"

mkdir -p "$NXF_WORK"

# to NXF_WORK on $BIG_STORAGE.
export TMPDIR="$TMP/scratch"
mkdir -p "$TMPDIR"

# ══════════════════════════════════════════════════
#  SINGULARITY CACHES (RAM DISK + BIG STORAGE)
# ══════════════════════════════════════════════════
export NXF_SINGULARITY_CACHEDIR="$HOME/singularity_cache"
mkdir -p "$NXF_SINGULARITY_CACHEDIR"

export APPTAINER_TMPDIR="/dev/shm/apptainer_tmp_${USER}_${SLURM_JOB_ID}"
export SINGULARITY_TMPDIR="$APPTAINER_TMPDIR"
export APPTAINER_CACHEDIR="/dev/shm/apptainer_cache_${USER}_${SLURM_JOB_ID}"
export SINGULARITY_CACHEDIR="$APPTAINER_CACHEDIR"
mkdir -p "$APPTAINER_TMPDIR" "$APPTAINER_CACHEDIR"

# ══════════════════════════════════════════════════
#  CREATE HYBRID PASSWD/GROUP FILES FOR CONTAINER
# ══════════════════════════════════════════════════
MY_PASSWD="$TMP/tmp_passwd_$SLURM_JOB_ID"
MY_GROUP="$TMP/tmp_group_$SLURM_JOB_ID"

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
FAKE_BIN="$TMP/fake_bin_$SLURM_JOB_ID"
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
#  DYNAMIC SLURM & PLUGIN CONFIGURATION
# ══════════════════════════════════════════════════
REPO_DIR="$(cd "$(dirname "$0")/.." && pwd)"

# 1. Resolve the true physical directory for Slurm config
HOST_SLURM_CONF=$(realpath "${SLURM_CONF:-/etc/slurm/slurm.conf}")
HOST_SLURM_DIR=$(dirname "$HOST_SLURM_CONF")

# 2. Dynamically locate the AUKS plugin to avoid mounting all of /usr/lib64
AUKS_SO=$(find /usr/lib64 -name "auks.so" 2>/dev/null | head -n 1)
if [ -n "$AUKS_SO" ]; then
    AUKS_BIND="--bind $AUKS_SO:$AUKS_SO"
    AUKS_DIR=$(dirname "$AUKS_SO")
    DYNAMIC_LD_LIB="/usr/lib64/slurm:$AUKS_DIR"
else
    AUKS_BIND=""
    DYNAMIC_LD_LIB="/usr/lib64/slurm"
fi

function clean_up {
# ══════════════════════════════════════════════════
#  CLEANUP
# ══════════════════════════════════════════════════
echo "Cleaning up all temporary directories..."
    kill "$MONITOR_PID" 2>/dev/null || true
    rm -f "$MY_PASSWD" "$MY_GROUP"
    rm -rf "$FAKE_BIN"
    rm -rf "$APPTAINER_TMPDIR"
    rm -rf "$APPTAINER_CACHEDIR"
    rm --recursive --force "$TMP" && echo "Clean up of $TMP completed successfully."
    exit
}
trap 'clean_up' EXIT # Setup clean_up to run on exit

# ══════════════════════════════════════════════════
#  /tmp KILL SWITCH
# ══════════════════════════════════════════════════
monitor_tmp_usage() {
    TMP_FS_LIMIT=90          # kill if /tmp filesystem >90% full
    TMP_JOB_LIMIT_GB=700     # kill if this job uses >700 GB

    while true; do
        sleep 60

        # Filesystem usage
        TMP_PCT=$(df -P /tmp | awk 'NR==2 {gsub("%","",$5); print $5}')
        TMP_PCT=${TMP_PCT:-0}

        # This job's usage
        USAGE_KB=$(du -sk "$TMP" 2>/dev/null | cut -f1)
        USAGE_KB=${USAGE_KB:-0}
        USAGE_GB=$((USAGE_KB / 1024 / 1024))

        if [ "$TMP_PCT" -ge "$TMP_FS_LIMIT" ]; then
            echo "$(date) KILL SWITCH: /tmp filesystem is ${TMP_PCT}% full. Cancelling job ${SLURM_JOB_ID}" >&2

            scancel "$SLURM_JOB_ID"

            sleep 5

            pkill -TERM -P $$ 2>/dev/null || true

            exit 1
        fi

        if [ "$USAGE_GB" -ge "$TMP_JOB_LIMIT_GB" ]; then
            echo "$(date) KILL SWITCH: job scratch usage is ${USAGE_GB}GB (limit ${TMP_JOB_LIMIT_GB}GB). Cancelling job ${SLURM_JOB_ID}" >&2

            scancel "$SLURM_JOB_ID"

            sleep 5

            pkill -TERM -P $$ 2>/dev/null || true

            exit 1
        fi
    done
}
monitor_tmp_usage &
MONITOR_PID=$!

# ══════════════════════════════════════════════════
#  RUN THE PIPELINE INSIDE THE CONTAINER
# ══════════════════════════════════════════════════
apptainer exec \
    --bind /usr/bin/apptainer:/usr/bin/apptainer \
    --bind "$FAKE_BIN/singularity":/usr/bin/singularity \
    --bind /etc/apptainer:/etc/apptainer \
    --bind /usr/bin/sbatch:/usr/bin/sbatch \
    --bind /usr/bin/squeue:/usr/bin/squeue \
    --bind /usr/bin/scancel:/usr/bin/scancel \
    --bind /usr/bin/scontrol:/usr/bin/scontrol \
    --bind "$HOST_SLURM_DIR":/etc/slurm \
    --bind /var/run/munge:/var/run/munge \
    --bind /etc/munge:/etc/munge \
    --bind /usr/lib64/slurm:/usr/lib64/slurm \
    --bind /etc/resolv.conf:/etc/resolv.conf \
    --bind /etc/hosts:/etc/hosts \
    --bind /etc/auks:/etc/auks \
    --bind /usr/lib64/libmunge.so.2:/usr/lib64/libmunge.so.2 \
    --bind /usr/lib64/libmunge.so.2.0.0:/usr/lib64/libmunge.so.2.0.0 \
    $AUKS_BIND \
    --bind "$MY_PASSWD":/etc/passwd \
    --bind "$MY_GROUP":/etc/group \
    --bind "$HOME/bin/nextflow":/usr/local/bin/nextflow \
    --bind "$HOME/java":"$HOME/java" \
    --bind /tudelft.net/staff-umbrella:/tudelft.net/staff-umbrella \
    --bind /home/nfs/alexdelossanto/Dataset_fusion_Microarray:/repo \
    --env LD_LIBRARY_PATH="/usr/lib64:/usr/lib64/slurm" \
    --env SLURM_CONF="/etc/slurm/slurm.conf" \
    --env JAVA_HOME="$HOME/java" \
    --env NXF_HOME="$NXF_HOME" \
    --env NXF_WORK="$NXF_WORK" \
    --env NXF_SINGULARITY_CACHEDIR="$NXF_SINGULARITY_CACHEDIR" \
    --env APPTAINER_CACHEDIR="$APPTAINER_CACHEDIR" \
    --env SINGULARITY_CACHEDIR="$SINGULARITY_CACHEDIR" \
    --env APPTAINER_TMPDIR="$APPTAINER_TMPDIR" \
    --env SINGULARITY_TMPDIR="$SINGULARITY_TMPDIR" \
    --env TMPDIR="$TMPDIR" \
    --env SLURM_JOB_ID="$SLURM_JOB_ID" \
    --env SLURM_ARRAY_TASK_ID="$SLURM_ARRAY_TASK_ID" \
    --env PATH="$HOME/java/bin:/usr/local/bin:/usr/bin:/bin" \
    --env SLURM_ACCOUNT="ewi-insy-prb" \
    --pwd /repo \
    "$SIF" \
    python3 -u /repo/src/data_importing/new_download.py \
        -o "$BIG_STORAGE/" \
        --array_index "$SLURM_ARRAY_TASK_ID" \
        --batch_size 1 --rna --container --old



echo "=================================================="
echo "Job Completed"
echo "=================================================="
echo "Job ID:          $SLURM_JOB_ID"
echo "Array Task ID:   $SLURM_ARRAY_TASK_ID"
echo "Node:            $(hostname)"
echo "Date: $(date)"
echo "SIF:             ${SIF}"
echo "=================================================="