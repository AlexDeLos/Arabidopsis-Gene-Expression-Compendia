#!/bin/bash
#SBATCH --account=ewi-insy-prb
#SBATCH --partition=general
#SBATCH --time=10:00:00
#SBATCH --qos=medium
#SBATCH --mem=100G
#SBATCH --cpus-per-task=1
#SBATCH --array=2%1
#SBATCH --output=./logs_slurm/stdout-%x-%A_%a.txt

# ══════════════════════════════════════════════════
#  PATHS
# ══════════════════════════════════════════════════
BIG_STORAGE="/tudelft.net/staff-umbrella/GeneExpressionStorage"
HOME="/home/nfs/alexdelossanto/"
TMP="/tmp/${USER}/${SLURM_JOB_ID}"
export TMP
mkdir -p "$TMP"
SIF="/home/nfs/alexdelossanto/Dataset_fusion_Microarray/dataset_fusion.sif"


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
export NXF_WORK="$BIG_STORAGE/work_alexdelossatos_nextflow/${SLURM_ARRAY_TASK_ID}"

mkdir -p "$NXF_WORK"

# ══════════════════════════════════════════════════
#  SINGULARITY CACHES (RAM DISK + BIG STORAGE)
# ══════════════════════════════════════════════════
# export NXF_SINGULARITY_CACHEDIR="$HOME/singularity_cache"
# mkdir -p "$NXF_SINGULARITY_CACHEDIR"
export NXF_SINGULARITY_CACHEDIR="/tudelft.net/staff-umbrella/GeneExpressionStorage/singularity_cache"
mkdir -p "$NXF_SINGULARITY_CACHEDIR"

export APPTAINER_TMPDIR="/dev/shm/apptainer_tmp_${USER}_${SLURM_JOB_ID}"
export SINGULARITY_TMPDIR="$APPTAINER_TMPDIR"
export APPTAINER_CACHEDIR="/dev/shm/apptainer_cache_${USER}_${SLURM_JOB_ID}"
export SINGULARITY_CACHEDIR="$APPTAINER_CACHEDIR"
mkdir -p "$APPTAINER_TMPDIR" "$APPTAINER_CACHEDIR"

# ══════════════════════════════════════════════════
#  CREATE HYBRID PASSWD/GROUP FILES FOR CONTAINER
# ══════════════════════════════════════════════════
MY_PASSWD="$TMP/passwd"
MY_GROUP="$TMP/group"

cat /etc/passwd > "$MY_PASSWD"
getent passwd $USER >> "$MY_PASSWD"

cat /etc/group > "$MY_GROUP"
getent group "$(id -gn)" >> "$MY_GROUP"
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
        shopt -s dotglob
        mv ./* "$ORIG_DIR/" 2>/dev/null || true
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
#  CLEANUP
# ══════════════════════════════════════════════════

function clean_up {
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
    while true; do
        sleep 60

        TMP_PCT=$(df -P /tmp 2>/dev/null | awk 'NR==2 {gsub("%","",$5); print $5}')
        if [ -z "$TMP_PCT" ]; then
            echo "$(date) WARNING: failed to determine /tmp usage" >&2
            continue
        fi

        if [ "$TMP_PCT" -ge 90 ]; then
            echo "$(date) KILL SWITCH: /tmp filesystem is ${TMP_PCT}% full. Cancelling job ${SLURM_JOB_ID}" >&2
            scancel "$SLURM_JOB_ID"
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
    --bind /etc/slurm:/etc/slurm \
    --bind /var/run/munge:/var/run/munge \
    --bind /etc/munge:/etc/munge \
    --bind /usr/lib64/slurm:/usr/lib64/slurm \
    --bind /etc/resolv.conf:/etc/resolv.conf \
    --bind /etc/hosts:/etc/hosts \
    --bind /etc/auks:/etc/auks \
    --bind /usr/lib64/libmunge.so.2:/usr/lib64/libmunge.so.2 \
    --bind /usr/lib64/libmunge.so.2.0.0:/usr/lib64/libmunge.so.2.0.0 \
    --bind "$MY_PASSWD":/etc/passwd \
    --bind "$MY_GROUP":/etc/group \
    --bind "$HOME/bin/nextflow":/usr/local/bin/nextflow \
    --bind "$HOME/java":"$HOME/java" \
    --bind /tudelft.net/staff-umbrella/GeneExpressionStorage:/tudelft.net/staff-umbrella/GeneExpressionStorage \
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