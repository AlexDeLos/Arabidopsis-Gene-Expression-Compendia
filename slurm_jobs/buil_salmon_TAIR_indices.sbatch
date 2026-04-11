#!/bin/bash
#SBATCH --job-name=build_salmon_tair12
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --account=testusers
#SBATCH --partition=all
#SBATCH --output=/tudelft.net/staff-homes-linux/d/alexdelossanto/Dataset_fusion_Microarray/logs_slurm/build_indices_%j.out
#SBATCH --error=/tudelft.net/staff-homes-linux/d/alexdelossanto/Dataset_fusion_Microarray/logs_slurm/build_indices_%j.err

set -eo pipefail

# --- Environment: use conda tools directly ---
module purge
module use /opt/insy/modulefiles
module load miniconda/3.9
eval "$(conda shell.bash hook)"
conda activate download_env

# Verify all required tools are present
echo "Verifying tools..."
for tool in gffread salmon; do
    if ! command -v "$tool" &>/dev/null; then
        echo "ERROR: '$tool' not found in download_env."
        echo "       Install with: conda install -c bioconda $tool"
        exit 1
    fi
    echo "  OK: $tool -> $(which $tool)"
done
echo ""

STORAGE="/tudelft.net/staff-umbrella/GeneExpressionStorage"
REF_BASE="$STORAGE/files_for_rna_seq"
TAIR12_DIR="$REF_BASE/tair12"

echo "============================================="
echo " Arabidopsis TAIR12 Salmon Index Builder"
echo " Start: $(date)"
echo "============================================="

# =============================================================================
# HELPER: build a Salmon index
# =============================================================================
build_index() {
    local ECOTYPE="$1"
    local GENOME_GZ="$2"
    local GFF_GZ="$3"
    local INDEX_DIR="$4"
    local WORK_DIR
    WORK_DIR="$(dirname "$INDEX_DIR")"

    local TRANSCRIPT_FA="${WORK_DIR}/${ECOTYPE}_transcripts.fasta"
    local TMP_FA="${WORK_DIR}/tmp_genome.fasta"
    local TMP_GFF="${WORK_DIR}/tmp_annotation.gff"

    echo ""
    echo "--- Building index for: $ECOTYPE ---"
    echo "    Genome (gz): $GENOME_GZ"
    echo "    GFF (gz):    $GFF_GZ"
    echo "    Index dir:   $INDEX_DIR"

    if [ -d "$INDEX_DIR" ] && [ "$(ls -A "$INDEX_DIR")" ]; then
        echo "Index already exists and is non-empty — skipping."
        return 0
    fi

    mkdir -p "$INDEX_DIR"

    # Step 1: Temporarily decompress files for gffread compatibility
    echo "[1/3] Temporarily uncompressing genome and GFF..."
    gunzip -c "$GENOME_GZ" > "$TMP_FA"
    gunzip -c "$GFF_GZ" > "$TMP_GFF"

    # Step 2: Extract transcripts
    echo "[2/3] Extracting transcripts with gffread..."
    gffread "$TMP_GFF" -g "$TMP_FA" -w "$TRANSCRIPT_FA"
    
    if [ ! -s "$TRANSCRIPT_FA" ]; then
        echo "ERROR: gffread produced an empty transcript FASTA."
        exit 1
    fi
    echo "      $(grep -c '^>' "$TRANSCRIPT_FA") transcripts extracted."

    # Step 3: Build index directly from transcript FASTA
    echo "[3/3] Running salmon index (8 threads)..."
    salmon index \
        --threads 8 \
        --transcripts "$TRANSCRIPT_FA" \
        --index "$INDEX_DIR"

    # Clean up temporary files
    echo "Cleaning up temporary uncompressed files..."
    rm -f "$TRANSCRIPT_FA" "$TMP_FA" "$TMP_GFF"

    echo "Index for '$ECOTYPE' complete."
    echo "Contents of $INDEX_DIR:"
    ls -lh "$INDEX_DIR"
}

# =============================================================================
# SECTION 1: TAIR12 index
# =============================================================================
echo ""
echo "=== SECTION 1: TAIR12 Salmon index ==="

[ -f "$TAIR12_DIR/tair12.fasta.gz" ] || { echo "ERROR: $TAIR12_DIR/tair12.fasta.gz not found."; exit 1; }
[ -f "$TAIR12_DIR/tair12_annotation.gff.gz" ] || { echo "ERROR: $TAIR12_DIR/tair12_annotation.gff.gz not found."; exit 1; }

build_index "tair12" \
    "$TAIR12_DIR/tair12.fasta.gz" \
    "$TAIR12_DIR/tair12_annotation.gff.gz" \
    "$TAIR12_DIR/salmon_index"

# =============================================================================
# Done
# =============================================================================
echo ""
echo "============================================="
echo " All indices complete."
echo " End: $(date)"
echo "============================================="
echo ""
echo "Final structure of TAIR12 directory:"
find "$TAIR12_DIR" -maxdepth 2 | sort