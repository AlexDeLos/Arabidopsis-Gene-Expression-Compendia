#!/bin/bash
#SBATCH --job-name=build_salmon_indices
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --account=ewi-insy-prb
#SBATCH --partition=general
#SBATCH --output=/home/nfs/alexdelossanto/Dataset_fusion_Microarray/logs_slurm/build_indices_%j.out
#SBATCH --error=/home/nfs/alexdelossanto/Dataset_fusion_Microarray/logs_slurm/build_indices_%j.err

set -euo pipefail

# --- Environment: use conda tools directly, no Singularity needed ---
module purge
module use /opt/insy/modulefiles
module load miniconda/3.9
eval "$(conda shell.bash hook)"
conda activate download_env

# EDirect for Ler URL resolution
export PATH="${HOME}/edirect:${PATH}"

# Verify all required tools are present before doing any work
echo "Verifying tools..."
for tool in gffread salmon wget; do
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

echo "============================================="
echo " Arabidopsis Salmon Index Builder"
echo " Start: $(date)"
echo "============================================="

# =============================================================================
# HELPER: build a decoy-aware Salmon index for one ecotype
# Args: $1=ecotype  $2=genome_fasta  $3=gtf  $4=index_output_dir
# =============================================================================
build_index() {
    local ECOTYPE="$1"
    local GENOME_FA="$2"
    local GTF="$3"
    local INDEX_DIR="$4"
    local WORK_DIR
    WORK_DIR="$(dirname "$INDEX_DIR")"

    local TRANSCRIPT_FA="${WORK_DIR}/${ECOTYPE}_transcripts.fasta"
    local DECOY_TXT="${WORK_DIR}/${ECOTYPE}_decoys.txt"
    local GENTROME="${WORK_DIR}/${ECOTYPE}_gentrome.fasta"

    echo ""
    echo "--- Building index for: $ECOTYPE ---"
    echo "    Genome:    $GENOME_FA"
    echo "    GTF:       $GTF"
    echo "    Index dir: $INDEX_DIR"

    if [ -d "$INDEX_DIR" ] && [ "$(ls -A "$INDEX_DIR")" ]; then
        echo "Index already exists and is non-empty — skipping."
        return 0
    fi

    mkdir -p "$INDEX_DIR"

    # Step 1: Extract transcripts
    echo "[1/4] Extracting transcripts with gffread..."
    gffread "$GTF" -g "$GENOME_FA" -w "$TRANSCRIPT_FA"
    if [ ! -s "$TRANSCRIPT_FA" ]; then
        echo "ERROR: gffread produced an empty transcript FASTA. Check GTF/FASTA compatibility."
        exit 1
    fi
    echo "      $(grep -c '^>' "$TRANSCRIPT_FA") transcripts extracted."

    # Step 2: Build decoy list from genome chromosome names
    echo "[2/4] Building decoy list..."
    grep "^>" "$GENOME_FA" | cut -d " " -f 1 | tr -d ">" > "$DECOY_TXT"
    echo "      $(wc -l < "$DECOY_TXT") decoy sequences."

    # Step 3: Concatenate transcripts + genome into gentrome
    echo "[3/4] Concatenating gentrome..."
    cat "$TRANSCRIPT_FA" "$GENOME_FA" > "$GENTROME"

    # Step 4: Build index
    # --gencode intentionally omitted: TAIR/custom GTF format, not GENCODE
    echo "[4/4] Running salmon index (8 threads)..."
    salmon index \
        --threads 8 \
        --transcripts "$GENTROME" \
        --decoys "$DECOY_TXT" \
        --index "$INDEX_DIR"

    # Clean up large intermediates, keep decoy list for reference
    rm -f "$GENTROME" "$TRANSCRIPT_FA"

    echo "Index for '$ECOTYPE' complete."
    echo "Contents of $INDEX_DIR:"
    ls -lh "$INDEX_DIR"
}

# =============================================================================
# SECTION 1: Col-0 index
# =============================================================================
echo ""
echo "=== SECTION 1: Col-0 Salmon index ==="

COL0_DIR="$REF_BASE/col-0"

[ -f "$COL0_DIR/col-0.fasta" ] || { echo "ERROR: $COL0_DIR/col-0.fasta not found."; exit 1; }
[ -f "$COL0_DIR/col-0.gtf"   ] || { echo "ERROR: $COL0_DIR/col-0.gtf not found.";   exit 1; }

build_index "col-0" \
    "$COL0_DIR/col-0.fasta" \
    "$COL0_DIR/col-0.gtf" \
    "$COL0_DIR/salmon_index"

# =============================================================================
# SECTION 2: Download Ler-0 genome
# =============================================================================
echo ""
echo "=== SECTION 2: Ler-0 genome download ==="

LER_DIR="$REF_BASE/ler"
mkdir -p "$LER_DIR"

if [ ! -f "$LER_DIR/ler.fasta" ]; then
    echo "Resolving Ler-0 FTP path via NCBI Entrez Direct..."
    LER_FTP_DIR=$(esearch -db assembly -query "GCA_001651475.1" \
        | esummary \
        | xtract -pattern DocumentSummary -element FtpPath_GenBank)

    if [ -z "$LER_FTP_DIR" ] || [ "$LER_FTP_DIR" = "na" ]; then
        echo "ERROR: Could not resolve FTP path for GCA_001651475.1."
        echo "       Test manually: esearch -db assembly -query GCA_001651475.1 | esummary | xtract -pattern DocumentSummary -element FtpPath_GenBank"
        exit 1
    fi

    LER_ASSEMBLY_NAME=$(basename "$LER_FTP_DIR")
    LER_FASTA_URL="${LER_FTP_DIR}/${LER_ASSEMBLY_NAME}_genomic.fna.gz"
    echo "Resolved URL: $LER_FASTA_URL"

    echo "Downloading (~130MB)..."
    wget --progress=bar:force -O "$LER_DIR/ler.fasta.gz" "$LER_FASTA_URL"
    gunzip "$LER_DIR/ler.fasta.gz"
    echo "Ler-0 FASTA ready."
else
    echo "Ler-0 FASTA already exists — skipping download."
fi

# No Ler-specific GTF available in Ensembl Plants yet.
# Col-0 GTF is used as fallback — gene coordinates are highly conserved across ecotypes.
if [ ! -f "$LER_DIR/ler.gtf" ]; then
    echo "No Ler-specific GTF available. Symlinking Col-0 GTF as fallback."
    echo "WARNING: Mapping bias possible in polymorphic regions."
    ln -s "$COL0_DIR/col-0.gtf" "$LER_DIR/ler.gtf"
fi

# =============================================================================
# SECTION 3: Ler index
# =============================================================================
echo ""
echo "=== SECTION 3: Ler Salmon index ==="

build_index "ler" \
    "$LER_DIR/ler.fasta" \
    "$LER_DIR/ler.gtf" \
    "$LER_DIR/salmon_index"

# =============================================================================
# Done
# =============================================================================
echo ""
echo "============================================="
echo " All indices complete."
echo " End: $(date)"
echo "============================================="
echo ""
echo "Final structure:"
find "$REF_BASE" -maxdepth 3 | sort
