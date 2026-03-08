#!/bin/bash
#SBATCH --job-name=build_salmon_index
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --account=ewi-insy-prb
#SBATCH --partition=general

STORAGE="/tudelft.net/staff-umbrella/GeneExpressionStorage"
REF_DIR="$STORAGE/files_for_rna_seq/col-0"
INDEX_DIR="$REF_DIR/salmon_index"

# nf-core/rnaseq uses a specific transcript FASTA as the index input.
# Step 1: Extract transcript sequences from genome + GTF using gffread
module load singularity

singularity exec $SINGULARITY_CACHEDIR/depot.galaxyproject.org-singularity-gffread-0.12.7.sif \
    gffread \
    "$REF_DIR/col-0.fasta" \
    -g "$REF_DIR/col-0.fasta" \
    -w "$REF_DIR/col-0_transcripts.fasta" \
    -C  # CDS only — use without -C for full transcriptome

# Step 2: Build the Salmon index
singularity exec $SINGULARITY_CACHEDIR/depot.galaxyproject.org-singularity-salmon-1.10.1.sif \
    salmon index \
    -t "$REF_DIR/col-0_transcripts.fasta" \
    -d <(grep "^>" "$REF_DIR/col-0.fasta" | cut -d " " -f 1 | tr -d ">") \
    -i "$INDEX_DIR" \
    -p 8