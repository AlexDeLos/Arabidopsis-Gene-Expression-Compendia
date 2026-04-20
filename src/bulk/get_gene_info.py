"""
build_gene_info.py
Produces arabidopsis_gene_info.csv with columns:
  tair_id, gene_length, in_rnaseq, in_microarray

Gene length = sum of all unique exon base-pairs (union of exons per gene),
which is the standard "effective gene length" used in RNA-seq tools like DESeq2.

Source: Araport11 GTF — same annotation used by the nf-core/rnaseq pipeline.
"""

import gzip
import re
# import numpy as np
import pandas as pd

# ── Paths — adjust as needed ──────────────────────────────────────────────────
RNA_PATH   = './new_storage/final_data/rnaseq_processed/filter.csv'
ARRAY_PATH = './new_storage/final_data/filter.csv'   # samples × genes
GTF_PATH   = './Araport11.gtf.gz'
OUT_PATH   = './src/bulk/metadata/arabidopsis_gene_info.csv'

# ── Step 1: gene universe from your data ─────────────────────────────────────
# RNA-seq matrix is genes × samples (index_col=0 gives gene rows)
# If yours is already samples × genes, remove the .T
rna_genes   = set(pd.read_csv(RNA_PATH,   index_col=0).index.tolist())
array_genes = set(pd.read_csv(ARRAY_PATH, index_col=0).index.tolist())

# Adjust if your matrix orientation differs — check with:
# print(pd.read_csv(RNA_PATH, index_col=0).shape)

union = sorted(rna_genes | array_genes)
print(f"RNA-seq:    {len(rna_genes):,} genes")
print(f"Microarray: {len(array_genes):,} genes")
print(f"Union:      {len(union):,} genes")
print(f"Overlap:    {len(rna_genes & array_genes):,} genes")

# ── Step 2: compute exon-union lengths from Araport11 GTF ────────────────────
# For each gene, collect all exon intervals, merge overlapping ones, sum lengths.
# This matches what tools like DESeq2/tximport use for length normalisation.

print("\nParsing GTF for exon lengths...")
exons = {}   # gene_id -> list of (start, end) tuples, 0-based half-open

open_fn = gzip.open if GTF_PATH.endswith('.gz') else open

with open_fn(GTF_PATH, 'rt') as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.rstrip('\n').split('\t')
        if len(fields) < 9 or fields[2] != 'exon':
            continue
        m = re.search(r'gene_id "([^"]+)"', fields[8])
        if not m:
            continue
        gene_id = m.group(1)
        start, end = int(fields[3]) - 1, int(fields[4])  # GTF is 1-based
        exons.setdefault(gene_id, []).append((start, end))

def sum_merged_intervals(intervals):
    """Merge overlapping intervals and return total length."""
    merged, _total = [], 0
    for s, e in sorted(intervals):
        if merged and s <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append([s, e])
    return sum(e - s for s, e in merged)

gene_lengths = {g: sum_merged_intervals(ivs) for g, ivs in exons.items()}
print(f"Gene lengths computed for {len(gene_lengths):,} genes in GTF")

# ── Step 3: assemble final dataframe ─────────────────────────────────────────
records = []
missing = []
for g in union:
    length = gene_lengths.get(g)
    if length is None:
        missing.append(g)
    records.append({
        'tair_id':        g,
        'gene_length':    length,   # None for genes absent from GTF
        'in_rnaseq':      g in rna_genes,
        'in_microarray':  g in array_genes,
    })

df = pd.DataFrame(records)
print(f"\nGenes with length found: {df['gene_length'].notna().sum():,}")
print(f"Genes missing from GTF:  {len(missing):,}")
if missing:
    print("  First 10:", missing[:10])

df.to_csv(OUT_PATH, index=False)
print(f"\nSaved → {OUT_PATH}")
print(df.head())