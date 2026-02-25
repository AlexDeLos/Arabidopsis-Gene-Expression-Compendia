import os
import subprocess
import sys
import gzip
import re

def patch_araport_gff3(input_gff, output_gff):
    """
    Injects an explicit 'gene_id' tag into every GFF3 line.
    This prevents gffread/RSEM from crashing on Araport11 miRNAs and 
    transposons that lack standard parent 'gene' hierarchies.
    """
    if os.path.exists(output_gff):
        print(f"Patched GFF3 already exists at {output_gff}. Skipping patching.")
        return output_gff
        
    print("Patching Araport11 GFF3 to ensure nf-core GTF compatibility...")
    with open(input_gff, 'r') as fin, open(output_gff, 'w') as fout:
        for line in fin:
            if line.startswith('#'):
                fout.write(line)
                continue
                
            parts = line.strip().split('\t')
            if len(parts) == 9:
                attrs = parts[8]
                
                # Deduce a valid gene_id from the Parent or ID tag
                parent_match = re.search(r'Parent=([^;]+)', attrs)
                if parent_match:
                    # e.g., Parent=AT1G01010.1 -> gene_id=AT1G01010
                    gene_id = parent_match.group(1).split('.')[0] 
                else:
                    id_match = re.search(r'ID=([^;]+)', attrs)
                    if id_match:
                        gene_id = id_match.group(1).split('.')[0]
                    else:
                        gene_id = "unassigned_gene"
                        
                # Inject gene_id into the attributes column if it is missing
                if 'gene_id=' not in attrs:
                    parts[8] = f"{attrs};gene_id={gene_id}"
                    
                fout.write('\t'.join(parts) + '\n')
            else:
                fout.write(line)
                
    print("  -> GFF3 successfully patched and saved!")
    return output_gff

def create_dummy_input():
    print("Creating dummy input files to bypass validation...")
    fq1 = "dummy_1.fastq.gz"
    fq2 = "dummy_2.fastq.gz"
    samplesheet = "dummy_samplesheet.csv"
    
    for fq in [fq1, fq2]:
        with gzip.open(fq, "wt") as f:
            f.write("@dummy_read\nACGT\n+\nIIII\n")
            
    with open(samplesheet, "w") as f:
        f.write("sample,fastq_1,fastq_2,strandedness\n")
        f.write(f"dummy_sample,{os.path.abspath(fq1)},{os.path.abspath(fq2)},unstranded\n")
        
    return samplesheet, [fq1, fq2, samplesheet]

def build_genome_indices(fasta_path, original_gff_path, out_dir):
    print("===================================================")
    print("Starting Nextflow Reference Genome Index Build")
    print("===================================================")
    
    os.makedirs(out_dir, exist_ok=True)
    
    # 1. Patch the GFF3 file and save it in the out_dir
    patched_gff_path = os.path.join(out_dir, "Araport11_patched_for_Nextflow.gff")
    patched_gff_path = patch_araport_gff3(original_gff_path, patched_gff_path)
    
    # 2. Generate dummy FASTQ files
    dummy_csv, files_to_cleanup = create_dummy_input()
    
    # 3. Construct the Nextflow command (Using the PATCHED GFF)
    cmd = [
        "nextflow", "run", "nf-core/rnaseq",
        "-profile", "singularity", 
        "-revision", "3.14.0",
        "--input", dummy_csv,          
        "--fasta", fasta_path,
        "--gff", patched_gff_path,     # <-- Handing Nextflow the fixed GFF
        "--save_reference",
        "--outdir", out_dir,
        "--skip_qc",
        "--skip_trimming",             
        "--skip_alignment",            
        "--skip_pseudo_alignment"      
    ]
    
    print(f"Running command:\n{' '.join(cmd)}\n")
    
    try:
        subprocess.run(cmd, check=True)
        print("\n[SUCCESS] Genome indices built successfully!")
        
        final_index_path = os.path.join(out_dir, "reference")
        print(f"Your pre-built indices are now saved in: {final_index_path}")
        
    except subprocess.CalledProcessError as e:
        print(f"\n[ERROR] Nextflow pipeline failed with exit code {e.returncode}")
        sys.exit(1)
    finally:
        for f in files_to_cleanup:
            if os.path.exists(f):
                os.remove(f)

if __name__ == "__main__":
    FASTA_FILE = "/tudelft.net/staff-umbrella/GeneExpressionStorage/files_for_rna_seq/GCA_000001735.2_TAIR10.1_genomic_renamed.fna"
    GFF_FILE = "/tudelft.net/staff-umbrella/GeneExpressionStorage/files_for_rna_seq/Araport11_GFF3_genes_transposons.20250813.gff"
    
    SHARED_INDEX_DIR = "/tudelft.net/staff-umbrella/GeneExpressionStorage/rnaseq_data/genome_indices"
    
    build_genome_indices(FASTA_FILE, GFF_FILE, SHARED_INDEX_DIR)