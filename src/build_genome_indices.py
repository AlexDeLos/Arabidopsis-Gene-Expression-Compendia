import os
import subprocess
import sys

def build_genome_indices(fasta_path, gff_path, out_dir):
    print("===================================================")
    print("Starting Nextflow Reference Genome Index Build")
    print("===================================================")
    
    # Ensure output directory exists
    os.makedirs(out_dir, exist_ok=True)
    
    # Construct the Nextflow command
    cmd = [
        "nextflow", "run", "nf-core/rnaseq",
        "-profile", "singularity", 
        "-revision", "3.14.0",
        "--fasta", fasta_path,
        "--gff", gff_path,
        "--save_reference",
        "--outdir", out_dir,
        "--skip_alignment",
        "--skip_pseudo_alignment",
        "--skip_qc"
    ]
    
    print(f"Running command:\n{' '.join(cmd)}\n")
    
    try:
        # Execute the command and stream output directly to the slurm log
        subprocess.run(cmd, check=True)
        print("\n[SUCCESS] Genome indices built successfully!")
        
        # Nextflow puts references in a 'reference' subfolder inside the outdir
        final_index_path = os.path.join(out_dir, "reference")
        print(f"Your pre-built indices are now saved in: {final_index_path}")
        print("You can now point your batch pipeline to this folder!")
        
    except subprocess.CalledProcessError as e:
        print(f"\n[ERROR] Nextflow pipeline failed with exit code {e.returncode}")
        sys.exit(1)

if __name__ == "__main__":
    # --- UPDATE THESE PATHS ---
    # Point these to the shared reference folder you mentioned
    FASTA_FILE = "/tudelft.net/staff-umbrella/GeneExpressionStorage/files_for_rna_seq/GCA_000001735.2_TAIR10.1_genomic_renamed.fna"
    GFF_FILE = "/tudelft.net/staff-umbrella/GeneExpressionStorage/files_for_rna_seq/Araport11_GFF3_genes_transposons.20250813.gff"
    
    # This is the central location where you want the compiled indices saved forever
    SHARED_INDEX_DIR = "/tudelft.net/staff-umbrella/GeneExpressionStorage/rnaseq_data/genome_indices"
    
    build_genome_indices(FASTA_FILE, GFF_FILE, SHARED_INDEX_DIR)