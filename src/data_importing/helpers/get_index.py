import os
import subprocess
import urllib.request
import shutil

# CONFIGURATION
# Where to store the genome
genome_dir = "new_storage/genome_index"
# Ensembl FTP link for Arabidopsis thaliana (TAIR10) DNA
fasta_url = "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-56/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz"
fasta_gz = os.path.join(genome_dir, "arabidopsis.fa.gz")
fasta_file = os.path.join(genome_dir, "arabidopsis.fa")
index_prefix = os.path.join(genome_dir, "tair10")

def setup_genome():
    if not os.path.exists(genome_dir):
        os.makedirs(genome_dir)

    # 1. Download FASTA
    if not os.path.exists(fasta_file):
        print(f"Downloading Genome from Ensembl...")
        try:
            urllib.request.urlretrieve(fasta_url, fasta_gz)
        except Exception as e:
            print(f"Error downloading: {e}")
            return
        
        print("Extracting FASTA...")
        # Unzip the file
        subprocess.run(["gunzip", "-f", fasta_gz], check=True)
    else:
        print("FASTA file already exists. Skipping download.")

    # 2. Build HISAT2 Index
    # Check if index already exists to avoid rebuilding
    if os.path.exists(f"{index_prefix}.1.ht2"):
        print("Index already exists. Skipping build.")
        return

    print("Building HISAT2 Index (This may take 15-30 mins)...")
    # Command: hisat2-build <input_fasta> <output_prefix>
    cmd = ["hisat2-build", "-p", "4", fasta_file, index_prefix]
    
    try:
        subprocess.run(cmd, check=True)
        print("\nSUCCESS! Genome Index built.")
        print(f"Your PATH_TO_INDEX is: {os.path.abspath(index_prefix)}")
    except FileNotFoundError:
        print("ERROR: 'hisat2-build' command not found.")
        print("Please install it via: conda install -c bioconda hisat2")
    except subprocess.CalledProcessError as e:
        print(f"Index build failed: {e}")

if __name__ == "__main__":
    setup_genome()