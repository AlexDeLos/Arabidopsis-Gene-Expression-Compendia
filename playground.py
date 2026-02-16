import GEOparse
import pandas as pd
import numpy as np
import subprocess
import os
import glob
import rpy2
from rpy2.robjects import r, pandas2ri
from rpy2.robjects.packages import importr
import tarfile


from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import StrVector
import rpy2.robjects as robjects

def install_r_dependencies():
    """
    Checks for required Bioconductor packages and installs them if missing.
    """
    print("Checking R dependencies...")
    
    # List of required R packages for this study
    required_pkgs = ["BiocManager", "affy", "genefilter", "hgu133plus2.db", "DESeq2"]
    
    # R utilities
    utils = importr('utils')
    base = importr('base')
    
    # 1. Check/Install BiocManager first (the package manager for BioC)
    print("Installing BiocManager...")
    utils.install_packages("BiocManager", repos="http://cran.us.r-project.org")
    
    # 2. Use BiocManager to install the rest
    bioc_manager = importr('BiocManager')
    
    for pkg in required_pkgs:
        if pkg == "BiocManager": continue # Already handled
        
        if not robjects.r.installed_packages().rx(True, "Package").count(pkg):
            print(f"Installing missing R package: {pkg} ...")
            # Equivalent to R: BiocManager::install("pkg")
            bioc_manager.install(StrVector([pkg]), ask=False)
        else:
            print(f" - {pkg} is already installed.")
            
    print("All R dependencies are ready.\n")

# Activate automatic conversion for rpy2
# pandas2ri.activate()

class StudyPipeline:
    def __init__(self, gse_id, working_dir="./study_data_test"):
        self.gse_id = gse_id
        self.working_dir = working_dir
        self.raw_dir = os.path.join(working_dir, gse_id, "raw")
        self.processed_dir = os.path.join(working_dir, gse_id, "processed")
        
        # Create directories
        os.makedirs(self.raw_dir, exist_ok=True)
        os.makedirs(self.processed_dir, exist_ok=True)
        
        print(f"Loading Metadata for {gse_id}...")
        self.gse = GEOparse.get_GEO(geo=gse_id, destdir=self.working_dir)

    def run_microarray_pipeline(self):
        """
        Replicates Microarray Pipeline with added .tar extraction logic.
        """
        print(f"\n--- Starting Microarray Pipeline for {self.gse_id} ---")
        
        # 1. Download Raw Data
        print("Downloading supplementary files...")
        self.gse.download_supplementary_files(directory=self.raw_dir)
        
        # --- NEW STEP: Auto-Extract TAR Files ---
        # GEO usually provides a single "GSEXXXX_RAW.tar" file
        tar_files = glob.glob(os.path.join(self.raw_dir, "*.tar"))
        
        if tar_files:
            print(f"Found archive: {tar_files[0]}. Extracting...")
            with tarfile.open(tar_files[0]) as tar:
                tar.extractall(path=self.raw_dir)
                print("Extraction complete.")
        
        # Check for CEL files (they might be .CEL or .CEL.gz)
        # Note: 'affy' package in R can read .CEL.gz directly, so we don't need to gunzip them.
        cel_files = []
        for root, dirs, files in os.walk(self.raw_dir):
            for file in files:
                if file.lower().endswith(".cel") or file.lower().endswith(".cel.gz"):
                    cel_files.append(os.path.join(root, file))
        
        if not cel_files:
            print("ERROR: No .CEL files found after extraction.")
            print(f"Check the folder: {self.raw_dir}")
            return

        print(f"Found {len(cel_files)} CEL files. Running R pipeline...")
        
        # R Script: Now robust to .CEL.gz files
        r_script = f"""
        library(affy)
        library(genefilter)
        library(hgu133plus2.db)

        # 2 & 3. Read CEL files (ReadAffy handles .gz automatically)
        # We specify celfile.path to ensure it looks in the correct dir
        raw_data <- ReadAffy(celfile.path="{self.raw_dir}")
        
        # RMA Normalization
        eset <- rma(raw_data)
        
        # 4. Filter lower 25% IQR
        # We use varFilter with var.func=IQR and var.cutoff=0.25
        filtered_eset <- varFilter(eset, var.func=IQR, var.cutoff=0.25, filterByQuantile=TRUE)
        
        # Extract expression matrix
        filtered_exprs <- exprs(filtered_eset)
        
        # Return as dataframe
        as.data.frame(filtered_exprs)
        """
        
        # Execute R code
        try:
            norm_data = r(r_script)
            output_file = os.path.join(self.processed_dir, f"{self.gse_id}_microarray_processed.csv")
            norm_data.to_csv(output_file)
            print(f"Microarray processing complete. Saved to: {output_file}")
        except Exception as e:
            print("R execution failed. Error details:")
            print(e)

    def run_rnaseq_pipeline(self, path_to_reference_genome):
        """
        Replicates RNA-seq Pipeline (Paper Section 2.3):
        1. Download FASTQ (via SRA).
        2. QC (FastQC) & Trim (Trimmomatic).
        3. Align to UCSC reference.
        4. Calculate TPM.
        5. Filter (Remove genes with sum=0).
        6. VST Transformation (DESeq2).
        """
        print(f"\n--- Starting RNA-seq Pipeline for {self.gse_id} ---")
        
        # 1. Download SRA/FASTQ
        # Requires SRA Toolkit installed
        print("Downloading FASTQ files via SRA...")
        self.gse.download_SRA(directory=self.raw_dir, filetype='fastq', keep_sra=False)
        
        fastq_files = glob.glob(os.path.join(self.raw_dir, "*.fastq"))
        if not fastq_files:
            print("No FASTQ files found.")
            return

        # 2-3. Processing Loop (QC -> Trim -> Align -> Count)
        # Note: This assumes external tools are in your PATH.
        print("Starting QC, Trimming, and Alignment (This may take a long time)...")
        
        for fq in fastq_files:
            base_name = os.path.basename(fq).replace(".fastq", "")
            
            # A. FASTQC
            subprocess.run(["fastqc", fq, "-o", self.processed_dir])
            
            # B. Trimmomatic (Paper: "Low-quality reads and residual adaptor sequences were trimmed")
            # Example command structure (user must provide correct adapters path)
            trimmed_fq = os.path.join(self.processed_dir, f"{base_name}_trimmed.fastq")
            cmd_trim = [
                "java", "-jar", "trimmomatic.jar", "SE", 
                fq, trimmed_fq, 
                "ILLUMINACLIP:adapters.fa:2:30:10", "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:15", "MINLEN:36"
            ]
            # subprocess.run(cmd_trim) # Uncomment to run if configured
            
            # C. Alignment (Paper: "aligned to the USCS reference transcriptome")
            # Usually implies STAR, HISAT2, or TopHat.
            # cmd_align = ["hisat2", "-x", path_to_reference_genome, "-U", trimmed_fq, "-S", f"{base_name}.sam"]
            # subprocess.run(cmd_align)
            
            # D. Quantification (Counts) -> This results in a count matrix
            pass 

        # 4. Load Count Matrix (Simulation of step)
        # Assuming you generated 'counts.csv' from the steps above
        # counts = pd.read_csv("counts.csv")
        
        print("Note: Actual alignment/counting skipped in script (requires huge reference files).")
        print("Mocking count data for demonstration of Normalization/Filtering steps...")
        # Create dummy count data for demonstration
        counts = pd.DataFrame(np.random.randint(0, 100, size=(1000, 10)), columns=[f"Sample_{i}" for i in range(10)])
        
        # 4. Calculate TPM (Transcripts Per Million)
        # (Simplified TPM calculation: reads / gene_length_kb / scaling_factor)
        # Paper says: "read counts... obtained, and transcripts per million (TPM) values were calculated."
        print("Calculating TPM...")
        # In real scenario, divide by gene length first. Here we assume lengths are equal for demo.
        rpk = counts # / gene_length_kb
        tpm = rpk.div(rpk.sum(axis=0), axis=1) * 1e6
        
        # 5. Filter (Paper: "removing all genes with a sum of 0 across all samples")
        print("Filtering Genes (Sum > 0)...")
        counts_filtered = counts[counts.sum(axis=1) > 0]
        
        # 6. VST Transformation via DESeq2
        # Paper: "subjected to variance-stabilizing transformation (VST) using the DESeq2"
        print("Running VST transformation via DESeq2 (R)...")
        
        with (rpy2.robjects.default_converter + pandas2ri.converter).context():
            r_counts = pandas2ri.py2rpy(counts_filtered)
            
        r_vst_script = """
        library(DESeq2)
        run_vst <- function(count_df) {
            # Create minimal colData
            colData <- data.frame(condition=factor(rep("A", ncol(count_df))))
            dds <- DESeqDataSetFromMatrix(countData = round(count_df), colData = colData, design = ~ 1)
            
            # Run VST
            vsd <- vst(dds, blind=TRUE)
            return(assay(vsd))
        }
        """
        r(r_vst_script)
        run_vst = r['run_vst']
        vst_matrix = run_vst(r_counts)
        
        # Save
        vst_df = pd.DataFrame(vst_matrix, index=counts_filtered.index, columns=counts_filtered.columns)
        output_vst = os.path.join(self.processed_dir, f"{self.gse_id}_rnaseq_vst.csv")
        vst_df.to_csv(output_vst)
        print(f"RNA-seq processing complete. Saved VST data to: {output_vst}")


# --- Usage Example ---
if __name__ == "__main__":
    install_r_dependencies() 
    # Example GSE ID from the prompt/paper context or user list
    series_id = "GSE12345" 
    
    pipeline = StudyPipeline(series_id)
    
    
    # Example logic TODO fix this later:
    if "affymetrix" in str(pipeline.gse.metadata).lower():
        pipeline.run_microarray_pipeline()
    else:
        # Note: You need a path to a genome reference (e.g., hg19/hg38) for RNA-seq alignment
        pipeline.run_rnaseq_pipeline(path_to_reference_genome="/path/to/hg38_index")