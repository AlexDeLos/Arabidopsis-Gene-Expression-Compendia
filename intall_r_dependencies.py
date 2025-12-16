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
    if not robjects.r.installed_packages().rx(True, "Package").count("BiocManager"):
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

# --- PLACE THIS AT THE VERY BOTTOM, INSIDE "if __name__ ==..." ---
if __name__ == "__main__":
    # 1. Install R packages first
    install_r_dependencies() 
    
    # 2. Then run your pipeline
    series_id = "GSE12345"
    pipeline = StudyPipeline(series_id)
    # ... rest of your code