#!/bin/bash
#SBATCH --job-name=install_packages
#SBATCH --partition=general
#SBATCH --time=0:45:00       # 45 mins is plenty for installation
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2    # Give the solver a little breathing room
#SBATCH --mem=16G            # Good choice. The solver can be memory hungry.
#SBATCH --output=install_log.txt

echo "------------------------------------------------"
echo "Active Conda Environment: $CONDA_DEFAULT_ENV"
echo "Python Location: $(which python)"
echo "------------------------------------------------"

# 4. Install
conda install -c conda-forge mamba -y

echo "mama installed ---------------------------------------------------------"
mamba install -y -c bioconda -c conda-forge bioconductor-oligo bioconductor-biobase bioconductor-geoquery r-rsqlite r-essentials