#!/bin/bash
#SBATCH --job-name=install_packages
#SBATCH --partition=general
#SBATCH --time=2:00:00       # 45 mins is plenty for installation
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2    # Give the solver a little breathing room
#SBATCH --mem=25G            # Good choice. The solver can be memory hungry.
#SBATCH --output=install_log.txt

echo "------------------------------------------------"
echo "Active Conda Environment: $CONDA_DEFAULT_ENV"
echo "Python Location: $(which python)"
echo "------------------------------------------------"

conda env create -f environment.yml -p ./my_cluster_env