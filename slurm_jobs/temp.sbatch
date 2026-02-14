#!/bin/sh

#SBATCH --partition=general

#SBATCH --qos=short

#SBATCH --mail-type=END

#SBATCH --time=3:00:00

#SBATCH --ntasks=1

#SBATCH --cpus-per-task=1

#SBATCH --mem=12G


# source ~/.bashrc

conda install -c bioconda -c conda-forge bioconductor-oligo bioconductor-biobase bioconductor-geoquery r-rsqlite r-essentials