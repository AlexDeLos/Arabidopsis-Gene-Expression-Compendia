#!/bin/sh

#SBATCH --partition=general

#SBATCH --qos=medium

#SBATCH --mail-type=END

#SBATCH --time=8:00:00

#SBATCH --ntasks=1

#SBATCH --cpus-per-task=8

#SBATCH --mem=90G

#SBATCH --gres=gpu

# source ~/.bashrc

module use /opt/insy/modulefiles
module load miniconda/3.9
conda activate download_env
python3 src/data_importing/new_download.py -o '/tudelft.net/staff-umbrella/GeneExpressionStorage/'
echo "done with main.py"