#!/bin/bash
#SBATCH --account=testusers
#SBATCH --partition=all
#SBATCH --time=240:00:00
#SBATCH --mem=46M
#SBATCH --cpus-per-task=1
#SBATCH --output=/home/alexdelossanto/linuxhome/Dataset_fusion_Microarray/monitor.log

./storage_monitor.sh