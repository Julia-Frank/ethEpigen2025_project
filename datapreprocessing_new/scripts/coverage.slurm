#!/bin/bash
#SBATCH --job-name=bigwig_gen
#SBATCH --output=logs/bigwig_gen_%j.out
#SBATCH --error=logs/bigwig_gen_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=4

# Load required modules or activate conda env
module load deeptools

# If using conda
source ~/.bashrc
conda activate ethepigen

# Run R script
./coverage.sh
