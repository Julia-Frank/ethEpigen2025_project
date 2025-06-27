#!/bin/bash
#SBATCH --job-name=atac_pipeline
#SBATCH --output=logs/atac_pipeline_%j.out
#SBATCH --error=logs/atac_pipeline_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=6G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=saleka@ethz.ch  

source ~/.bashrc
conda activate atac_env

# Run R Markdown document
echo "Running alignment and index building..."
Rscript alignment_indexing.R

echo "Pipeline complete."

