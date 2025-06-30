#!/bin/bash
#SBATCH --job-name=peak
#SBATCH --output=logs/peak_%A_%a.out
#SBATCH --error=logs/peak_%A_%a.err
#SBATCH --time=8:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-16%4

mkdir -p peaks

samples=( $(find aligned -maxdepth 1 -type f -name "*.bam" ! -name "*.bam.*" ! -name "*.rg.bam" | xargs -n1 basename | sed 's/.bam//') )
base=${samples[$SLURM_ARRAY_TASK_ID-1]}

echo "Starting task $SLURM_ARRAY_TASK_ID for sample $base"
bam=aligned/"$base".bam

if [ ! -f "$bam" ]; then
  echo "Missing BAM file: $bam"
  exit 1
fi

if [ -f peaks/"$base"_peaks.narrowPeak ]; then
  echo "$base already processed, skipping"
  exit 0
fi

macs3 callpeak --outdir peaks -n $base --gsize mm -t $bam --format BAMPE

echo "Peak calling complete for $base"
