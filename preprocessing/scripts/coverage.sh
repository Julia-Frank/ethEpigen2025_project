#!/bin/bash
#SBATCH --job-name=coverage
#SBATCH --output=log_coverage/%A_%a.out
#SBATCH --error=log_coverage/%A_%a.err
#SBATCH --time=8:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-16%4

samples=( $(find aligned -maxdepth 1 -name "*.bam" -not -name "*rg*" -not -name "*.bam.2" -not -name "*.bai" | xargs -n1 basename | sed 's/.bam//') )
base=${samples[$SLURM_ARRAY_TASK_ID-1]}

bamdir=aligned
bwdir=tracks
mkdir -p $bwdir log_coverage

echo "Starting task $SLURM_ARRAY_TASK_ID for sample $base"
if [ -f "tracks/$base.bw" ]; then
    echo "tracks/$base.bw found; skipping"
else
  bamCoverage -p 4 --binSize 20 --effectiveGenomeSize 2652783500 --normalizeUsing CPM -b $bamdir/"$base".bam -o $bwdir/"$base".bw
fi
