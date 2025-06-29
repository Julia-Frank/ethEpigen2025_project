#!/bin/bash
#SBATCH --job-name=align
#SBATCH --output=logs/align_%A_%a.out
#SBATCH --error=logs/align_%A_%a.err
#SBATCH --time=30:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-12%4

set -x  # print commands

adapters=/cluster/home/jjesslyn/miniconda3/envs/epigenomics/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa
ncpus=4

mkdir -p trimmed


samples=( $(ls rawdata/fastq/*_1.fastq.gz | xargs -n1 basename | sed 's/_1.fastq.gz//') )
base=${samples[$SLURM_ARRAY_TASK_ID-1]}

echo "Processing sample: $base"

# TRIMMING
 
trimdir=trimmed

if [ -f "$trimdir/"$base"_1.paired.fastq.gz" ]; then
  echo $trimdir/$base"_*fastq.gz found; skipping"
else
  echo "Running Trimmomatic on $base..."
  trimmomatic PE -threads $ncpus -summary $trimdir/"$base".stats -phred33 \
   rawdata/fastq/"$base"_1.fastq.gz rawdata/fastq/"$base"_2.fastq.gz \
   $trimdir/"$base"_1.paired.fastq.gz $trimdir/"$base"_1.unpaired.fastq.gz \
   $trimdir/"$base"_2.paired.fastq.gz $trimdir/"$base"_2.unpaired.fastq.gz \
   ILLUMINACLIP:$adapters:2:15:4:4:true LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:25
  echo "Finished trimming $base"
fi
