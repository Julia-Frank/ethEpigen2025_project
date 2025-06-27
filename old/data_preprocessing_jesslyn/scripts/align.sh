#!/bin/bash
#SBATCH --job-name=align
#SBATCH --output=logs/align_%A_%a.out
#SBATCH --error=logs/align_%A_%a.err
#SBATCH --time=3:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-12%4

set -x  # print commands
module load stack/2024-06 gcc/12.2.0
module load picard/2.26.2

ref=/cluster/scratch/jjesslyn/genome/mm10/mm10
ncpus=4

mkdir -p aligned logs

samples=( $(ls rawdata/fastq/*_1.fastq.gz | xargs -n1 basename | sed 's/_1.fastq.gz//') )
base=${samples[$SLURM_ARRAY_TASK_ID-1]}

echo "Processing sample: $base"
trimdir=trimmed
aligndir=aligned

# ALIGNMENT

if [ -f "$aligndir/$base.bam" ]; then
    echo $aligndir/$base".bam found; skipping"
else
  # Alignment
  echo "[$base] Starting alignment..."
  bowtie2 -p $ncpus --dovetail --no-mixed --no-discordant -I 15 -X 2000 -x $ref \
    -1 $trimdir/"$base"_1.paired.fastq.gz -2 $trimdir/"$base"_2.paired.fastq.gz 2> $aligndir/"$base".bowtie2 | \
  samtools view -bS - | \
  samtools sort -@ $ncpus -m 4G -o $aligndir/"$base".bam
  
  # Check if alignment and sorting succeeded
  if [ ! -f "$aligndir/"$base".bam" ]; then
      echo "Alignment failed, exiting."
      exit 1
  fi
  
  echo "[$base] Filtering reads ..."
  samtools view -b -q 20 $aligndir/"$base".bam > $aligndir/"$base".mapq20.bam

  # Replace original BAM only after successful filtering
  if [ -f "$aligndir/"$base".mapq20.bam" ]; then
      mv $aligndir/"$base".mapq20.bam $aligndir/"$base".bam
  else
      echo "Filtering failed, original BAM preserved."
      exit 1
  fi

  echo "[$base] Removing duplicates ..."
  picard AddOrReplaceReadGroups \
        I=$aligndir/"$base".bam \
        O=$aligndir/"$base".rg.bam \
        RGID=1 \
        RGLB=lib1 \
        RGPL=ILLUMINA \
        RGPU=unit1 \
        RGSM=$base

  if [ -f "$aligndir/"$base".rg.bam" ]; then
      mv $aligndir/"$base".rg.bam $aligndir/"$base".bam
  else
      echo "AddOrReplaceReadGroups failed, exiting."
      exit 1
  fi

  picard MarkDuplicates I=$aligndir/"$base".bam O=$aligndir/"$base".bam.2 \
    M=$aligndir/"$base".picard.dupMetrics.txt

  # Only replace original file if Picard was successful
  if [ -f "$aligndir/"$base".bam.2" ]; then
      mv $aligndir/"$base".bam.2 $aligndir/"$base".bam
  else
      echo "Duplicate removal failed. Original BAM is preserved."
      exit 1
  fi

 echo "[$base] Indexing BAM..."
 samtools index $aligndir/"$base".bam

fi
