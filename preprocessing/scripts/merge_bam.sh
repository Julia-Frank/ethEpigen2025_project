#!/bin/bash
#SBATCH --job-name=merge_bam
#SBATCH --output=logs/merge_bam.out
#SBATCH --error=logs/merge_bam.err
#SBATCH --time=05:00:00              # max runtime hh:mm:ss
#SBATCH --cpus-per-task=4             # number of CPU cores
#SBATCH --mem-per-cpu=16

bamdir=aligned

echo "Merging uninjured BAMs..."
samtools merge -@ 4 $bamdir/uninjured_merged.bam \
  $bamdir/SRR14371953.bam \
  $bamdir/SRR14371954.bam

echo "Merging injured BAMs..."
samtools merge -@ 4 $bamdir/injured_merged.bam \
  $bamdir/SRR14371955.bam \
  $bamdir/SRR14371956.bam \
  $bamdir/SRR14371957.bam

echo "Merging mdx BAMs..."
samtools merge -@ 4 $bamdir/mdx_merged.bam \
  $bamdir/SRR15343189.bam \
  $bamdir/SRR15343190.bam \
  $bamdir/SRR15343191.bam

echo "Merging wt BAMs..."
samtools merge -@ 4 $bamdir/wt_merged.bam \
  $bamdir/SRR15343192.bam \
  $bamdir/SRR15343193.bam \
  $bamdir/SRR15343194.bam \
  $bamdir/SRR15343195.bam

echo "Indexing merged BAM files..."
for bam in $bamdir/uninjured_merged.bam $bamdir/injured_merged.bam $bamdir/mdx_merged.bam $bamdir/wt_merged.bam; do
  samtools index $bam
done

echo "Done!"
