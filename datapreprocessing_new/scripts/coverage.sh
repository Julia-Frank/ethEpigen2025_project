#!/bin/bash

BAMDIR="bam_new"
SORTEDDIR="bam_new/sorted"
BWDIR="tracks"

mkdir -p "$SORTEDDIR"
mkdir -p "$BWDIR"

for bam in "$BAMDIR"/*.bam; do
  base=$(basename "$bam" .bam)
  sorted_bam="$SORTEDDIR/${base}.sorted.bam"

  if [[ -f "$sorted_bam" ]]; then
    echo "Skipping sort: $sorted_bam already exists."
  else
    echo "Sorting $bam -> $sorted_bam"
    samtools sort -@ 8 -o "$sorted_bam" "$bam"
    echo "Indexing $sorted_bam"
    samtools index "$sorted_bam"
  fi
done  

for sorted_bam in "$SORTEDDIR"/*.sorted.bam; do
  sample=$(basename "$sorted_bam" .sorted.bam)
  echo "Generating bigWig for $sample"
  
  bamCoverage \
    -b "$sorted_bam" \
    -o "$BWDIR/${sample}.bw" \
    --binSize 20 \
    --extendReads 100 \
    --normalizeUsing CPM \
    --numberOfProcessors 8 
done

