#!/bin/bash

mkdir -p peaks_new

for bam in bam_new/*.bam; do
  sample=$(basename "$bam" .bam)
  echo "Calling peaks for $sample..."
  
  macs3 callpeak \
    -t "$bam" \
    -f BAMPE \
    -n "$sample" \
    -g mm \
    --outdir "peaks_new/$sample" 
done

echo "Peak calling finished."
