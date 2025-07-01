#!/bin/bash
#SBATCH --job-name=common_peaks.sh
#SBATCH --output=logs/common_peaks.out
#SBATCH --error=logs/common_peaks.err
#SBATCH --time=8:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=4

peakdir=peaks_corrected

echo "Finding replicate peaks for UI"
intersectBed -wo -a $peakdir/uninjured_merged_peaks.narrowPeak -b $peakdir/SRR14371953_peaks.narrowPeak | 
awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | cut -f 1-10 | sort | uniq | 
intersectBed -wo -a stdin -b $peakdir/SRR14371954_peaks.narrowPeak | 
awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | cut -f 1-10 | sort | uniq > $peakdir/ui_final_peaks.narrowPeak

echo "Finding replicate peaks for I"
intersectBed -wo -a $peakdir/injured_merged_peaks.narrowPeak -b $peakdir/SRR14371955_peaks.narrowPeak | 
awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | cut -f 1-10 | sort | uniq | 
intersectBed -wo -a stdin -b $peakdir/SRR14371956_peaks.narrowPeak | 
awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | cut -f 1-10 | sort | uniq | 
intersectBed -wo -a stdin -b $peakdir/SRR14371957_peaks.narrowPeak | 
awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | cut -f 1-10 | sort | uniq > $peakdir/i_final_peaks.narrowPeak

echo "Finding replicate peaks for MDX"
intersectBed -wo -a $peakdir/mdx_merged_peaks.narrowPeak -b $peakdir/SRR15343189_peaks.narrowPeak |
awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | cut -f 1-10 | sort | uniq | 
intersectBed -wo -a stdin -b $peakdir/SRR15343190_peaks.narrowPeak | 
awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | cut -f 1-10 | sort | uniq | 
intersectBed -wo -a stdin -b $peakdir/SRR15343191_peaks.narrowPeak | 
awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | cut -f 1-10 | sort | uniq > $peakdir/mdx_final_peaks.narrowPeak

echo "Finding replicate peaks for WT"
intersectBed -wo -a $peakdir/wt_merged_peaks.narrowPeak -b $peakdir/SRR15343192_peaks.narrowPeak |
awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | cut -f 1-10 | sort | uniq | 
intersectBed -wo -a stdin -b $peakdir/SRR15343193_peaks.narrowPeak | 
awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | cut -f 1-10 | sort | uniq | 
intersectBed -wo -a stdin -b $peakdir/SRR15343194_peaks.narrowPeak | 
awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | cut -f 1-10 | sort | uniq | 
intersectBed -wo -a stdin -b $peakdir/SRR15343195_peaks.narrowPeak | 
awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | cut -f 1-10 | sort | uniq > $peakdir/wt_final_peaks.narrowPeak


