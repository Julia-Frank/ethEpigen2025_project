#!/bin/bash
#SBATCH --job-name=generate_count_matrix
#SBATCH --output=logs/matrix_gen.out
#SBATCH --error=logs/matrix_gen.err
#SBATCH --time=05:00:00              # max runtime hh:mm:ss
#SBATCH --cpus-per-task=4             # number of CPU cores
#SBATCH --mem-per-cpu=16
#SBATCH --mail-type=END,FAIL #Email to me when finished

# Initialize Conda for bash shell 
source /cluster/home/jfrank/miniconda3/etc/profile.d/conda.sh
# Activate the conda environment
conda activate r_env 
# Make sure you have bedtools & subread installed

# Navigate to the directory containing the peaks
cd "/cluster/scratch/jfrank/peaks" || { echo "Directory not found"; exit 1; }

# Check that filepath to the bam files is correct
bamdir="/cluster/scratch/jfrank/aligned"

# Define output folder
outputfolder="/cluster/scratch/jfrank/count_matrices"
mkdir -p "$outputfolder"

# Concatenate all peak files
cat SRR14371953_peaks.narrowPeak SRR14371954_peaks.narrowPeak SRR14371955_peaks.narrowPeak SRR14371956_peaks.narrowPeak SRR14371957_peaks.narrowPeak > all_peaks_UI_I.bed
cat SRR15343189_peaks.narrowPeak SRR15343190_peaks.narrowPeak SRR15343191_peaks.narrowPeak SRR15343192_peaks.narrowPeak SRR15343193_peaks.narrowPeak SRR15343194_peaks.narrowPeak SRR15343195_peaks.narrowPeak > all_peaks_MDX_WT.bed
# Sort the peaks
sort -k1,1 -k2,2n all_peaks_UI_I.bed > all_peaks_UI_I.sorted.bed
sort -k1,1 -k2,2n all_peaks_MDX_WT.bed > all_peaks_MDX_WT.sorted.bed
# Merge overlapping regions to get consensus peaks
bedtools merge -i all_peaks_UI_I.sorted.bed > consensus_peaks_UI_I.bed
bedtools merge -i all_peaks_MDX_WT.sorted.bed > consensus_peaks_MDX_WT.bed
# Convert to SAF
awk 'BEGIN {OFS="\t"} {print "peak"NR, $1, $2+1, $3, "."}' consensus_peaks_UI_I.bed > consensus_peaks_UI_I.saf
awk 'BEGIN {OFS="\t"} {print "peak"NR, $1, $2+1, $3, "."}' consensus_peaks_MDX_WT.bed > consensus_peaks_MDX_WT.saf

# Define SRR IDs fo both experiments explicitly
ui_i_srrs=(
  SRR14371953
  SRR14371954
  SRR14371955
  SRR14371956
  SRR14371957
)
mdx_wt_srrs=(
  SRR15343189
  SRR15343190
  SRR15343191
  SRR15343192
  SRR15343193
  SRR15343194
  SRR15343195
)

# Select BAM files matching SRRs
bamfiles_ui_i=()
for srr in "${ui_i_srrs[@]}"; do
  bamfile="${bamdir}/${srr}.bam"
  if [[ -f "$bamfile" ]]; then
    bamfiles_ui_i+=("$bamfile")
  else
    echo "Warning: BAM file $bamfile not found"
  fi
done

bamfiles_mdx_wt=()
for srr in "${mdx_wt_srrs[@]}"; do
  bamfile="${bamdir}/${srr}.bam"
  if [[ -f "$bamfile" ]]; then
    bamfiles_mdx_wt+=("$bamfile")
  else
    echo "Warning: BAM file $bamfile not found"
  fi
done

# Count how many sequencing reads map to the consensus peaks
featureCounts -a consensus_peaks_MDX_WT.saf -F SAF -o "${outputfolder}/counts_MDX_WT.txt" -T 4 "${bamfiles_mdx_wt[@]}"
echo "MDX_WT count matrix generated: ${outputfolder}/counts_MDX_WT.txt"

featureCounts -a consensus_peaks_UI_I.saf -F SAF -o "${outputfolder}/counts_UI_I.txt" -T 4 "${bamfiles_ui_i[@]}"
echo "UI_I count matrix generated: ${outputfolder}/counts_UI_I.txt"

echo "Done!"
