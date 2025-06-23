#!/bin/bash
#SBATCH --job-name=trimming_qc
#SBATCH --output=logs/trimming_qc_%j.out
#SBATCH --error=logs/trimming_qc_%j.err
#SBATCH --time=10:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=5
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=saleka@ethz.ch

# Activating conda environment
source ~/.bashrc
conda activate ethepigen

# Define paths
TRIM_JAR="/cluster/scratch/saleka/Trimmomatic-0.39/trimmomatic-0.39.jar"
ADAPTERS="/cluster/scratch/saleka/Trimmomatic-0.39/adapters/TruSeq3-PE.fa"
RAW_DIR="rawdata"
OUTDIR="trimmed_new"
QC_BEFORE="qc_raw_new"
QC_AFTER="qc_trimmed_new"

mkdir -p $OUTDIR $QC_BEFORE $QC_AFTER

# Check MultiQC availability
if ! command -v multiqc &> /dev/null; then
    echo "Error: multiqc not found in environment. Exiting."
    exit 1
fi

# QC Before trimming
if [ -z "$(ls -A $QC_BEFORE/*.html 2>/dev/null)" ]; then
    echo "▶ Running FastQC on raw reads..."
    for fq in ${RAW_DIR}/*.fastq.gz; do
        fastqc -t 4 -o $QC_BEFORE "$fq"
    done
    multiqc -o $QC_BEFORE $QC_BEFORE
else
    echo "Raw FastQC already done — skipping."
fi


# Trimming 
echo "Starting Trimmomatic trimming..."
while read sample; do
    echo "Processing $sample"
    java -jar $TRIM_JAR PE -threads 4 \
        ${RAW_DIR}/${sample}_1.fastq.gz ${RAW_DIR}/${sample}_2.fastq.gz \
        ${OUTDIR}/${sample}_1P.fastq.gz ${OUTDIR}/${sample}_1U.fastq.gz \
        ${OUTDIR}/${sample}_2P.fastq.gz ${OUTDIR}/${sample}_2U.fastq.gz \
        ILLUMINACLIP:${ADAPTERS}:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36
done < /cluster/scratch/saleka/rawdata/sample_ids.txt

# QC after trimming 
if [ -z "$(ls -A $QC_AFTER/*.html 2>/dev/null)" ]; then
    echo "▶ Running FastQC on trimmed reads..."
    for fq in ${OUTDIR}/*P.fastq.gz; do
        fastqc -t 4 -o $QC_AFTER "$fq"
    done
    multiqc -o $QC_AFTER $QC_AFTER
else
    echo "Trimmed FastQC already done — skipping."
fi

echo "QC and trimming completed successfully."
