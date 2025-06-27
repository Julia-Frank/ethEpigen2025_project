# ATAC-seq Paired-End Alignment & Indexing (mm10)


suppressPackageStartupMessages({
  library(AnnotationHub)
  library(Rsubread)
  library(rtracklayer)
  library(Biostrings)
})

raw_dir <- "/cluster/scratch/saleka/rawdata"
samples <- readLines(file.path(raw_dir, "sample_ids.txt"))

fastq_pairs <- lapply(samples, function(s) {
  fq1 <- file.path(raw_dir, paste0(s, "_1.fastq.gz"))
  fq2 <- file.path(raw_dir, paste0(s, "_2.fastq.gz"))
  c(R1 = fq1, R2 = fq2)
})
names(fastq_pairs) <- samples

genome_fasta <- "mm10_genome/mm10.fa.gz"

if (!file.exists("mm10_genome/rsubread.00.b.array")) { 
  Rsubread::buildindex("mm10_genome/rsubread", reference = genome_fasta)
}

dir.create("bam_new", showWarnings = FALSE)

bam_files <- mapply(function(sample, fq) {
  out_bam <- file.path("bam_new", paste0(sample, ".bam"))

  Rsubread::align(
    index = "mm10_genome/rsubread",
    readfile1 = fq["R1"],
    readfile2 = fq["R2"],
    output_file = out_bam,
    nthreads = 6,
    input_format = "gzFASTQ",
    phredOffset = 33,
    PE_orientation = "fr",
    type = "dna",
    unique = TRUE
  )

  return(out_bam)
}, sample = names(fastq_pairs), fq = fastq_pairs, SIMPLIFY = FALSE)

message("Alignment complete.")

