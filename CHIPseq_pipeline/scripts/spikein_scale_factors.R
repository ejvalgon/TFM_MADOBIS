#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript spikein_scale_factors.R <spike_reads_table.tsv> <scale_factors.tsv>")
}

input_tsv <- args[1]
output_tsv <- args[2]

df <- read.delim(
  input_tsv,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

required_cols <- c("sample", "spike_reads")
missing_cols <- setdiff(required_cols, colnames(df))

if (length(missing_cols) > 0) {
  stop("Missing columns in input file: ", paste(missing_cols, collapse = ", "))
}

if (any(duplicated(df$sample))) {
  duplicated_samples <- unique(df$sample[duplicated(df$sample)])
  stop("Duplicated samples found in input file: ", paste(duplicated_samples, collapse = ", "))
}

spike_reads <- as.numeric(df$spike_reads)

if (any(is.na(spike_reads))) {
  stop("NA or non-numeric values found in spike_reads")
}

if (any(spike_reads <= 0)) {
  bad <- df$sample[spike_reads <= 0]
  stop(
    "Samples with zero or negative spike-in reads found. Size factors cannot be calculated: ",
    paste(bad, collapse = ", ")
  )
}

names(spike_reads) <- df$sample

# 1) Original scale factor based on the median spike-in read count
median_count <- median(spike_reads)
scale_factor_median <- median_count / spike_reads

# 2) DESeq2 size factor using only spike-in reads
count_mat <- matrix(
  as.integer(round(spike_reads)),
  nrow = 1
)

rownames(count_mat) <- "spike_in"
colnames(count_mat) <- names(spike_reads)

deseq2_size_factor <- estimateSizeFactorsForMatrix(count_mat)

# bamCoverage requires the inverse of the DESeq2 size factor,
# because bamCoverage multiplies the signal by --scaleFactor.
scale_factor_deseq2 <- 1 / deseq2_size_factor

out <- data.frame(
  sample = names(spike_reads),
  spike_reads = as.integer(round(spike_reads)),
  scale_factor_median = as.numeric(scale_factor_median),
  deseq2_size_factor = as.numeric(deseq2_size_factor),
  scale_factor_deseq2 = as.numeric(scale_factor_deseq2)
)

write.table(
  out,
  file = output_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)