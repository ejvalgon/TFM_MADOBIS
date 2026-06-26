# ==============================================================================
# Count RNAPII ChIP-seq fragments in internal gene-body regions
# ==============================================================================
# This script quantifies paired-end RNAPII ChIP-seq fragments within internal
# gene-body regions defined by a strand-aware trimmed GTF annotation.
#
# Input annotation:
#   - protein-coding genes at least 3 kb long
#   - 2500 bp removed from the TSS side
#   - 500 bp removed from the TES side
#
# Main downstream output:
#   - RNAPII_ChIPseq_internal_gene_body_counts.tsv
#     Used as input for differential RNAPII occupancy analysis with DESeq2.
# ============================================================================== 

library(Rsubread)

# ------------------------------------------------------------------------------
# 1. Input and output paths
# ------------------------------------------------------------------------------

bam_dir <- "Data_final_CPT/split_bam/main"

gtf_file <- file.path(
  "../Script_Limpios",
  "gtf_processing",
  "gtf_trimmed",
  "genes_protein_coding_min3kb_trimmed_TSS2500_TES500.gtf"
)

output_dir <- file.path(
  "../Script_Limpios",
  "Count_Matrix_ChIPseq",
  "RNAPII_internal_gene_body"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 2. Locate and validate input files
# ------------------------------------------------------------------------------

bam_files <- list.files(
  path = bam_dir,
  pattern = "\\.bam$",
  full.names = TRUE
)

if (length(bam_files) == 0L) {
  stop("No BAM files were found in: ", bam_dir, call. = FALSE)
}

if (!file.exists(gtf_file)) {
  stop("The trimmed GTF file was not found: ", gtf_file, call. = FALSE)
}

# Assign standard sample names from the original BAM filenames.
sample_names <- vapply(
  basename(bam_files),
  FUN.VALUE = character(1),
  FUN = function(file_name) {
    if (grepl("^R1_.*_CPT_", file_name)) return("CPT_R1")
    if (grepl("^R2_.*_CPT_", file_name)) return("CPT_R2")
    if (grepl("^R1_.*_DM_", file_name))  return("DMSO_R1")
    if (grepl("^R2_.*_DM_", file_name))  return("DMSO_R2")

    NA_character_
  }
)

if (anyNA(sample_names)) {
  stop(
    "Sample names could not be assigned to the following BAM files: ",
    paste(basename(bam_files[is.na(sample_names)]), collapse = ", "),
    call. = FALSE
  )
}

if (anyDuplicated(sample_names)) {
  stop("More than one BAM file was assigned to the same sample.", call. = FALSE)
}

expected_samples <- c("DMSO_R1", "CPT_R1", "DMSO_R2", "CPT_R2")

if (!setequal(sample_names, expected_samples)) {
  stop(
    "The detected samples do not match the expected set: ",
    paste(expected_samples, collapse = ", "),
    call. = FALSE
  )
}

# Reorder BAM files before counting so all downstream outputs use the same order.
bam_files <- bam_files[match(expected_samples, sample_names)]
sample_names <- expected_samples

# ------------------------------------------------------------------------------
# 3. Count paired-end fragments with featureCounts
# ------------------------------------------------------------------------------

fc <- featureCounts(
  files = bam_files,
  annot.ext = gtf_file,
  isGTFAnnotationFile = TRUE,
  GTF.featureType = "gene",
  GTF.attrType = "gene_id",
  useMetaFeatures = TRUE,
  isPairedEnd = TRUE,
  requireBothEndsMapped = TRUE,
  countReadPairs = TRUE,
  allowMultiOverlap = FALSE,
  nthreads = 8
)

# ------------------------------------------------------------------------------
# 4. Prepare count and annotation tables
# ------------------------------------------------------------------------------

count_matrix <- fc$counts
colnames(count_matrix) <- sample_names

count_table <- cbind(fc$annotation, count_matrix)

# ------------------------------------------------------------------------------
# 5. Export final outputs
# ------------------------------------------------------------------------------

# Count matrix used by the downstream RNAPII DESeq2 analysis.
write.table(
  count_matrix,
  file = file.path(output_dir, "RNAPII_ChIPseq_internal_gene_body_counts.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

# Counts together with the genomic coordinates and region lengths reported by
# featureCounts. This table is retained for traceability and length-based analyses.
write.table(
  count_table,
  file = file.path(
    output_dir,
    "RNAPII_ChIPseq_internal_gene_body_counts_with_annotation.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# featureCounts read-assignment summary for quality control.
write.table(
  fc$stat,
  file = file.path(
    output_dir,
    "RNAPII_ChIPseq_internal_gene_body_assignment_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
