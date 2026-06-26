################################################################################
# Prepare protein-coding gene annotations for RNAPII ChIP-seq analyses
################################################################################
#
# This script generates the gene annotation files used in downstream RNAPII
# ChIP-seq analyses.
#
# Main outputs:
#   1) A filtered protein-coding gene annotation.
#   2) A strand-aware trimmed gene-body annotation.
#
# The main annotation keeps gene-level protein-coding features located on
# standard chromosomes, with defined strand and a configurable minimum gene
# length. The trimmed annotation removes 2.5 kb from the TSS side and 0.5 kb
# from the TES side of each gene, and is intended for internal gene-body signal
# analyses.
#
################################################################################

# ==============================================================================
# 1. Packages
# ==============================================================================

suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(GenomeInfoDb)
})

# ==============================================================================
# 2. Input, output and filtering parameters
# ==============================================================================

# Input GTF annotation.
gtf_file <- "Homo_sapiens.GRCh38.115.gtf"

# Output directories.
outdir <- "../Script_Limpios/gtf_processing"
outdir_trimmed <- file.path(outdir, "gtf_trimmed")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir_trimmed, recursive = TRUE, showWarnings = FALSE)

# Minimum gene length retained in the main protein-coding annotation.
# Set to 0 to disable length filtering.
# Current analysis: 3000 bp, to ensure consistency with downstream internal
# gene-body analyses after TSS/TES trimming
min_gene_length <- 3000

# Strand-aware trimming parameters for internal gene-body analyses.
tss_trim <- 2500
tes_trim <- 500

# Output file names.
length_tag <- if (min_gene_length > 0) {
  paste0("min", min_gene_length / 1000, "kb")
} else {
  "no_min_length"
}

protein_coding_gtf <- file.path(
  outdir,
  paste0("genes_protein_coding_", length_tag, ".gtf")
)

trimmed_gtf <- file.path(
  outdir_trimmed,
  paste0(
    "genes_protein_coding_", length_tag,
    "_trimmed_TSS", tss_trim,
    "_TES", tes_trim,
    ".gtf"
  )
)

summary_file <- file.path(outdir, "annotation_filtering_summary.tsv")

# ==============================================================================
# 3. Helper functions
# ==============================================================================

append_summary <- function(step, n_genes, output_file = NA_character_) {
  data.frame(
    step = step,
    n_genes = n_genes,
    output_file = output_file,
    stringsAsFactors = FALSE
  )
}

trim_gene_body <- function(gr, tss_trim = 2500, tes_trim = 500) {
  gr_trimmed <- gr
  gene_strand <- as.character(strand(gr_trimmed))

  plus_genes <- gene_strand == "+"
  minus_genes <- gene_strand == "-"

  # Genes on the + strand: TSS = start, TES = end.
  start(gr_trimmed)[plus_genes] <- start(gr_trimmed)[plus_genes] + tss_trim
  end(gr_trimmed)[plus_genes] <- end(gr_trimmed)[plus_genes] - tes_trim

  # Genes on the - strand: TSS = end, TES = start.
  start(gr_trimmed)[minus_genes] <- start(gr_trimmed)[minus_genes] + tes_trim
  end(gr_trimmed)[minus_genes] <- end(gr_trimmed)[minus_genes] - tss_trim

  # Remove any regions that may still become invalid after trimming.
  gr_trimmed[width(gr_trimmed) > 0]
}

# ==============================================================================
# 4. Build filtered protein-coding gene annotation
# ==============================================================================

message("Loading GTF annotation: ", gtf_file)
gtf <- rtracklayer::import(gtf_file)

summary_table <- list()

summary_table[[length(summary_table) + 1]] <- append_summary(
  step = "input_gtf_features",
  n_genes = length(gtf)
)

# Keep gene-level features only.
genes <- gtf[gtf$type == "gene"]

summary_table[[length(summary_table) + 1]] <- append_summary(
  step = "gene_level_features",
  n_genes = length(genes)
)

# Keep only standard chromosomes.
genes <- GenomeInfoDb::keepStandardChromosomes(
  genes,
  pruning.mode = "coarse"
)

summary_table[[length(summary_table) + 1]] <- append_summary(
  step = "standard_chromosomes",
  n_genes = length(genes)
)

# Keep only genes with defined strand.
genes <- genes[as.character(strand(genes)) %in% c("+", "-")]

summary_table[[length(summary_table) + 1]] <- append_summary(
  step = "defined_strand",
  n_genes = length(genes)
)

# Keep protein-coding genes.
if (!"gene_biotype" %in% colnames(as.data.frame(mcols(genes)))) {
  stop(
    "Column 'gene_biotype' was not found in the GTF metadata.",
    call. = FALSE
  )
}

genes <- genes[genes$gene_biotype == "protein_coding"]

summary_table[[length(summary_table) + 1]] <- append_summary(
  step = "protein_coding_genes",
  n_genes = length(genes)
)

# Optional minimum gene length filter.
if (min_gene_length > 0) {
  genes <- genes[width(genes) >= min_gene_length]
  
  summary_table[[length(summary_table) + 1]] <- append_summary(
    step = paste0("minimum_gene_length_", min_gene_length, "bp"),
    n_genes = length(genes),
    output_file = protein_coding_gtf
  )
} else {
  summary_table[[length(summary_table) + 1]] <- append_summary(
    step = "no_minimum_gene_length_filter",
    n_genes = length(genes),
    output_file = protein_coding_gtf
  )
}

message("Exporting protein-coding gene annotation: ", protein_coding_gtf)

rtracklayer::export(
  genes,
  protein_coding_gtf,
  format = "gtf"
)

# ==============================================================================
# 5. Build trimmed gene-body annotation
# ==============================================================================

# Genes shorter than or equal to the total trimmed length would become zero-length
# or invalid after trimming, so they are removed before coordinate adjustment.
minimum_length_for_trimming <- tss_trim + tes_trim

genes_for_trimming <- genes[width(genes) > minimum_length_for_trimming]
summary_table[[length(summary_table) + 1]] <- append_summary(
  step = paste0("long_enough_for_trimming_>", minimum_length_for_trimming, "bp"),
  n_genes = length(genes_for_trimming)
)

genes_trimmed <- trim_gene_body(
  gr = genes_for_trimming,
  tss_trim = tss_trim,
  tes_trim = tes_trim
)

summary_table[[length(summary_table) + 1]] <- append_summary(
  step = paste0("trimmed_gene_bodies_TSS", tss_trim, "_TES", tes_trim),
  n_genes = length(genes_trimmed),
  output_file = trimmed_gtf
)

message("Exporting trimmed gene-body annotation: ", trimmed_gtf)
rtracklayer::export(
  genes_trimmed,
  trimmed_gtf,
  format = "gtf"
)

# ==============================================================================
# 6. Save filtering summary
# ==============================================================================

summary_table <- do.call(rbind, summary_table)

write.table(
  summary_table,
  file = summary_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Annotation filtering summary saved in: ", summary_file)
message("Done.")
