# ==============================================================================
# Define the common active-gene universe for RNA-seq and RNAPII ChIP-seq
# ==============================================================================
# This script defines the gene universe used for downstream correspondence
# analyses between RNA-seq expression changes and RNAPII ChIP-seq occupancy.
#
# A gene is retained when it:
#   1. passed the low-count filter in both DESeq2 analyses;
#   2. belongs to the protein-coding, minimum-3-kb annotation used for ChIP-seq;
#   3. has mean DMSO RNA-seq expression >= 1 RPKM;
#   4. has mean DMSO RNAPII ChIP-seq signal above a matched percentile.
#
# The ChIP-seq percentile is calculated from the fraction of comparable RNA-seq
# genes with mean DMSO RPKM < 1. The same lower fraction is then excluded from
# the basal RNAPII ChIP-seq signal distribution.
#
# Main downstream output:
#   common_active_genes_DMSO_RNAseq_RNAPII_ChIPseq.tsv
# ==============================================================================

library(dplyr)
library(ggplot2)
library(tibble)
library(stringr)
library(ggVennDiagram)

# ------------------------------------------------------------------------------
# 1. Input and output paths
# ------------------------------------------------------------------------------

rna_featurecounts_file <- file.path(
  "../Datos_Bulk_RNAseq",
  "counts",
  "exon",
  "featureCounts_exon.txt"
)

rna_deseq2_results_file <- file.path(
  "../Script_Limpios",
  "DESeq2_RNAseq",
  "tables",
  "RNAseq_DESeq2_CPT_vs_DMSO_all_genes.tsv"
)

chip_normalized_counts_file <- file.path(
  "../Script_Limpios",
  "DESeq2_ChIPseq",
  "RNAPII_internal_gene_body",
  "tables",
  "RNAPII_ChIPseq_spikein_normalized_counts.tsv"
)

chip_deseq2_results_file <- file.path(
  "../Script_Limpios",
  "DESeq2_ChIPseq",
  "RNAPII_internal_gene_body",
  "tables",
  "RNAPII_ChIPseq_DESeq2_CPT_vs_DMSO_all_genes.tsv"
)

output_dir <- file.path(
  "../Script_Limpios",
  "RNAseq_ChIPseq_correspondence",
  "common_gene_universe"
)

table_dir <- file.path(output_dir, "tables")
plot_dir <- file.path(output_dir, "plots")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

input_files <- c(
  rna_featurecounts_file,
  rna_deseq2_results_file,
  chip_normalized_counts_file,
  chip_deseq2_results_file
)

missing_files <- input_files[!file.exists(input_files)]

if (length(missing_files) > 0L) {
  stop(
    "The following required input files were not found:\n",
    paste(missing_files, collapse = "\n"),
    call. = FALSE
  )
}

# ------------------------------------------------------------------------------
# 2. Analysis thresholds
# ------------------------------------------------------------------------------

rna_active_rpkm_threshold <- 1

# The ChIP-seq quantile is calculated later from the RNA-seq distribution after
# both datasets have been restricted to the comparable annotation universe.

# ------------------------------------------------------------------------------
# 3. Helper functions
# ------------------------------------------------------------------------------

clean_gene_id <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_remove("^hs_") %>%
    stringr::str_remove("\\.[0-9]+$")
}

# PNG export resolution.
# The DPI is calculated automatically to maximise image resolution while keeping
# each PNG below the 100-megapixel limit of the downstream figure-assembly software.
maximum_png_pixels <- 100000000

safe_png_dpi <- function(width, height, maximum_pixels = maximum_png_pixels) {
  floor(
    sqrt(
      maximum_pixels / (width * height)
    )
  )
}

save_plot <- function(plot, filename, width, height) {
  ggsave(
    filename = file.path(plot_dir, paste0(filename, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    bg = "white"
  )
  
  ggsave(
    filename = file.path(plot_dir, paste0(filename, ".png")),
    plot = plot,
    width = width,
    height = height,
    dpi = safe_png_dpi(width, height),
    bg = "white"
  )
}

# ------------------------------------------------------------------------------
# 4. Import exon-level RNA-seq counts and calculate DMSO RPKM
# ------------------------------------------------------------------------------

rna_featurecounts <- read.delim(
  file = rna_featurecounts_file,
  header = TRUE,
  skip = 1,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_annotation_columns <- c(
  "Geneid", "Chr", "Start", "End", "Strand", "Length"
)

if (!all(required_annotation_columns %in% colnames(rna_featurecounts))) {
  stop(
    "The RNA-seq featureCounts file does not contain the expected annotation ",
    "columns.",
    call. = FALSE
  )
}

rna_sample_columns <- setdiff(
  colnames(rna_featurecounts),
  required_annotation_columns
)

if (length(rna_sample_columns) == 0L) {
  stop("No RNA-seq sample columns were detected.", call. = FALSE)
}

rna_original_sample_names <- basename(rna_sample_columns)
rna_original_sample_names <- sub("\\.bam$", "", rna_original_sample_names)

rna_sample_name_map <- c(
  "DM_rep1" = "DMSO_R1",
  "DM_rep2" = "DMSO_R2",
  "DM_rep3" = "DMSO_R3",
  "DM_rep4" = "DMSO_R4",
  "CPT_rep1" = "CPT_R1",
  "CPT_rep2" = "CPT_R2",
  "CPT_rep3" = "CPT_R3",
  "CPT_rep4" = "CPT_R4"
)

rna_sample_names <- unname(
  rna_sample_name_map[rna_original_sample_names]
)

if (anyNA(rna_sample_names)) {
  stop(
    "RNA-seq sample names could not be assigned for: ",
    paste(rna_original_sample_names[is.na(rna_sample_names)], collapse = ", "),
    call. = FALSE
  )
}

rna_counts <- as.matrix(
  rna_featurecounts[, rna_sample_columns, drop = FALSE]
)

storage.mode(rna_counts) <- "numeric"
colnames(rna_counts) <- rna_sample_names
rownames(rna_counts) <- clean_gene_id(rna_featurecounts$Geneid)

if (anyNA(rna_counts) || any(rna_counts < 0)) {
  stop("The RNA-seq count matrix contains invalid values.", call. = FALSE)
}

rna_gene_length_kb <- rna_featurecounts$Length / 1000
names(rna_gene_length_kb) <- clean_gene_id(rna_featurecounts$Geneid)

valid_length <- (
  !is.na(rna_gene_length_kb) &
  rna_gene_length_kb > 0
)

rna_gene_length_kb <- rna_gene_length_kb[valid_length]
rna_counts <- rna_counts[
  rownames(rna_counts) %in% names(rna_gene_length_kb),
  ,
  drop = FALSE
]
rna_gene_length_kb <- rna_gene_length_kb[rownames(rna_counts)]

rna_library_size_millions <- colSums(rna_counts) / 1e6

rna_rpkm <- sweep(
  rna_counts,
  2,
  rna_library_size_millions,
  FUN = "/"
)

rna_rpkm <- sweep(
  rna_rpkm,
  1,
  rna_gene_length_kb,
  FUN = "/"
)

rna_dm_samples <- paste0("DMSO_R", 1:4)

if (!all(rna_dm_samples %in% colnames(rna_rpkm))) {
  stop(
    "Not all four expected DMSO RNA-seq samples were found.",
    call. = FALSE
  )
}

rna_mean_dm_rpkm <- rowMeans(
  rna_rpkm[, rna_dm_samples, drop = FALSE],
  na.rm = TRUE
)

# Restrict RNA-seq genes to those retained after the DESeq2 total-count filter.
rna_deseq2_results <- read.delim(
  rna_deseq2_results_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (!"gene_id" %in% colnames(rna_deseq2_results)) {
  stop("The RNA-seq DESeq2 table does not contain gene_id.", call. = FALSE)
}

rna_deseq2_gene_ids <- unique(
  clean_gene_id(rna_deseq2_results$gene_id)
)

rna_mean_dm_rpkm <- rna_mean_dm_rpkm[
  names(rna_mean_dm_rpkm) %in% rna_deseq2_gene_ids
]

# ------------------------------------------------------------------------------
# 5. Import spike-in-normalized RNAPII ChIP-seq counts
# ------------------------------------------------------------------------------

chip_normalized_counts <- read.delim(
  chip_normalized_counts_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_chip_columns <- c(
  "gene_id", "DMSO_R1", "CPT_R1", "DMSO_R2", "CPT_R2"
)

if (!all(required_chip_columns %in% colnames(chip_normalized_counts))) {
  stop(
    "The RNAPII ChIP-seq normalized-count table does not contain the expected ",
    "columns.",
    call. = FALSE
  )
}

chip_gene_ids <- clean_gene_id(chip_normalized_counts$gene_id)

chip_count_matrix <- as.matrix(
  chip_normalized_counts[
    ,
    c("DMSO_R1", "CPT_R1", "DMSO_R2", "CPT_R2"),
    drop = FALSE
  ]
)

storage.mode(chip_count_matrix) <- "numeric"
rownames(chip_count_matrix) <- chip_gene_ids

if (anyNA(chip_count_matrix) || any(chip_count_matrix < 0)) {
  stop(
    "The RNAPII ChIP-seq normalized-count matrix contains invalid values.",
    call. = FALSE
  )
}

chip_mean_dm_normalized_count <- rowMeans(
  chip_count_matrix[, c("DMSO_R1", "DMSO_R2"), drop = FALSE],
  na.rm = TRUE
)

# Confirm that the ChIP-seq genes correspond to those retained after the DESeq2
# total-count filter.
chip_deseq2_results <- read.delim(
  chip_deseq2_results_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (!"gene_id" %in% colnames(chip_deseq2_results)) {
  stop("The ChIP-seq DESeq2 table does not contain gene_id.", call. = FALSE)
}

chip_deseq2_gene_ids <- unique(
  clean_gene_id(chip_deseq2_results$gene_id)
)

chip_mean_dm_normalized_count <- chip_mean_dm_normalized_count[
  names(chip_mean_dm_normalized_count) %in% chip_deseq2_gene_ids
]

# The ChIP-seq count matrix was generated from the final protein-coding,
# minimum-3-kb annotation. It therefore defines the eligible annotation universe.
annotation_gene_universe <- names(chip_mean_dm_normalized_count)

rna_mean_dm_rpkm <- rna_mean_dm_rpkm[
  names(rna_mean_dm_rpkm) %in% annotation_gene_universe
]

# ------------------------------------------------------------------------------
# 6. Define active genes in DMSO
# ------------------------------------------------------------------------------

# Calculate the fraction of comparable RNA-seq genes below 1 RPKM. This fraction
# defines the matched lower percentile excluded from the ChIP-seq distribution.
rna_fraction_below_active_threshold <- mean(
  rna_mean_dm_rpkm < rna_active_rpkm_threshold,
  na.rm = TRUE
)

chip_active_quantile <- rna_fraction_below_active_threshold

active_rna_genes <- names(rna_mean_dm_rpkm)[
  rna_mean_dm_rpkm >= rna_active_rpkm_threshold
]

chip_active_threshold <- as.numeric(
  quantile(
    chip_mean_dm_normalized_count,
    probs = chip_active_quantile,
    na.rm = TRUE,
    names = FALSE
  )
)

active_chip_genes <- names(chip_mean_dm_normalized_count)[
  chip_mean_dm_normalized_count >= chip_active_threshold
]

common_active_genes <- intersect(
  active_rna_genes,
  active_chip_genes
)

matched_quantile_label <- paste0(
  "Q",
  format(
    round(chip_active_quantile * 100, 1),
    trim = TRUE,
    nsmall = 1
  )
)

# ------------------------------------------------------------------------------
# 7. Export threshold and filtering summaries
# ------------------------------------------------------------------------------

filtering_summary <- data.frame(
  metric = c(
    "RNA-seq genes retained after DESeq2 low-count filtering",
    "ChIP-seq genes retained after DESeq2 low-count filtering",
    "Eligible protein-coding minimum-3-kb genes represented in ChIP-seq",
    "RNA-seq active genes in DMSO",
    "RNAPII ChIP-seq active genes in DMSO",
    "Genes active in both datasets in DMSO"
  ),
  value = c(
    length(rna_deseq2_gene_ids),
    length(chip_deseq2_gene_ids),
    length(annotation_gene_universe),
    length(active_rna_genes),
    length(active_chip_genes),
    length(common_active_genes)
  )
)

write.table(
  filtering_summary,
  file = file.path(table_dir, "common_gene_universe_filtering_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

threshold_summary <- data.frame(
  dataset = c("RNA-seq", "RNAPII ChIP-seq"),
  metric = c("Mean DMSO RPKM", "Mean DMSO spike-in-normalized count"),
  threshold_type = c(
    paste0("RPKM >= ", rna_active_rpkm_threshold),
    paste0(
      matched_quantile_label,
      " matched to the RNA-seq fraction below 1 RPKM"
    )
  ),
  threshold_value = c(
    rna_active_rpkm_threshold,
    chip_active_threshold
  ),
  fraction_removed = c(
    rna_fraction_below_active_threshold,
    chip_active_quantile
  )
)

write.table(
  threshold_summary,
  file = file.path(table_dir, "common_gene_universe_thresholds.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

quantile_probabilities <- sort(
  unique(
    c(
      0,
      0.10,
      chip_active_quantile,
      0.25,
      0.50,
      0.75,
      0.90
    )
  )
)

chip_quantile_summary <- data.frame(
  quantile = paste0(
    "Q",
    format(
      round(quantile_probabilities * 100, 1),
      trim = TRUE,
      nsmall = 1
    )
  ),
  probability = quantile_probabilities,
  threshold_value = as.numeric(
    quantile(
      chip_mean_dm_normalized_count,
      probs = quantile_probabilities,
      na.rm = TRUE
    )
  )
) %>%
  dplyr::mutate(
    genes_retained = vapply(
      threshold_value,
      function(x) {
        sum(chip_mean_dm_normalized_count >= x, na.rm = TRUE)
      },
      integer(1)
    ),
    fraction_retained = genes_retained /
      length(chip_mean_dm_normalized_count)
  )

write.table(
  chip_quantile_summary,
  file = file.path(table_dir, "RNAPII_ChIPseq_DMSO_quantile_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 8. Export active-gene tables
# ------------------------------------------------------------------------------

rna_active_table <- tibble::tibble(
  gene_id = active_rna_genes,
  RNAseq_mean_DMSO_RPKM = unname(
    rna_mean_dm_rpkm[active_rna_genes]
  )
) %>%
  dplyr::arrange(gene_id)

chip_active_table <- tibble::tibble(
  gene_id = active_chip_genes,
  RNAPII_ChIPseq_mean_DMSO_normalized_count = unname(
    chip_mean_dm_normalized_count[active_chip_genes]
  )
) %>%
  dplyr::arrange(gene_id)

common_active_table <- tibble::tibble(
  gene_id = common_active_genes,
  RNAseq_mean_DMSO_RPKM = unname(
    rna_mean_dm_rpkm[common_active_genes]
  ),
  RNAPII_ChIPseq_mean_DMSO_normalized_count = unname(
    chip_mean_dm_normalized_count[common_active_genes]
  )
) %>%
  dplyr::arrange(gene_id)

write.table(
  rna_active_table,
  file = file.path(table_dir, "active_genes_DMSO_RNAseq.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  chip_active_table,
  file = file.path(table_dir, "active_genes_DMSO_RNAPII_ChIPseq.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  common_active_table,
  file = file.path(
    table_dir,
    "common_active_genes_DMSO_RNAseq_RNAPII_ChIPseq.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  common_active_genes,
  file = file.path(
    table_dir,
    "common_active_genes_DMSO_gene_list.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# ------------------------------------------------------------------------------
# 9. Diagnostic histograms
# ------------------------------------------------------------------------------

rna_histogram <- ggplot(
  data.frame(mean_DMSO_RPKM = rna_mean_dm_rpkm),
  aes(x = log10(mean_DMSO_RPKM + 1))
) +
  geom_histogram(bins = 100, fill = "grey40") +
  geom_vline(
    xintercept = log10(rna_active_rpkm_threshold + 1),
    linetype = "dashed",
    linewidth = 0.6
  ) +
  labs(
    x = "log10(mean DMSO RPKM + 1)",
    y = "Number of genes"
  ) +
  theme_classic(base_size = 20)

chip_histogram <- ggplot(
  data.frame(
    mean_DMSO_normalized_count = chip_mean_dm_normalized_count
  ),
  aes(x = log10(mean_DMSO_normalized_count + 1))
) +
  geom_histogram(bins = 100, fill = "grey40") +
  geom_vline(
    xintercept = log10(chip_active_threshold + 1),
    linetype = "dashed",
    linewidth = 0.6
  ) +
  labs(
    x = "log10(mean DMSO spike-in-normalized count + 1)",
    y = "Number of genes"
  ) +
  theme_classic(base_size = 20)

save_plot(
  rna_histogram,
  "RNAseq_mean_DMSO_RPKM_distribution",
  width = 7.5,
  height = 5
)

save_plot(
  chip_histogram,
  "RNAPII_ChIPseq_mean_DMSO_normalized_count_distribution_matched_percentile",
  width = 7.5,
  height = 5
)

# ------------------------------------------------------------------------------
# 10. Venn diagram of active genes
# ------------------------------------------------------------------------------

venn_data <- list(
  "RNA-seq" = active_rna_genes,
  "RNAPII ChIP-seq" = active_chip_genes
)

# Manual set-label positions retained from the original script.
label_positions <- tibble::tibble(
  label = c("RNA-seq", "ChIP-seq"),
  x = c(-1.5, -1.5),
  y = c(-3, 7)
)

venn_plot <- ggVennDiagram(
  venn_data,
  label_alpha = 0,
  label_size = 7,
  set_size = 0
) +
  scale_fill_gradient(
    low = "#FFF7BC",
    high = "#EA580C",
    name = "Gene count"
  ) +
  coord_equal(clip = "off") +
  geom_text(
    data = label_positions,
    aes(
      x = x,
      y = y,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 0,
    fontface = "bold",
    size = 7
  ) +
  theme_void(base_size = 14) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 12,
      margin = margin(b = 12)
    ),
    plot.margin = margin(40, 30, 30, 90),
    legend.position = "right"
  )

ggsave(
  filename = file.path(
    plot_dir,
    "common_active_genes_DMSO_venn.pdf"
  ),
  plot = venn_plot,
  width = 8,
  height = 6,
  bg = "white"
)

ggsave(
  filename = file.path(
    plot_dir,
    "common_active_genes_DMSO_venn.png"
  ),
  plot = venn_plot,
  width = 8,
  height = 6,
  dpi = safe_png_dpi(8, 6),
  bg = "white"
)