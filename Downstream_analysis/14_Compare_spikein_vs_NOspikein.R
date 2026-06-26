# ==============================================================================
# Compare RNAPII ChIP-seq differential occupancy with and without spike-in scaling
# ==============================================================================
# This independent script evaluates how the use of spike-in-derived size factors
# affects the detection of CPT-induced changes in RNAPII gene-body occupancy.
# It uses the same input files as the main DESeq2 ChIP-seq analysis, but writes
# all tables and figures to a separate output directory to avoid overwriting
# previous results.
#
# Two analyses are compared:
#   1) DESeq2 internal size-factor normalization
#   2) DESeq2 using externally derived spike-in size factors
#
# The script generates a category count barplot and a transition matrix showing
# how gene classifications change between both normalization strategies.
# Positive log2 fold changes indicate increased RNAPII occupancy in CPT.
# ==============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(ggplot2)
  library(tibble)
  library(cowplot)
})

# ------------------------------------------------------------------------------
# 1. Input and independent output paths
# ------------------------------------------------------------------------------

count_file <- file.path(
  "../Script_Limpios",
  "Count_Matrix_ChIPseq",
  "RNAPII_internal_gene_body",
  "RNAPII_ChIPseq_internal_gene_body_counts.tsv"
)

spike_factor_file <- file.path(
  "Data_final_CPT",
  "spikein_counts",
  "scale_factors.tsv"
)

# IMPORTANT: this is a new folder and does not overwrite the original DESeq2 output
output_dir <- file.path(
  "../Script_Limpios",
  "DESeq2_ChIPseq",
  "RNAPII_internal_gene_body_spikein_vs_internal_two_plots_v3"
)

table_dir <- file.path(output_dir, "tables")
plot_dir <- file.path(output_dir, "plots")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 2. Import count matrix
# ------------------------------------------------------------------------------

if (!file.exists(count_file)) {
  stop("The RNAPII ChIP-seq count matrix was not found: ", count_file,
       call. = FALSE)
}

if (!file.exists(spike_factor_file)) {
  stop("The spike-in size-factor file was not found: ", spike_factor_file,
       call. = FALSE)
}

count_matrix <- read.delim(
  file = count_file,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

count_matrix <- as.matrix(count_matrix)
storage.mode(count_matrix) <- "numeric"

if (anyNA(count_matrix)) {
  stop("The count matrix contains missing or non-numeric values.", call. = FALSE)
}

if (any(count_matrix < 0)) {
  stop("The count matrix contains negative values.", call. = FALSE)
}

count_matrix <- round(count_matrix)

expected_samples <- c("DMSO_R1", "CPT_R1", "DMSO_R2", "CPT_R2")

if (!setequal(colnames(count_matrix), expected_samples)) {
  stop(
    "The count-matrix columns do not match the expected samples: ",
    paste(expected_samples, collapse = ", "),
    call. = FALSE
  )
}

count_matrix <- count_matrix[, expected_samples, drop = FALSE]

# ------------------------------------------------------------------------------
# 3. Sample metadata
# ------------------------------------------------------------------------------

sample_metadata <- data.frame(
  sample = expected_samples,
  treatment = factor(
    c("DMSO", "CPT", "DMSO", "CPT"),
    levels = c("DMSO", "CPT")
  ),
  replicate = factor(c("R1", "R1", "R2", "R2")),
  row.names = expected_samples
)

# ------------------------------------------------------------------------------
# 4. Import and match spike-in size factors
# ------------------------------------------------------------------------------

spike_factors <- read.delim(
  file = spike_factor_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

required_spike_columns <- c("sample", "deseq2_size_factor")

if (!all(required_spike_columns %in% colnames(spike_factors))) {
  stop(
    "The spike-in file must contain the columns: ",
    paste(required_spike_columns, collapse = ", "),
    call. = FALSE
  )
}

spike_factors <- spike_factors %>%
  mutate(
    sample_clean = case_when(
      sample == "R1_RNAPII_AV_DM_S1"   ~ "DMSO_R1",
      sample == "R1_RNAPII_AV_CPT_S4"  ~ "CPT_R1",
      sample == "R2_RNAPII_AV_DM_S7"   ~ "DMSO_R2",
      sample == "R2_RNAPII_AV_CPT_S10" ~ "CPT_R2",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(sample_clean)) %>%
  arrange(match(sample_clean, expected_samples))

if (!identical(spike_factors$sample_clean, expected_samples)) {
  stop(
    "The spike-in samples could not be matched uniquely to the count matrix.",
    call. = FALSE
  )
}

if (
  anyNA(spike_factors$deseq2_size_factor) ||
  any(spike_factors$deseq2_size_factor <= 0)
) {
  stop("Spike-in DESeq2 size factors must be positive and non-missing.",
       call. = FALSE)
}

spike_size_factors <- spike_factors$deseq2_size_factor
names(spike_size_factors) <- spike_factors$sample_clean

# ------------------------------------------------------------------------------
# 5. Construct and filter DESeq2 datasets
# ------------------------------------------------------------------------------

dds_base <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_metadata,
  design = ~ treatment
)

genes_before_filtering <- nrow(dds_base)
keep <- rowSums(counts(dds_base)) >= 100
dds_base <- dds_base[keep, ]
genes_after_filtering <- nrow(dds_base)

if (genes_after_filtering == 0L) {
  stop("No genes remained after applying the total-count filter.", call. = FALSE)
}

filtering_summary <- data.frame(
  metric = c(
    "Genes before low-count filtering",
    "Genes retained with total counts >= 100",
    "Genes removed with total counts < 100"
  ),
  value = c(
    genes_before_filtering,
    genes_after_filtering,
    genes_before_filtering - genes_after_filtering
  )
)

write.table(
  filtering_summary,
  file = file.path(table_dir, "filtering_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Two independent copies: one internally normalized, one spike-in normalized
dds_internal <- dds_base
dds_spikein <- dds_base

# Internal DESeq2 normalization is estimated automatically by DESeq()
dds_internal <- DESeq(dds_internal)

# Spike-in normalization uses the externally derived size factors
deseq2_spike_size_factors <- spike_size_factors[colnames(dds_spikein)]
sizeFactors(dds_spikein) <- deseq2_spike_size_factors
dds_spikein <- DESeq(dds_spikein)

# Save size factors used by each approach
size_factor_table <- data.frame(
  sample = expected_samples,
  treatment = sample_metadata$treatment,
  replicate = sample_metadata$replicate,
  internal_DESeq2_size_factor = sizeFactors(dds_internal)[expected_samples],
  spikein_DESeq2_size_factor = sizeFactors(dds_spikein)[expected_samples]
)

write.table(
  size_factor_table,
  file = file.path(table_dir, "size_factors_internal_vs_spikein.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 6. Extract results from both analyses
# ------------------------------------------------------------------------------

classify_result <- function(padj, log2fc) {
  case_when(
    is.na(padj) ~ "padj_NA",
    padj <= 0.05 & log2fc > 1 ~ "Up_in_CPT",
    padj <= 0.05 & log2fc < -1 ~ "Down_in_CPT",
    TRUE ~ "Not_changed"
  )
}

res_internal <- results(
  dds_internal,
  contrast = c("treatment", "CPT", "DMSO"),
  alpha = 0.05
)

res_spikein <- results(
  dds_spikein,
  contrast = c("treatment", "CPT", "DMSO"),
  alpha = 0.05
)

internal_table <- as.data.frame(res_internal) %>%
  rownames_to_column("gene_id") %>%
  mutate(significance_internal = classify_result(padj, log2FoldChange)) %>%
  rename(
    baseMean_internal = baseMean,
    log2FoldChange_internal = log2FoldChange,
    lfcSE_internal = lfcSE,
    stat_internal = stat,
    pvalue_internal = pvalue,
    padj_internal = padj
  )

spikein_table <- as.data.frame(res_spikein) %>%
  rownames_to_column("gene_id") %>%
  mutate(significance_spikein = classify_result(padj, log2FoldChange)) %>%
  rename(
    baseMean_spikein = baseMean,
    log2FoldChange_spikein = log2FoldChange,
    lfcSE_spikein = lfcSE,
    stat_spikein = stat,
    pvalue_spikein = pvalue,
    padj_spikein = padj
  )

comparison_table <- internal_table %>%
  inner_join(spikein_table, by = "gene_id") %>%
  mutate(
    delta_log2FC_spikein_minus_internal =
      log2FoldChange_spikein - log2FoldChange_internal,
    category_change = paste(
      significance_internal,
      significance_spikein,
      sep = "_to_"
    )
  )

write.table(
  internal_table,
  file = file.path(table_dir, "DESeq2_internal_normalization_all_genes.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  spikein_table,
  file = file.path(table_dir, "DESeq2_spikein_normalization_all_genes.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  comparison_table,
  file = file.path(table_dir, "DESeq2_internal_vs_spikein_comparison_all_genes.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 7. Summary tables
# ------------------------------------------------------------------------------

plot_breaks <- c("Up_in_CPT", "Down_in_CPT", "Not_changed")

category_counts_long <- bind_rows(
  comparison_table %>%
    filter(significance_internal %in% plot_breaks) %>%
    dplyr::count(significance_internal, name = "n") %>%
    rename(significance = significance_internal) %>%
    mutate(normalization = "Internal DESeq2"),
  comparison_table %>%
    filter(significance_spikein %in% plot_breaks) %>%
    dplyr::count(significance_spikein, name = "n") %>%
    rename(significance = significance_spikein) %>%
    mutate(normalization = "Spike-in")
) %>%
  mutate(
    significance = factor(significance, levels = plot_breaks),
    normalization = factor(normalization, levels = c("Internal DESeq2", "Spike-in")),
    normalization_label = factor(
      recode(
        as.character(normalization),
        "Internal DESeq2" = "No spike-in",
        "Spike-in" = "Spike-in"
      ),
      levels = c("No spike-in", "Spike-in")
    )
  )

category_cross_table <- comparison_table %>%
  dplyr::count(significance_internal, significance_spikein, name = "n") %>%
  arrange(significance_internal, significance_spikein)

summary_metrics <- data.frame(
  metric = c(
    "Genes analysed after filtering",
    "Up in CPT with internal DESeq2 normalization",
    "Down in CPT with internal DESeq2 normalization",
    "Up in CPT with spike-in normalization",
    "Down in CPT with spike-in normalization",
    "Genes Up only with spike-in normalization",
    "Genes Down only with spike-in normalization",
    "Median delta log2FC spike-in minus internal"
  ),
  value = c(
    nrow(comparison_table),
    sum(comparison_table$significance_internal == "Up_in_CPT", na.rm = TRUE),
    sum(comparison_table$significance_internal == "Down_in_CPT", na.rm = TRUE),
    sum(comparison_table$significance_spikein == "Up_in_CPT", na.rm = TRUE),
    sum(comparison_table$significance_spikein == "Down_in_CPT", na.rm = TRUE),
    sum(comparison_table$significance_internal != "Up_in_CPT" &
          comparison_table$significance_spikein == "Up_in_CPT", na.rm = TRUE),
    sum(comparison_table$significance_internal != "Down_in_CPT" &
          comparison_table$significance_spikein == "Down_in_CPT", na.rm = TRUE),
    median(comparison_table$delta_log2FC_spikein_minus_internal, na.rm = TRUE)
  )
)

write.table(
  category_counts_long,
  file = file.path(table_dir, "category_counts_internal_vs_spikein_long.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  category_cross_table,
  file = file.path(table_dir, "category_cross_table_internal_vs_spikein.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  summary_metrics,
  file = file.path(table_dir, "summary_metrics_internal_vs_spikein.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Gene lists useful for inspection in IGV or downstream checks
write.table(
  comparison_table %>%
    filter(significance_internal != "Up_in_CPT", significance_spikein == "Up_in_CPT") %>%
    arrange(padj_spikein, desc(log2FoldChange_spikein)),
  file = file.path(table_dir, "genes_up_only_with_spikein.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  comparison_table %>%
    filter(significance_internal == "Up_in_CPT", significance_spikein != "Up_in_CPT") %>%
    arrange(padj_internal, desc(log2FoldChange_internal)),
  file = file.path(table_dir, "genes_up_only_with_internal.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 8. Plot settings
# ------------------------------------------------------------------------------

occupancy_colours <- c(
  "Up_in_CPT" = "#C44E8C",
  "Down_in_CPT" = "#2C7FB8",
  "Not_changed" = "#D9D9D9",
  "padj_NA" = "#D9D9D9"
)

plot_labels <- c(
  "Up_in_CPT" = "Increased\n in CPT",
  "Down_in_CPT" = "Decreased\n in CPT",
  "Not_changed" = "Not changed"
)

# High-resolution PNG settings.
# The chosen plot sizes at 1200 dpi stay below the 100-megapixel limit.
png_dpi <- 1200

save_plot <- function(plot, filename, width, height) {
  megapixels <- (width * png_dpi) * (height * png_dpi) / 1e6

  if (megapixels > 100) {
    stop(
      "Requested PNG size exceeds 100 megapixels: ",
      round(megapixels, 1), " MP. Reduce width, height or dpi.",
      call. = FALSE
    )
  }

  ggsave(
    filename = file.path(plot_dir, paste0(filename, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    device = cairo_pdf,
    bg = "white"
  )

  ggsave(
    filename = file.path(plot_dir, paste0(filename, ".png")),
    plot = plot,
    width = width,
    height = height,
    dpi = png_dpi,
    limitsize = FALSE,
    bg = "white"
  )
}

plot_theme <- function(show_legend = FALSE) {
  theme_classic(base_size = 20) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      axis.title = element_text(size = 22),
      axis.text = element_text(size = 20, colour = "black"),
      axis.line = element_line(linewidth = 0.7),
      axis.ticks = element_line(linewidth = 0.7),
      axis.ticks.length = unit(0.22, "cm"),
      legend.position = if (show_legend) "right" else "none",
      legend.text = element_text(size = 18),
      legend.title = element_blank(),
      legend.key.size = unit(0.7, "cm"),
      plot.margin = margin(12, 16, 12, 12)
    )
}

# ------------------------------------------------------------------------------
# 9. Plot: category counts with/without spike-in
# ------------------------------------------------------------------------------

category_counts_plot <- ggplot(
  category_counts_long,
  aes(x = normalization_label, y = n, fill = significance)
) +
  geom_col(colour = "white", linewidth = 0.7, width = 0.72) +
  geom_text(
    aes(label = n),
    position = position_stack(vjust = 0.5),
    size = 6,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = occupancy_colours,
    breaks = plot_breaks,
    labels = plot_labels,
    drop = FALSE
  ) +
  labs(
    x = NULL,
    y = "Number of genes",
    fill = NULL
  ) +
  plot_theme(show_legend = FALSE) +
  theme(
    axis.text.x = element_text(size = 18, colour = "black", angle = 0, hjust = 0.5),
    axis.title.y = element_text(size = 22),
    legend.position = "none",
    plot.margin = margin(12, 16, 12, 16)
  )

save_plot(
  category_counts_plot,
  "category_counts_internal_vs_spikein",
  width = 6.2,
  height = 5.0
)

# ------------------------------------------------------------------------------
# 10. Plot: category transition heatmap
# ------------------------------------------------------------------------------

transition_data <- category_cross_table %>%
  filter(
    significance_internal %in% plot_breaks,
    significance_spikein %in% plot_breaks
  )

# Add zero-count combinations so that absent transitions are shown explicitly as 0
transition_grid <- expand.grid(
  significance_internal = plot_breaks,
  significance_spikein = plot_breaks,
  stringsAsFactors = FALSE
)

transition_data <- transition_grid %>%
  left_join(transition_data, by = c("significance_internal", "significance_spikein")) %>%
  mutate(
    n = ifelse(is.na(n), 0L, n),
    significance_internal = factor(significance_internal, levels = plot_breaks),
    significance_spikein = factor(significance_spikein, levels = plot_breaks),
    label_colour = ifelse(n > max(n, na.rm = TRUE) * 0.45, "white", "black")
  )

transition_heatmap <- ggplot(
  transition_data,
  aes(x = significance_internal, y = significance_spikein, fill = n)
) +
  geom_tile(colour = "grey85", linewidth = 0.8) +
  geom_text(aes(label = n, colour = label_colour), size = 6, fontface = "bold") +
  scale_colour_identity() +
  scale_fill_gradient(low = "white", high = "#08306B") +
  scale_x_discrete(labels = plot_labels) +
  scale_y_discrete(labels = plot_labels) +
  labs(
    x = "No spike-in normalization",
    y = "Spike-in normalization",
    fill = "Genes"
  ) +
  theme_classic(base_size = 18) +
  theme(
    axis.title = element_text(size = 22),
    axis.text.x = element_text(size = 16, angle = 30, hjust = 1, colour = "black"),
    axis.text.y = element_text(size = 16, colour = "black"),
    legend.position = "right",
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    plot.margin = margin(12, 16, 12, 12)
  )

save_plot(
  transition_heatmap,
  "category_transition_internal_vs_spikein_heatmap",
  width = 7,
  height = 6
)
