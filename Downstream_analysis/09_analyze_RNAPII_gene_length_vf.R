# ==============================================================================
# Analyse the relationship between RNAPII occupancy changes and gene length
# ==============================================================================
# This script combines differential RNAPII ChIP-seq results with internal
# gene-body length information and generates the final gene-length analyses.
#
# Input data:
#   - RNAPII ChIP-seq DESeq2 results for CPT versus DMSO
#   - featureCounts annotation table for strand-aware internal gene-body regions
#
# Analysis-specific filter:
#   - internal gene-body regions shorter than 1 kb are excluded
#
# Main analyses:
#   - descriptive gene-length statistics by RNAPII response category
#   - Shapiro-Wilk normality diagnostics within each comparison group
#   - Kruskal-Wallis and pairwise Wilcoxon tests
#   - internal gene-body length distributions by RNAPII category
#   - RNAPII log2FC distributions in the shortest and longest 10% of genes
#   - ranked RNAPII log2FC plot coloured by gene-length group
#   - RNAPII category composition across gene-length deciles
#   - sensitivity analysis across progressively stronger RNAPII-gain genes
#
# ==============================================================================

library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(cowplot)

# ------------------------------------------------------------------------------
# 1. Input and output paths
# ------------------------------------------------------------------------------

chip_deseq2_results_file <- file.path(
  "../Script_Limpios",
  "DESeq2_ChIPseq",
  "RNAPII_internal_gene_body",
  "tables",
  "RNAPII_ChIPseq_DESeq2_CPT_vs_DMSO_all_genes.tsv"
)

chip_annotation_file <- file.path(
  "../Script_Limpios",
  "Count_Matrix_ChIPseq",
  "RNAPII_internal_gene_body",
  "RNAPII_ChIPseq_internal_gene_body_counts_with_annotation.tsv"
)

output_dir <- file.path(
  "../Script_Limpios",
  "RNAPII_gene_length_analysis_vf"
)

table_dir <- file.path(output_dir, "tables")
plot_dir <- file.path(output_dir, "plots")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

input_files <- c(
  chip_deseq2_results_file,
  chip_annotation_file
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
# 2. Analysis settings
# ------------------------------------------------------------------------------

# Minimum length of the strand-aware trimmed internal gene-body region.
minimum_internal_gene_body_length_bp <- 1000

# Fraction used to define the shortest and longest gene groups.
extreme_length_fraction <- 0.10

# Shapiro-Wilk normality diagnostics.
# stats::shapiro.test() accepts a maximum of 5,000 observations. Groups larger
# than this limit are reproducibly subsampled without replacement.
shapiro_maximum_n <- 5000L
shapiro_random_seed <- 123L

# Thresholds used only for the final ranked RNAPII-gain sensitivity analysis.
rnapii_gain_log2fc_threshold <- 1
rnapii_gain_padj_threshold <- 0.05

# Nested fractions retained after ranking significant RNAPII-gain genes from
# highest to lowest log2FC.
rnapii_gain_retained_fractions <- c(
  1.00,
  0.90,
  0.75,
  0.50,
  0.25,
  0.10
)

rnapii_gain_retained_labels <- c(
  "All",
  "Top 90%",
  "Top 75%",
  "Top 50%",
  "Top 25%",
  "Top 10%"
)

# RNAPII differential-occupancy categories.
category_levels <- c(
  "Up_in_CPT",
  "Down_in_CPT",
  "Not_changed"
)

category_labels <- c(
  "Up_in_CPT" = "UP",
  "Down_in_CPT" = "DOWN",
  "Not_changed" = "NC"
)

category_plot_levels <- c(
  "UP",
  "DOWN",
  "NC"
)

category_colours <- c(
  "UP" = "#C44E8C",
  "DOWN" = "#2C7FB8",
  "NC" = "#D9D9D9"
)

# Extreme gene-length groups.
length_group_levels <- c(
  "Shortest 10%",
  "Middle 80%",
  "Longest 10%"
)

length_group_colours <- c(
  "Shortest 10%" = "#009E9A",
  "Middle 80%" = "#D9D9D9",
  "Longest 10%" = "#D89C00"
)

sensitivity_plot_colour <- "#3B102F"

# PNG export resolution.
plot_dpi <- 1200

# ------------------------------------------------------------------------------
# 3. Helper functions
# ------------------------------------------------------------------------------

clean_gene_id <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_remove("^hs_") %>%
    stringr::str_remove("\\.[0-9]+$")
}

save_plot <- function(plot, filename, width, height) {
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
    dpi = plot_dpi,
    bg = "white"
  )
}

save_plot_legend <- function(plot, output_base, width, height) {
  legend_only <- cowplot::get_legend(plot)
  legend_plot <- cowplot::ggdraw(legend_only)

  ggsave(
    filename = paste0(output_base, ".pdf"),
    plot = legend_plot,
    width = width,
    height = height,
    device = cairo_pdf,
    bg = "white"
  )

  ggsave(
    filename = paste0(output_base, ".png"),
    plot = legend_plot,
    width = width,
    height = height,
    dpi = plot_dpi,
    bg = "white"
  )
}

run_shapiro_test <- function(x, maximum_n = shapiro_maximum_n) {

  # Retain only finite observations because shapiro.test() does not accept
  # missing or infinite values.
  x <- x[is.finite(x)]
  n_total <- length(x)

  if (n_total < 3L) {
    stop(
      "Shapiro-Wilk test requires at least three finite observations.",
      call. = FALSE
    )
  }

  n_tested <- min(n_total, maximum_n)
  was_subsampled <- n_total > maximum_n

  if (was_subsampled) {
    x <- sample(
      x,
      size = maximum_n,
      replace = FALSE
    )
  }

  test_result <- shapiro.test(x)

  tibble::tibble(
    n_total = n_total,
    n_tested = n_tested,
    subsampled = was_subsampled,
    shapiro_W = unname(test_result$statistic),
    shapiro_p_value = test_result$p.value
  )
}

# ------------------------------------------------------------------------------
# 4. Import RNAPII ChIP-seq DESeq2 results
# ------------------------------------------------------------------------------

chip_deseq2_results <- read.delim(
  file = chip_deseq2_results_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_deseq2_columns <- c(
  "gene_id",
  "baseMean",
  "log2FoldChange",
  "padj",
  "significance"
)

if (!all(required_deseq2_columns %in% colnames(chip_deseq2_results))) {
  stop(
    "The RNAPII ChIP-seq DESeq2 table does not contain all required columns.",
    call. = FALSE
  )
}

chip_deseq2_results <- chip_deseq2_results %>%
  dplyr::mutate(
    gene_id = clean_gene_id(gene_id),
    significance = dplyr::if_else(
      significance == "NC",
      "Not_changed",
      significance
    )
  ) %>%
  dplyr::filter(
    significance %in% category_levels
  )

if (anyDuplicated(chip_deseq2_results$gene_id)) {
  stop(
    "Duplicated gene IDs were found in the RNAPII ChIP-seq DESeq2 table.",
    call. = FALSE
  )
}

# ------------------------------------------------------------------------------
# 5. Import internal gene-body lengths
# ------------------------------------------------------------------------------

chip_annotation <- read.delim(
  file = chip_annotation_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_annotation_columns <- c(
  "GeneID",
  "Length"
)

if (!all(required_annotation_columns %in% colnames(chip_annotation))) {
  stop(
    "The annotated featureCounts table does not contain GeneID and Length.",
    call. = FALSE
  )
}

gene_lengths <- chip_annotation %>%
  dplyr::transmute(
    gene_id = clean_gene_id(GeneID),
    internal_gene_body_length_bp = as.numeric(Length)
  ) %>%
  dplyr::filter(
    !is.na(internal_gene_body_length_bp),
    internal_gene_body_length_bp >= minimum_internal_gene_body_length_bp
  ) %>%
  dplyr::mutate(
    internal_gene_body_length_kb = internal_gene_body_length_bp / 1000
  )

if (anyDuplicated(gene_lengths$gene_id)) {
  stop(
    "Duplicated gene IDs were found in the gene-length annotation table.",
    call. = FALSE
  )
}

# ------------------------------------------------------------------------------
# 6. Build the master gene-length analysis table
# ------------------------------------------------------------------------------

gene_length_analysis <- chip_deseq2_results %>%
  dplyr::select(
    gene_id,
    baseMean,
    log2FoldChange,
    padj,
    significance
  ) %>%
  dplyr::inner_join(
    gene_lengths,
    by = "gene_id"
  ) %>%
  dplyr::mutate(
    category = factor(
      unname(category_labels[significance]),
      levels = category_plot_levels
    )
  ) %>%
  dplyr::arrange(
    category,
    gene_id
  )

if (nrow(gene_length_analysis) == 0L) {
  stop(
    "No genes remained after joining DESeq2 results with gene lengths.",
    call. = FALSE
  )
}

if (anyNA(gene_length_analysis$category)) {
  stop(
    "One or more RNAPII response categories could not be assigned.",
    call. = FALSE
  )
}

# ------------------------------------------------------------------------------
# 7. Define extreme gene-length groups
# ------------------------------------------------------------------------------

short_length_threshold <- as.numeric(
  quantile(
    gene_length_analysis$internal_gene_body_length_kb,
    probs = extreme_length_fraction,
    na.rm = TRUE,
    names = FALSE
  )
)

long_length_threshold <- as.numeric(
  quantile(
    gene_length_analysis$internal_gene_body_length_kb,
    probs = 1 - extreme_length_fraction,
    na.rm = TRUE,
    names = FALSE
  )
)

gene_length_analysis <- gene_length_analysis %>%
  dplyr::mutate(
    length_group = dplyr::case_when(
      internal_gene_body_length_kb <= short_length_threshold ~
        "Shortest 10%",
      internal_gene_body_length_kb >= long_length_threshold ~
        "Longest 10%",
      TRUE ~ "Middle 80%"
    ),
    length_group = factor(
      length_group,
      levels = length_group_levels
    )
  )

extreme_length_genes <- gene_length_analysis %>%
  dplyr::filter(
    length_group %in% c(
      "Shortest 10%",
      "Longest 10%"
    )
  ) %>%
  droplevels()

length_threshold_table <- data.frame(
  group = c(
    "Shortest 10%",
    "Longest 10%"
  ),
  percentile = c(
    extreme_length_fraction,
    1 - extreme_length_fraction
  ),
  threshold_length_kb = c(
    short_length_threshold,
    long_length_threshold
  )
)

write.table(
  gene_length_analysis,
  file = file.path(
    table_dir,
    "RNAPII_gene_length_analysis_table.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  length_threshold_table,
  file = file.path(
    table_dir,
    "RNAPII_gene_length_extreme_group_thresholds.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 8. Descriptive gene-length statistics
# ------------------------------------------------------------------------------

length_summary <- gene_length_analysis %>%
  dplyr::group_by(category) %>%
  dplyr::summarise(
    n_genes = dplyr::n(),
    median_length_kb = median(
      internal_gene_body_length_kb,
      na.rm = TRUE
    ),
    mean_length_kb = mean(
      internal_gene_body_length_kb,
      na.rm = TRUE
    ),
    q1_length_kb = quantile(
      internal_gene_body_length_kb,
      probs = 0.25,
      na.rm = TRUE,
      names = FALSE
    ),
    q3_length_kb = quantile(
      internal_gene_body_length_kb,
      probs = 0.75,
      na.rm = TRUE,
      names = FALSE
    ),
    iqr_length_kb = IQR(
      internal_gene_body_length_kb,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

write.table(
  length_summary,
  file = file.path(
    table_dir,
    "RNAPII_gene_length_summary_by_DESeq2_category.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 9. Shapiro-Wilk diagnostics for gene length within DESeq2 categories
# ------------------------------------------------------------------------------

# Normality is assessed separately within UP, DOWN and NC because the subsequent
# statistical comparison is performed across these three independent groups.
# The test is applied to internal gene-body length, which is the response
# variable used in the Kruskal-Wallis and pairwise Wilcoxon analyses.
set.seed(shapiro_random_seed)

gene_length_shapiro_by_category <- gene_length_analysis %>%
  dplyr::group_by(category) %>%
  dplyr::group_modify(
    ~ run_shapiro_test(
      .x$internal_gene_body_length_kb
    )
  ) %>%
  dplyr::ungroup()

write.table(
  gene_length_shapiro_by_category,
  file = file.path(
    table_dir,
    "RNAPII_gene_length_Shapiro_Wilk_by_DESeq2_category.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 10. Statistical comparison of gene-length distributions
# ------------------------------------------------------------------------------

kruskal_result <- kruskal.test(
  internal_gene_body_length_kb ~ category,
  data = gene_length_analysis
)

kruskal_table <- data.frame(
  test = "Kruskal-Wallis",
  variable = "internal_gene_body_length_kb",
  statistic = unname(kruskal_result$statistic),
  degrees_of_freedom = unname(kruskal_result$parameter),
  p_value = kruskal_result$p.value
)

pairwise_comparisons <- list(
  c("UP", "NC"),
  c("UP", "DOWN"),
  c("DOWN", "NC")
)

pairwise_wilcoxon <- lapply(
  pairwise_comparisons,
  function(comparison) {
    group_1 <- comparison[[1]]
    group_2 <- comparison[[2]]

    values_1 <- gene_length_analysis %>%
      dplyr::filter(category == group_1) %>%
      dplyr::pull(internal_gene_body_length_kb)

    values_2 <- gene_length_analysis %>%
      dplyr::filter(category == group_2) %>%
      dplyr::pull(internal_gene_body_length_kb)

    test_result <- suppressWarnings(
      wilcox.test(
        values_1,
        values_2,
        exact = FALSE
      )
    )

    median_1 <- median(values_1, na.rm = TRUE)
    median_2 <- median(values_2, na.rm = TRUE)

    data.frame(
      group_1 = group_1,
      group_2 = group_2,
      n_group_1 = length(values_1),
      n_group_2 = length(values_2),
      median_group_1_kb = median_1,
      median_group_2_kb = median_2,
      median_difference_kb = median_1 - median_2,
      median_ratio_group1_over_group2 = median_1 / median_2,
      wilcoxon_W = unname(test_result$statistic),
      p_value = test_result$p.value
    )
  }
) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(
    p_adjusted_BH = p.adjust(
      p_value,
      method = "BH"
    ),
    significance_label = dplyr::case_when(
      p_adjusted_BH < 0.001 ~ "***",
      p_adjusted_BH < 0.01 ~ "**",
      p_adjusted_BH < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

write.table(
  kruskal_table,
  file = file.path(
    table_dir,
    "RNAPII_gene_length_Kruskal_Wallis.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  pairwise_wilcoxon,
  file = file.path(
    table_dir,
    "RNAPII_gene_length_pairwise_Wilcoxon_BH.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 11. RNAPII log2FC statistics for the shortest and longest genes
# ------------------------------------------------------------------------------

extreme_length_summary <- extreme_length_genes %>%
  dplyr::group_by(length_group) %>%
  dplyr::summarise(
    n_genes = dplyr::n(),
    median_log2FC = median(
      log2FoldChange,
      na.rm = TRUE
    ),
    mean_log2FC = mean(
      log2FoldChange,
      na.rm = TRUE
    ),
    q1_log2FC = quantile(
      log2FoldChange,
      probs = 0.25,
      na.rm = TRUE,
      names = FALSE
    ),
    q3_log2FC = quantile(
      log2FoldChange,
      probs = 0.75,
      na.rm = TRUE,
      names = FALSE
    ),
    .groups = "drop"
  )

# Normality is assessed separately within the shortest and longest groups.
# Here the tested response variable is log2FoldChange, because gene length is
# used only to define the two groups compared by the Wilcoxon rank-sum test.
set.seed(shapiro_random_seed)

log2fc_shapiro_by_extreme_length_group <- extreme_length_genes %>%
  dplyr::group_by(length_group) %>%
  dplyr::group_modify(
    ~ run_shapiro_test(
      .x$log2FoldChange
    )
  ) %>%
  dplyr::ungroup()

write.table(
  log2fc_shapiro_by_extreme_length_group,
  file = file.path(
    table_dir,
    "RNAPII_log2FC_Shapiro_Wilk_shortest_vs_longest_10_percent.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

extreme_length_wilcoxon <- wilcox.test(
  log2FoldChange ~ length_group,
  data = extreme_length_genes,
  exact = FALSE
)

shortest_median_log2FC <- extreme_length_summary %>%
  dplyr::filter(length_group == "Shortest 10%") %>%
  dplyr::pull(median_log2FC)

longest_median_log2FC <- extreme_length_summary %>%
  dplyr::filter(length_group == "Longest 10%") %>%
  dplyr::pull(median_log2FC)

extreme_length_wilcoxon_table <- data.frame(
  comparison = "Shortest 10% vs Longest 10%",
  shortest_median_log2FC = shortest_median_log2FC,
  longest_median_log2FC = longest_median_log2FC,
  median_difference = shortest_median_log2FC - longest_median_log2FC,
  wilcoxon_W = unname(extreme_length_wilcoxon$statistic),
  p_value = extreme_length_wilcoxon$p.value
)

write.table(
  extreme_length_summary,
  file = file.path(
    table_dir,
    "RNAPII_log2FC_summary_shortest_vs_longest_10_percent.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  extreme_length_wilcoxon_table,
  file = file.path(
    table_dir,
    "RNAPII_log2FC_shortest_vs_longest_10_percent_Wilcoxon.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 12. Violin plot: gene length by RNAPII DESeq2 category
# ------------------------------------------------------------------------------

category_axis_labels <- length_summary %>%
  dplyr::mutate(
    axis_label = paste0(
      as.character(category),
      "\n(n = ",
      n_genes,
      ")"
    )
  ) %>%
  {
    stats::setNames(
      .$axis_label,
      .$category
    )
  }

sig_annotations_category <- pairwise_wilcoxon %>%
  dplyr::filter(
    group_1 == "UP",
    group_2 %in% c("DOWN", "NC"),
    significance_label != "ns"
  ) %>%
  dplyr::mutate(
    xstart = match(group_1, category_plot_levels),
    xend = match(group_2, category_plot_levels),
    label = significance_label
  ) %>%
  dplyr::arrange(xend)

if (nrow(sig_annotations_category) > 0L) {
  y_max_length <- max(
    gene_length_analysis$internal_gene_body_length_kb,
    na.rm = TRUE
  )

  sig_annotations_category <- sig_annotations_category %>%
    dplyr::mutate(
      y = y_max_length * c(1.60, 2.3)[seq_len(n())],
      y_text = y * 1.08
    )
} else {
  sig_annotations_category <- data.frame(
    xstart = numeric(0),
    xend = numeric(0),
    y = numeric(0),
    y_text = numeric(0),
    label = character(0)
  )
}

p_gene_length_by_category <- ggplot(
  gene_length_analysis,
  aes(
    x = category,
    y = internal_gene_body_length_kb,
    fill = category
  )
) +
  geom_violin(
    width = 0.75,
    trim = TRUE,
    alpha = 0.75,
    colour = NA
  ) +
  geom_boxplot(
    width = 0.14,
    outlier.shape = NA,
    alpha = 0.95,
    colour = "black",
    linewidth = 0.35,
    coef = 0
  ) +
  geom_segment(
    data = sig_annotations_category,
    aes(x = xstart, xend = xend, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  geom_segment(
    data = sig_annotations_category,
    aes(x = xstart, xend = xstart, y = y / 1.05, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  geom_segment(
    data = sig_annotations_category,
    aes(x = xend, xend = xend, y = y / 1.05, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  geom_text(
    data = sig_annotations_category,
    aes(x = (xstart + xend) / 2, y = y_text, label = label),
    inherit.aes = FALSE,
    size = 7,
    fontface = "bold"
  ) +
  scale_y_log10() +
  scale_x_discrete(
    labels = category_axis_labels
  ) +
  scale_fill_manual(
    values = category_colours,
    drop = FALSE
  ) +
  labs(
    x = NULL,
    y = "Internal gene-body length (kb, log scale)",
    fill = NULL
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.text = element_text(
      colour = "black"
    ),
    axis.text.x = element_text(
      size = 18,
      colour = "black"
    ),
    legend.position = "none"
  )

print(p_gene_length_by_category)

save_plot(
  p_gene_length_by_category,
  "RNAPII_gene_length_by_DESeq2_category_violin",
  width = 6,
  height = 5
)

# ------------------------------------------------------------------------------
# 13. Violin plot: RNAPII log2FC in the shortest and longest genes
# ------------------------------------------------------------------------------

extreme_axis_labels <- extreme_length_summary %>%
  dplyr::mutate(
    axis_label = paste0(
      length_group,
      "\n(n = ",
      n_genes,
      ")"
    )
  ) %>%
  {
    stats::setNames(
      .$axis_label,
      .$length_group
    )
  }

short_long_significance_label <- dplyr::case_when(
  extreme_length_wilcoxon_table$p_value < 0.001 ~ "***",
  extreme_length_wilcoxon_table$p_value < 0.01 ~ "**",
  extreme_length_wilcoxon_table$p_value < 0.05 ~ "*",
  TRUE ~ "ns"
)

sig_annotation_short_long <- data.frame(
  xstart = 1,
  xend = 2,
  y = max(extreme_length_genes$log2FoldChange, na.rm = TRUE) + 0.45,
  y_text = max(extreme_length_genes$log2FoldChange, na.rm = TRUE) + 0.75,
  label = short_long_significance_label
)

p_log2FC_short_vs_long <- ggplot(
  extreme_length_genes,
  aes(
    x = length_group,
    y = log2FoldChange,
    fill = length_group
  )
) +
  geom_violin(
    width = 0.75,
    trim = TRUE,
    alpha = 0.75,
    colour = NA
  ) +
  geom_boxplot(
    width = 0.14,
    outlier.shape = NA,
    alpha = 0.95,
    colour = "black",
    linewidth = 0.35,
    coef = 0
  ) +
  geom_segment(
    data = sig_annotation_short_long,
    aes(x = xstart, xend = xend, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  geom_segment(
    data = sig_annotation_short_long,
    aes(x = xstart, xend = xstart, y = y - 0.15, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  geom_segment(
    data = sig_annotation_short_long,
    aes(x = xend, xend = xend, y = y - 0.15, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  geom_text(
    data = sig_annotation_short_long,
    aes(x = (xstart + xend) / 2, y = y_text, label = label),
    inherit.aes = FALSE,
    size = 7,
    fontface = "bold"
  ) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.4
  ) +
  scale_x_discrete(
    labels = extreme_axis_labels
  ) +
  scale_fill_manual(
    values = length_group_colours[
      c("Shortest 10%", "Longest 10%")
    ],
    drop = FALSE
  ) +
  labs(
    x = NULL,
    y = "RNAPII ChIP-seq log2FC (CPT vs DMSO)",
    fill = NULL
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.text = element_text(
      colour = "black"
    ),
    axis.text.x = element_text(
      size = 18,
      colour = "black"
    ),
    legend.position = "none"
  )

print(p_log2FC_short_vs_long)

save_plot(
  p_log2FC_short_vs_long,
  "RNAPII_log2FC_shortest_vs_longest_10_percent_violin",
  width = 6,
  height = 5
)

# ------------------------------------------------------------------------------
# 14. Ranked RNAPII log2FC plot by gene-length group
# ------------------------------------------------------------------------------

ranked_log2FC_data <- gene_length_analysis %>%
  dplyr::filter(
    !is.na(log2FoldChange)
  ) %>%
  dplyr::arrange(log2FoldChange) %>%
  dplyr::mutate(
    log2FC_rank = dplyr::row_number()
  )

p_ranked_log2FC_by_length <- ggplot(
  ranked_log2FC_data,
  aes(
    x = log2FC_rank,
    y = log2FoldChange,
    colour = length_group
  )
) +
  geom_point(
    size = 3,
    alpha = 0.8
  ) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.4
  ) +
  geom_hline(
    yintercept = c(-1, 1),
    linetype = "dashed",
    linewidth = 0.35
  ) +
  scale_colour_manual(
    values = length_group_colours,
    breaks = length_group_levels,
    drop = FALSE
  ) +
  labs(
    x = "Genes ranked by RNAPII ChIP-seq log2FC",
    y = "RNAPII ChIP-seq log2FC (CPT vs DMSO)",
    colour = NULL
  ) +
  theme_classic(base_size = 20) +
  theme(
    axis.text = element_text(
      colour = "black"
    ),
    legend.position = "none"
  )

print(p_ranked_log2FC_by_length)

save_plot(
  p_ranked_log2FC_by_length,
  "RNAPII_log2FC_rank_by_gene_length_group",
  width = 9,
  height = 6
)

ranked_log2FC_legend_source <- p_ranked_log2FC_by_length +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 12)
  )

save_plot_legend(
  plot = ranked_log2FC_legend_source,
  output_base = file.path(
    plot_dir,
    "RNAPII_log2FC_rank_gene_length_group_common_legend"
  ),
  width = 3.2,
  height = 2.2
)

# ------------------------------------------------------------------------------
# 15. RNAPII category composition across gene-length deciles
# ------------------------------------------------------------------------------

length_decile_composition <- gene_length_analysis %>%
  dplyr::mutate(
    length_decile = dplyr::ntile(
      internal_gene_body_length_kb,
      10
    )
  ) %>%
  dplyr::count(
    length_decile,
    category,
    name = "n_genes"
  ) %>%
  tidyr::complete(
    length_decile = 1:10,
    category = factor(
      category_plot_levels,
      levels = category_plot_levels
    ),
    fill = list(n_genes = 0L)
  ) %>%
  dplyr::group_by(length_decile) %>%
  dplyr::mutate(
    total_genes_in_decile = sum(n_genes),
    percentage = 100 * n_genes / total_genes_in_decile
  ) %>%
  dplyr::ungroup()

write.table(
  length_decile_composition,
  file = file.path(
    table_dir,
    "RNAPII_category_composition_across_gene_length_deciles.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

p_category_by_length_decile <- ggplot(
  length_decile_composition,
  aes(
    x = length_decile,
    y = percentage,
    colour = category,
    group = category
  )
) +
  geom_line(
    linewidth = 1
  ) +
  geom_point(
    size = 2.5
  ) +
  scale_x_continuous(
    breaks = 1:10
  ) +
  scale_colour_manual(
    values = category_colours,
    breaks = category_plot_levels,
    drop = FALSE
  ) +
  labs(
    x = "Gene-length decile (shortest to longest)",
    y = "Genes in each RNAPII category (%)",
    colour = NULL
  ) +
  theme_classic(base_size = 18) +
  theme(
    axis.text = element_text(
      colour = "black"
    ),
    legend.position = "none"
  )

print(p_category_by_length_decile)

save_plot(
  p_category_by_length_decile,
  "RNAPII_category_composition_across_gene_length_deciles",
  width = 7,
  height = 5
)

category_decile_legend_source <- p_category_by_length_decile +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 12)
  )

save_plot_legend(
  plot = category_decile_legend_source,
  output_base = file.path(
    plot_dir,
    "RNAPII_category_composition_common_legend"
  ),
  width = 2.4,
  height = 2.0
)

# ==============================================================================
# OPTIONAL FINAL ANALYSIS
# Gene-length sensitivity across ranked RNAPII-gain genes
#
# This complete final section can be removed without affecting any preceding
# table, statistical test or plot.
# ==============================================================================

# ------------------------------------------------------------------------------
# 16. Select significant RNAPII-gain genes and rank them by log2FC
# ------------------------------------------------------------------------------

rnapii_gain_genes <- gene_length_analysis %>%
  dplyr::filter(
    !is.na(log2FoldChange),
    !is.na(padj),
    log2FoldChange > rnapii_gain_log2fc_threshold,
    padj <= rnapii_gain_padj_threshold
  ) %>%
  dplyr::arrange(
    dplyr::desc(log2FoldChange)
  )

if (nrow(rnapii_gain_genes) == 0L) {
  stop(
    "No genes met the RNAPII-gain thresholds used for the sensitivity analysis.",
    call. = FALSE
  )
}

# ------------------------------------------------------------------------------
# 17. Calculate gene length for nested ranked RNAPII-gain subsets
# ------------------------------------------------------------------------------

rnapii_gain_length_sensitivity <- lapply(
  seq_along(rnapii_gain_retained_fractions),
  function(i) {
    retained_fraction <- rnapii_gain_retained_fractions[[i]]

    n_retained <- max(
      1L,
      ceiling(nrow(rnapii_gain_genes) * retained_fraction)
    )

    retained_genes <- rnapii_gain_genes %>%
      dplyr::slice_head(n = n_retained)

    data.frame(
      retained_group = rnapii_gain_retained_labels[[i]],
      retained_fraction = retained_fraction,
      n_genes = nrow(retained_genes),
      mean_internal_gene_body_length_bp = mean(
        retained_genes$internal_gene_body_length_bp,
        na.rm = TRUE
      ),
      median_internal_gene_body_length_bp = median(
        retained_genes$internal_gene_body_length_bp,
        na.rm = TRUE
      )
    )
  }
) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(
    retained_group = factor(
      retained_group,
      levels = rnapii_gain_retained_labels
    )
  )

write.table(
  rnapii_gain_length_sensitivity,
  file = file.path(
    table_dir,
    "RNAPII_gain_gene_length_sensitivity_by_log2FC_rank.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

rnapii_gain_length_sensitivity_long <- rnapii_gain_length_sensitivity %>%
  tidyr::pivot_longer(
    cols = c(
      mean_internal_gene_body_length_bp,
      median_internal_gene_body_length_bp
    ),
    names_to = "statistic",
    values_to = "internal_gene_body_length_bp"
  ) %>%
  dplyr::mutate(
    statistic = dplyr::recode(
      statistic,
      mean_internal_gene_body_length_bp = "Mean",
      median_internal_gene_body_length_bp = "Median"
    ),
    statistic = factor(
      statistic,
      levels = c("Mean", "Median")
    )
  )

# ------------------------------------------------------------------------------
# 18. Plot gene length across ranked RNAPII-gain subsets
# ------------------------------------------------------------------------------

p_rnapii_gain_length_sensitivity <- ggplot(
  rnapii_gain_length_sensitivity_long,
  aes(
    x = retained_group,
    y = internal_gene_body_length_bp,
    group = statistic,
    linetype = statistic
  )
) +
  geom_line(
    linewidth = 0.9,
    colour = sensitivity_plot_colour
  ) +
  geom_point(
    size = 2.4,
    colour = sensitivity_plot_colour
  ) +
  scale_linetype_manual(
    values = c(
      "Mean" = "solid",
      "Median" = "dashed"
    )
  ) +
  scale_y_continuous(
    labels = scales::label_number(
      big.mark = ","
    )
  ) +
  labs(
    x = "Genes retained by highest RNAPII log2FC",
    y = "Internal gene-body length (bp)",
    linetype = "Statistic"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(
      colour = "black"
    ),
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    ),
    legend.position = "right"
  )

print(p_rnapii_gain_length_sensitivity)

save_plot(
  p_rnapii_gain_length_sensitivity,
  "RNAPII_gain_gene_length_sensitivity_by_log2FC_rank",
  width = 7.5,
  height = 5
)
