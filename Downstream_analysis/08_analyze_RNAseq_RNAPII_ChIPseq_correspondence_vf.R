# ==============================================================================
# Analyse RNA-seq and RNAPII ChIP-seq correspondence
# ==============================================================================
# This script compares differential RNA abundance with differential RNAPII
# occupancy within the common active-gene universe defined in the preceding
# filtering script.
#
# Inputs:
#   - RNA-seq DESeq2 results for CPT versus DMSO
#   - RNAPII ChIP-seq DESeq2 results for CPT versus DMSO
#   - common active-gene universe in DMSO
#
# Genes with padj = NA in either dataset are retained in the original DESeq2
# outputs but excluded here because they cannot be assigned to the three final
# comparison categories:
#   - Up_in_CPT
#   - Down_in_CPT
#   - Not_changed
#
# Final outputs:
#   - one integrated correspondence table
#   - differential-category barplot
#   - RNA-seq/RNAPII ChIP-seq transition heatmap
#   - RNA-seq composition by RNAPII ChIP-seq category
#   - RNA-seq versus RNAPII ChIP-seq log2FC scatterplot
# ==============================================================================

library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(tibble)
library(grid)
library(cowplot)

# ------------------------------------------------------------------------------
# 1. Input and output paths
# ------------------------------------------------------------------------------

rna_results_file <- file.path(
  "../Script_Limpios",
  "DESeq2_RNAseq",
  "tables",
  "RNAseq_DESeq2_CPT_vs_DMSO_all_genes.tsv"
)

chip_results_file <- file.path(
  "../Script_Limpios",
  "DESeq2_ChIPseq",
  "RNAPII_internal_gene_body",
  "tables",
  "RNAPII_ChIPseq_DESeq2_CPT_vs_DMSO_all_genes.tsv"
)

common_universe_file <- file.path(
  "../Script_Limpios",
  "RNAseq_ChIPseq_correspondence",
  "common_gene_universe",
  "tables",
  "common_active_genes_DMSO_RNAseq_RNAPII_ChIPseq.tsv"
)

output_dir <- file.path(
  "../Script_Limpios",
  "RNAseq_ChIPseq_correspondence",
  "correspondence_analysis"
)

table_dir <- file.path(output_dir, "tables")
plot_dir <- file.path(output_dir, "plots")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

input_files <- c(
  rna_results_file,
  chip_results_file,
  common_universe_file
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
# 2. Helper functions and plot settings
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
    device = cairo_pdf,
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

comparison_categories <- c(
  "Up_in_CPT",
  "Down_in_CPT",
  "Not_changed"
)

differential_colours <- c(
  "Up_in_CPT" = "#C44E8C",
  "Down_in_CPT" = "#2C7FB8",
  "Not_changed" = "#D9D9D9"
)

correspondence_levels <- c(
  "NC in both",
  "Concordant UP",
  "Concordant DOWN",
  "Opposite switch",
  "Changed in one"
)

correspondence_colours <- c(
  "NC in both" = "#8C8C8C",
  "Concordant UP" = differential_colours[["Up_in_CPT"]],
  "Concordant DOWN" = differential_colours[["Down_in_CPT"]],
  "Opposite switch" = "#6A3D9A",
  "Changed in one" = "#FDB863"
)

# ------------------------------------------------------------------------------
# 3. Import and validate input tables
# ------------------------------------------------------------------------------

rna_results <- read.delim(
  rna_results_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

chip_results <- read.delim(
  chip_results_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

common_universe <- read.delim(
  common_universe_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_result_columns <- c(
  "gene_id",
  "log2FoldChange",
  "padj",
  "significance"
)

if (!all(required_result_columns %in% colnames(rna_results))) {
  stop(
    "The RNA-seq DESeq2 table does not contain all required columns.",
    call. = FALSE
  )
}

if (!all(required_result_columns %in% colnames(chip_results))) {
  stop(
    "The RNAPII ChIP-seq DESeq2 table does not contain all required columns.",
    call. = FALSE
  )
}

if (!"gene_id" %in% colnames(common_universe)) {
  stop(
    "The common active-gene universe does not contain a gene_id column.",
    call. = FALSE
  )
}

rna_results <- rna_results %>%
  dplyr::mutate(
    gene_id = clean_gene_id(gene_id),
    significance = dplyr::if_else(
      significance == "NC",
      "Not_changed",
      significance
    )
  )

chip_results <- chip_results %>%
  dplyr::mutate(
    gene_id = clean_gene_id(gene_id),
    significance = dplyr::if_else(
      significance == "NC",
      "Not_changed",
      significance
    )
  )

common_gene_ids <- common_universe$gene_id %>%
  clean_gene_id() %>%
  unique()

if (anyDuplicated(rna_results$gene_id)) {
  stop("Duplicated gene IDs were found in the RNA-seq DESeq2 table.",
       call. = FALSE)
}

if (anyDuplicated(chip_results$gene_id)) {
  stop("Duplicated gene IDs were found in the RNAPII ChIP-seq DESeq2 table.",
       call. = FALSE)
}

valid_input_categories <- c(comparison_categories, "padj_NA")

if (!all(rna_results$significance %in% valid_input_categories)) {
  stop("Unexpected significance categories were found in RNA-seq.",
       call. = FALSE)
}

if (!all(chip_results$significance %in% valid_input_categories)) {
  stop("Unexpected significance categories were found in RNAPII ChIP-seq.",
       call. = FALSE)
}

# ------------------------------------------------------------------------------
# 4. Build the integrated correspondence table
# ------------------------------------------------------------------------------

rna_table <- rna_results %>%
  dplyr::select(
    gene_id,
    RNAseq_log2FC = log2FoldChange,
    RNAseq_padj = padj,
    RNAseq_significance = significance
  )

chip_table <- chip_results %>%
  dplyr::select(
    gene_id,
    RNAPII_ChIPseq_log2FC = log2FoldChange,
    RNAPII_ChIPseq_padj = padj,
    RNAPII_ChIPseq_significance = significance
  )

integrated_all <- dplyr::inner_join(
  rna_table,
  chip_table,
  by = "gene_id"
) %>%
  dplyr::filter(gene_id %in% common_gene_ids)

# Only genes assigned to the three final categories in both datasets are used in
# correspondence plots and in the final integrated correspondence table.
correspondence_table <- integrated_all %>%
  dplyr::filter(
    RNAseq_significance %in% comparison_categories,
    RNAPII_ChIPseq_significance %in% comparison_categories
  ) %>%
  dplyr::mutate(
    correspondence_class = dplyr::case_when(
      RNAseq_significance == "Not_changed" &
        RNAPII_ChIPseq_significance == "Not_changed" ~ "NC in both",

      RNAseq_significance == "Up_in_CPT" &
        RNAPII_ChIPseq_significance == "Up_in_CPT" ~ "Concordant UP",

      RNAseq_significance == "Down_in_CPT" &
        RNAPII_ChIPseq_significance == "Down_in_CPT" ~ "Concordant DOWN",

      RNAseq_significance %in% c("Up_in_CPT", "Down_in_CPT") &
        RNAPII_ChIPseq_significance %in% c("Up_in_CPT", "Down_in_CPT") &
        RNAseq_significance != RNAPII_ChIPseq_significance ~ "Opposite switch",

      TRUE ~ "Changed in one"
    ),
    correspondence_class = factor(
      correspondence_class,
      levels = correspondence_levels
    )
  ) %>%
  dplyr::arrange(gene_id)

if (nrow(correspondence_table) == 0L) {
  stop(
    "No genes remained after applying the common-universe and category filters.",
    call. = FALSE
  )
}

write.table(
  correspondence_table,
  file = file.path(
    table_dir,
    "RNAseq_RNAPII_ChIPseq_correspondence_common_active_genes.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  correspondence_table,
  file = file.path(
    table_dir,
    "RNAseq_RNAPII_ChIPseq_correspondence_common_active_genes.rds"
  )
)

# ------------------------------------------------------------------------------
# 5. Differential-category barplot
# ------------------------------------------------------------------------------

rna_category_summary <- correspondence_table %>%
  dplyr::count(significance = RNAseq_significance) %>%
  dplyr::mutate(dataset = "RNA-seq")

chip_category_summary <- correspondence_table %>%
  dplyr::count(significance = RNAPII_ChIPseq_significance) %>%
  dplyr::mutate(dataset = "RNAPII ChIP-seq")

category_summary <- dplyr::bind_rows(
  rna_category_summary,
  chip_category_summary
) %>%
  tidyr::complete(
    dataset = c("RNA-seq", "RNAPII ChIP-seq"),
    significance = comparison_categories,
    fill = list(n = 0L)
  ) %>%
  dplyr::mutate(
    dataset = factor(
      dataset,
      levels = c("RNA-seq", "RNAPII ChIP-seq")
    ),
    significance = factor(
      significance,
      levels = comparison_categories
    )
  )

category_barplot <- ggplot(
  category_summary,
  aes(
    x = dataset,
    y = n,
    fill = significance
  )
) +
  geom_col(width = 0.4) +
  geom_text(
    aes(label = ifelse(n > 0, n, "")),
    position = position_stack(vjust = 0.5),
    size = 6,
    fontface = "bold",
    colour = "black"
  ) +
  scale_x_discrete(
    labels = c(
      "RNA-seq" = "RNA-seq\n (HeLa)",
      "RNAPII ChIP-seq" = "RNAPII ChIP-seq\n (RPE-1)"
    )
  ) +
  scale_fill_manual(
    values = differential_colours,
    breaks = comparison_categories,
    labels = c(
      "Up in CPT",
      "Down in CPT",
      "NC"
    ),
    drop = FALSE
  ) +
  labs(
    x = NULL,
    y = "Number of genes",
    fill = NULL
  ) +
  theme_classic(base_size = 22) +
  theme(
    axis.text = element_text(colour = "black"),
    legend.position = "none"
  )

save_plot(
  category_barplot,
  "RNAseq_RNAPII_ChIPseq_differential_categories",
  width = 6,
  height = 5
)

# ------------------------------------------------------------------------------
# 5b. Differential-category barplot legend
# ------------------------------------------------------------------------------

category_barplot_legend_source <- ggplot(
  category_summary,
  aes(
    x = dataset,
    y = n,
    fill = significance
  )
) +
  geom_col(width = 0.4) +
  scale_fill_manual(
    values = differential_colours,
    breaks = comparison_categories,
    labels = c(
      "Up in CPT",
      "Down in CPT",
      "NC"
    ),
    drop = FALSE
  ) +
  guides(
    fill = guide_legend(
      ncol = 1,
      byrow = TRUE,
      override.aes = list(
        alpha = 1
      )
    )
  ) +
  labs(
    fill = NULL
  ) +
  theme_void(base_size = 14) +
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    legend.text = element_text(size = 14),
    legend.key.size = grid::unit(0.45, "cm"),
    legend.spacing.y = grid::unit(0.15, "cm")
  )

category_barplot_common_legend <- cowplot::get_legend(
  category_barplot_legend_source
)

category_barplot_common_legend_plot <- cowplot::ggdraw(
  category_barplot_common_legend
)

ggsave(
  filename = file.path(
    plot_dir,
    "RNAseq_RNAPII_ChIPseq_differential_categories_legend.pdf"
  ),
  plot = category_barplot_common_legend_plot,
  width = 2.2,
  height = 1.8,
  device = cairo_pdf,
  bg = "white"
)

ggsave(
  filename = file.path(
    plot_dir,
    "RNAseq_RNAPII_ChIPseq_differential_categories_legend.png"
  ),
  plot = category_barplot_common_legend_plot,
  width = 2.2,
  height = 1.8,
  dpi = safe_png_dpi(2.2, 1.8),
  bg = "white"
)

# ------------------------------------------------------------------------------
# 6. RNA-seq/RNAPII ChIP-seq transition heatmap
# ------------------------------------------------------------------------------

transition_table <- correspondence_table %>%
  dplyr::count(
    RNAseq_significance,
    RNAPII_ChIPseq_significance
  ) %>%
  tidyr::complete(
    RNAseq_significance = comparison_categories,
    RNAPII_ChIPseq_significance = comparison_categories,
    fill = list(n = 0L)
  ) %>%
  dplyr::mutate(
    RNAseq_significance = factor(
      RNAseq_significance,
      levels = comparison_categories
    ),
    RNAPII_ChIPseq_significance = factor(
      RNAPII_ChIPseq_significance,
      levels = comparison_categories
    )
  )

transition_heatmap <- ggplot(
  transition_table,
  aes(
    x = RNAPII_ChIPseq_significance,
    y = RNAseq_significance,
    fill = n
  )
) +
  geom_tile(
    colour = "white",
    linewidth = 1
  ) +
  geom_text(
    aes(label = n),
    size = 7,
    fontface = "bold",
    colour = "black"
  ) +
  scale_x_discrete(
    labels = c(
      "Up_in_CPT" = "RNAPII\n UP",
      "Down_in_CPT" = "RNAPII\n DOWN",
      "Not_changed" = "RNAPII\n NC"
    )
  ) +
  scale_y_discrete(
    labels = c(
      "Up_in_CPT" = "RNA\n UP",
      "Down_in_CPT" = "RNA\n DOWN",
      "Not_changed" = "RNA\n NC"
    )
  ) +
  scale_fill_gradientn(
    colours = c(
      "#FFFFCC",
      "#FFEDA0",
      "#FED976",
      "#FEB24C",
      "#FD8D3C",
      "#F03B20",
      "#BD0026"
    ),
    name = "Gene count"
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 18) +
  theme(
    axis.text.x = element_text(
      angle = 35,
      hjust = 1,
      colour = "black"
    ),
    axis.text.y = element_text(colour = "black"),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 14)
  )

save_plot(
  transition_heatmap,
  "RNAseq_RNAPII_ChIPseq_category_transition_heatmap",
  width = 6,
  height = 5
)

# ------------------------------------------------------------------------------
# 7. RNA-seq composition by RNAPII ChIP-seq category
# ------------------------------------------------------------------------------

composition_table <- correspondence_table %>%
  dplyr::count(
    RNAPII_ChIPseq_significance,
    correspondence_class
  ) %>%
  tidyr::complete(
    RNAPII_ChIPseq_significance = comparison_categories,
    correspondence_class = correspondence_levels,
    fill = list(n = 0L)
  ) %>%
  dplyr::group_by(RNAPII_ChIPseq_significance) %>%
  dplyr::mutate(
    percentage = 100 * n / sum(n),
    percentage_label = dplyr::if_else(
      percentage >= 1,
      paste0(round(percentage, 1), "%"),
      ""
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    RNAPII_ChIPseq_significance = factor(
      RNAPII_ChIPseq_significance,
      levels = comparison_categories
    ),
    correspondence_class = factor(
      correspondence_class,
      levels = correspondence_levels
    )
  )

composition_barplot <- ggplot(
  composition_table,
  aes(
    x = RNAPII_ChIPseq_significance,
    y = percentage,
    fill = correspondence_class
  )
) +
  geom_col(width = 0.7) +
  geom_text(
    aes(label = percentage_label),
    position = position_stack(vjust = 0.5),
    size = 5,
    fontface = "bold",
    colour = "black"
  ) +
  scale_x_discrete(
    labels = c(
      "Up_in_CPT" = "UP",
      "Down_in_CPT" = "DOWN",
      "Not_changed" = "NC"
    )
  ) +
  scale_fill_manual(
    values = correspondence_colours,
    breaks = correspondence_levels,
    drop = FALSE
  ) +
  scale_y_continuous(
    breaks = seq(0, 100, by = 25),
    limits = c(0, 100),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = "RNAPII ChIP-seq category",
    y = "Percentage of genes",
    fill = NULL
  ) +
  theme_classic(base_size = 18) +
  theme(
    axis.text = element_text(colour = "black"),
    legend.position = "none"
  )

save_plot(
  composition_barplot,
  "RNAseq_composition_by_RNAPII_ChIPseq_category",
  width = 7,
  height = 5
)

# ------------------------------------------------------------------------------
# 8. RNA-seq versus RNAPII ChIP-seq log2FC scatterplot
# ------------------------------------------------------------------------------

scatter_data <- correspondence_table %>%
  dplyr::filter(
    !is.na(RNAPII_ChIPseq_log2FC),
    !is.na(RNAseq_log2FC)
  )

if (nrow(scatter_data) < 3L) {
  stop(
    "Fewer than three genes are available for the scatterplot correlation.",
    call. = FALSE
  )
}

spearman_test <- suppressWarnings(
  cor.test(
    scatter_data$RNAPII_ChIPseq_log2FC,
    scatter_data$RNAseq_log2FC,
    method = "spearman",
    exact = FALSE
  )
)

scatter_data_RNAPII_UP <- scatter_data %>%
  dplyr::filter(RNAPII_ChIPseq_significance == "Up_in_CPT")

if (nrow(scatter_data_RNAPII_UP) < 3L) {
  stop(
    "Fewer than three RNAPII-UP genes are available for correlation analysis.",
    call. = FALSE
  )
}

spearman_test_RNAPII_UP <- suppressWarnings(
  cor.test(
    scatter_data_RNAPII_UP$RNAPII_ChIPseq_log2FC,
    scatter_data_RNAPII_UP$RNAseq_log2FC,
    method = "spearman",
    exact = FALSE
  )
)

spearman_summary <- tibble::tibble(
  gene_set = c(
    "All common active genes",
    "RNAPII UP genes"
  ),
  n_genes = c(
    nrow(scatter_data),
    nrow(scatter_data_RNAPII_UP)
  ),
  spearman_rho = c(
    unname(spearman_test$estimate),
    unname(spearman_test_RNAPII_UP$estimate)
  ),
  p_value = c(
    spearman_test$p.value,
    spearman_test_RNAPII_UP$p.value
  )
)

write.table(
  spearman_summary,
  file = file.path(
    table_dir,
    "RNAseq_RNAPII_ChIPseq_Spearman_correlations.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

spearman_rho <- round(
  unname(spearman_test$estimate),
  3
)

spearman_p_label <- if (
  spearman_test$p.value < 2.2e-16
) {
  "P < 2.2e-16"
} else {
  paste0("P = ", signif(spearman_test$p.value, 3))
}

spearman_label <- paste0(
  "Spearman rho = ",
  spearman_rho,
  "\n",
  spearman_p_label
)

scatterplot <- ggplot(
  scatter_data,
  aes(
    x = RNAPII_ChIPseq_log2FC,
    y = RNAseq_log2FC
  )
) +
  geom_point(
    data = scatter_data %>%
      dplyr::filter(correspondence_class == "NC in both"),
    aes(colour = correspondence_class),
    alpha = 0.18,
    size = 0.55
  ) +
  geom_point(
    data = scatter_data %>%
      dplyr::filter(correspondence_class != "NC in both"),
    aes(colour = correspondence_class),
    alpha = 0.75,
    size = 1.2
  ) +
  geom_vline(
    xintercept = c(-1, 0, 1),
    linetype = c("dashed", "solid", "dashed"),
    linewidth = c(0.35, 0.4, 0.35)
  ) +
  geom_hline(
    yintercept = c(-1, 0, 1),
    linetype = c("dashed", "solid", "dashed"),
    linewidth = c(0.35, 0.4, 0.35)
  ) +
  annotate(
    "text",
    x = Inf,
    y = -Inf,
    label = spearman_label,
    hjust = 1.1,
    vjust = -0.8,
    size = 6
  ) +
  scale_colour_manual(
    values = correspondence_colours,
    breaks = correspondence_levels,
    drop = FALSE
  ) +
  labs(
    x = "RNAPII ChIP-seq log2FC (CPT vs DMSO)",
    y = "RNA-seq log2FC (CPT vs DMSO)",
    colour = NULL
  ) +
  theme_classic(base_size = 20) +
  theme(
    axis.text = element_text(colour = "black"),
    legend.position = "none",
    plot.margin = margin(10, 15, 10, 15)
  )

save_plot(
  scatterplot,
  "RNAPII_ChIPseq_RNAseq_log2FC_correspondence_scatter",
  width = 7,
  height = 5
)

# ------------------------------------------------------------------------------
# 9. Common legend for RNA-seq/RNAPII correspondence plots
# ------------------------------------------------------------------------------

correspondence_legend_source <- ggplot(
  composition_table,
  aes(
    x = RNAPII_ChIPseq_significance,
    y = percentage,
    fill = correspondence_class
  )
) +
  geom_col(width = 0.7) +
  scale_fill_manual(
    values = correspondence_colours,
    breaks = correspondence_levels,
    drop = FALSE
  ) +
  guides(
    fill = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        alpha = 1
      )
    )
  ) +
  labs(
    fill = NULL
  ) +
  theme_void(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 12),
    legend.key.size = grid::unit(0.45, "cm"),
    legend.spacing.x = grid::unit(0.35, "cm")
  )

correspondence_common_legend <- cowplot::get_legend(
  correspondence_legend_source
)

correspondence_common_legend_plot <- cowplot::ggdraw(
  correspondence_common_legend
)

ggsave(
  filename = file.path(
    plot_dir,
    "RNAseq_RNAPII_correspondence_common_legend.pdf"
  ),
  plot = correspondence_common_legend_plot,
  width = 9,
  height = 1.0,
  device = cairo_pdf,
  bg = "white"
)

ggsave(
  filename = file.path(
    plot_dir,
    "RNAseq_RNAPII_correspondence_common_legend.png"
  ),
  plot = correspondence_common_legend_plot,
  width = 9,
  height = 1.0,
  dpi = safe_png_dpi(6.7, 1.0),
  bg = "white"
)