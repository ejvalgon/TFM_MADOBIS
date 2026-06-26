# ==============================================================================
# PCA from deepTools multiBamSummary and multiBigwigSummary
# RNA Pol II ChIP-seq - AV Control vs CPT
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load libraries
# ------------------------------------------------------------------------------

library(tidyverse)
library(cowplot)

# ------------------------------------------------------------------------------
# 2. Input files
# ------------------------------------------------------------------------------

multi_bam_file <- "Data_final_CPT/multiBamSummary/RNA_pol_II_CHIPseq_multiBamSummary.txt"

multi_bw_median_file <- "../bigwigs_CORRECTOS/multiBigwigSummary/RNA_pol_II_CHIPseq_Final_BW_multiBigwigSummary.spikein_median.txt"

multi_bw_deseq2_file <- "../bigwigs_CORRECTOS/multiBigwigSummary/RNA_pol_II_CHIPseq_Final_BW_multiBigwigSummary.spikein_deseq2.txt"

deseq2_spikein_counts_file <- "../Script_Limpios/DESeq2_ChIPseq/RNAPII_internal_gene_body/tables/RNAPII_ChIPseq_spikein_normalized_counts.tsv"
  

outdir <- "../Script_Limpios/PCA/plots"

dir.create(outdir, recursive = TRUE)

# ------------------------------------------------------------------------------
# 2B. Plot export settings
# ------------------------------------------------------------------------------

# PNG resolution is calculated automatically to maximize image quality while
# keeping each file below the 100-megapixel limit of the figure-assembly software.
maximum_png_pixels <- 100000000

safe_png_dpi <- function(width, height, maximum_pixels = maximum_png_pixels) {
  floor(
    sqrt(
      maximum_pixels / (width * height)
    )
  )
}

save_png <- function(plot, filename, width, height) {
  ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    dpi = safe_png_dpi(width, height),
    bg = "white"
  )
}

# ------------------------------------------------------------------------------
# 3. Define selected samples and sample metadata
# ------------------------------------------------------------------------------

samples_av <- c(
  "R1_RNAPII_AV_CPT_S4",
  "R1_RNAPII_AV_DM_S1",
  "R2_RNAPII_AV_CPT_S10",
  "R2_RNAPII_AV_DM_S7"
)

metadata_av <- tibble(
  sample = samples_av,
  condition = c("CPT", "Control", "CPT", "Control"),
  replicate = c("R1", "R1", "R2", "R2"),
  sample_id = c("CPT_R1", "DM_R1", "CPT_R2", "DM_R2")
)

# Make condition a factor to control plotting order
metadata_av <- metadata_av %>%
  mutate(
    condition = factor(condition, levels = c("Control", "CPT")),
    replicate = factor(replicate, levels = c("R1", "R2")),
    sample_id = factor(sample_id, levels = c("DM_R1", "DM_R2", "CPT_R1", "CPT_R2"))
  )

# Colors for each sample
chipseq_colors <- c(
  "DM_R1"  = "#007A3D",  # dark green
  "DM_R2"  = "#7AD151",  # light green
  "CPT_R1" = "#6A00A8",  # dark purple
  "CPT_R2" = "#C77CFF"   # light purple
)

# Labels displayed in the PCA legend
chipseq_sample_labels <- c(
  "DM_R1"  = "DMSO Rep1",
  "DM_R2"  = "DMSO Rep2",
  "CPT_R1" = "CPT Rep1",
  "CPT_R2" = "CPT Rep2"
)

# ------------------------------------------------------------------------------
# Metadata for DESeq2 spike-in normalized count matrix
# ------------------------------------------------------------------------------

samples_deseq2 <- c(
  "DMSO_R1",
  "DMSO_R2",
  "CPT_R1",
  "CPT_R2"
)

metadata_deseq2 <- tibble(
  sample = samples_deseq2,
  condition = c("Control", "Control", "CPT", "CPT"),
  replicate = c("R1", "R2", "R1", "R2"),
  sample_id = c("DM_R1", "DM_R2", "CPT_R1", "CPT_R2")
)

metadata_deseq2 <- metadata_deseq2 %>%
  mutate(
    condition = factor(
      condition,
      levels = c("Control", "CPT")
    ),
    replicate = factor(
      replicate,
      levels = c("R1", "R2")
    ),
    sample_id = factor(
      sample_id,
      levels = c("DM_R1", "DM_R2", "CPT_R1", "CPT_R2")
    )
  )

# ------------------------------------------------------------------------------
# 4. Function to read deepTools raw count matrix
# ------------------------------------------------------------------------------

read_deeptools_matrix <- function(file, keep_canonical = TRUE) {
  
  message("Reading file: ", file)
  
  mat <- read.delim(
    file,
    header = TRUE,
    comment.char = "",
    check.names = FALSE,
    na.strings = c("NA", "NaN", "nan")
  )
  
  # Clean deepTools column names
  colnames(mat) <- gsub("'", "", colnames(mat))
  colnames(mat)[1] <- gsub("^#", "", colnames(mat)[1])
  
  # Standardize coordinate column names
  colnames(mat)[1:3] <- c("chr", "start", "end")
  
  message("Initial number of bins: ", nrow(mat))
  
  # Keep only canonical chromosomes: 1-22, X and Y
  if (keep_canonical) {
    canonical_chr <- c(as.character(1:22), "X", "Y")
    
    mat <- mat %>%
      filter(chr %in% canonical_chr)
    
    message("Bins after keeping canonical chromosomes: ", nrow(mat))
  }
  
  coordinates <- mat[, 1:3]
  signal_matrix <- mat[, -c(1:3)]
  
  signal_matrix <- as.matrix(signal_matrix)
  mode(signal_matrix) <- "numeric"
  
  rownames(signal_matrix) <- paste(
    coordinates$chr,
    coordinates$start,
    coordinates$end,
    sep = "_"
  )
  
  return(signal_matrix)
}

# ------------------------------------------------------------------------------
# 4B. Function to read DESeq2 spike-in normalized count matrix
# ------------------------------------------------------------------------------

read_deseq2_normalized_counts <- function(file) {
  
  message("Reading DESeq2 normalized count matrix: ", file)
  
  mat <- read.delim(
    file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE
  )
  
  stopifnot("gene_id" %in% colnames(mat))
  
  signal_matrix <- mat %>%
    tibble::column_to_rownames("gene_id") %>%
    as.matrix()
  
  mode(signal_matrix) <- "numeric"
  
  # Clean possible Ensembl version suffixes
  rownames(signal_matrix) <- stringr::str_remove(
    rownames(signal_matrix),
    "\\.\\d+$"
  )
  
  message("Genes in DESeq2 normalized matrix: ", nrow(signal_matrix))
  
  return(signal_matrix)
}

# ------------------------------------------------------------------------------
# 5. Function to prepare matrix for PCA
# ------------------------------------------------------------------------------

prepare_matrix_for_pca <- function(signal_matrix, samples, metadata) {
  
  # Check that selected samples exist in the matrix
  missing_samples <- setdiff(samples, colnames(signal_matrix))
  
  if (length(missing_samples) > 0) {
    stop(
      "The following samples are missing from the matrix: ",
      paste(missing_samples, collapse = ", ")
    )
  }
  
  # Select only AV samples
  signal_matrix <- signal_matrix[, samples, drop = FALSE]
  
  # Match metadata order to matrix column order
  metadata <- metadata %>%
    filter(sample %in% samples) %>%
    arrange(match(sample, colnames(signal_matrix)))
  
  stopifnot(all(metadata$sample == colnames(signal_matrix)))
  
  message("Bins before filtering: ", nrow(signal_matrix))
  
  # Remove bins with NA/nan values
  signal_matrix <- signal_matrix[complete.cases(signal_matrix), ]
  message("Bins after removing NA/nan bins: ", nrow(signal_matrix))
  
  # Remove bins with zero signal across all selected samples
  signal_matrix <- signal_matrix[rowSums(signal_matrix) > 0, ]
  message("Bins after removing all-zero bins: ", nrow(signal_matrix))
  
  # Remove bins with no variance across samples
  bin_variance <- apply(signal_matrix, 1, var)
  signal_matrix <- signal_matrix[bin_variance > 0, ]
  message("Bins after removing zero-variance bins: ", nrow(signal_matrix))
  
  # Log2 transformation
  signal_matrix_log <- log2(signal_matrix + 1)
  
  return(list(
    matrix = signal_matrix_log,
    metadata = metadata
  ))
}

# ------------------------------------------------------------------------------
# 6. Function to run PCA
# ------------------------------------------------------------------------------

run_pca <- function(signal_matrix_log) {
  
  pca <- prcomp(
    t(signal_matrix_log),
    center = TRUE,
    scale. = TRUE
  )
  
  return(pca)
}

# ------------------------------------------------------------------------------
# 7. Functions to plot PCA and export the common legend
# ------------------------------------------------------------------------------

plot_pca <- function(pca, metadata, title, output_file, colors) {
  
  pca_df <- as.data.frame(pca$x) %>%
    rownames_to_column("sample") %>%
    left_join(metadata, by = "sample")
  
  variance_explained <- summary(pca)$importance[2, ] * 100
  
  p <- ggplot(
    pca_df,
    aes(
      x = PC1,
      y = PC2,
      color = sample_id,
      shape = condition
    )
  ) +
    geom_point(size = 5) +
    scale_color_manual(
      values = colors,
      breaks = names(colors),
      labels = chipseq_sample_labels[names(colors)]
    ) +
    scale_shape_manual(
      values = c(
        "Control" = 16,
        "CPT" = 17
      )
    ) +
    theme_classic(base_size = 24) +
    labs(
      title = title,
      x = paste0("PC1: ", round(variance_explained[1], 1), "% variance"),
      y = paste0("PC2: ", round(variance_explained[2], 1), "% variance"),
      color = "Sample",
      shape = "Condition"
    ) +
    theme(
      legend.position = "none"
    )
  
  save_png(
    plot = p,
    filename = output_file,
    width = 7,
    height = 5
  )
  
  return(p)
}

plot_pca_common_legend <- function(metadata, colors, output_file) {
  
  legend_data <- metadata %>%
    distinct(sample_id, condition) %>%
    mutate(
      PC1 = seq_len(n()),
      PC2 = seq_len(n())
    )
  
  sample_shapes <- c(
    "DM_R1" = 16,
    "DM_R2" = 16,
    "CPT_R1" = 17,
    "CPT_R2" = 17
  )
  
  legend_source <- ggplot(
    legend_data,
    aes(
      x = PC1,
      y = PC2,
      color = sample_id,
      shape = sample_id
    )
  ) +
    geom_point(size = 4) +
    scale_color_manual(
      values = colors,
      breaks = names(colors),
      labels = chipseq_sample_labels[names(colors)],
      guide = "none"
    ) +
    scale_shape_manual(
      values = sample_shapes,
      breaks = names(colors),
      labels = chipseq_sample_labels[names(colors)],
      name = NULL
    ) +
    guides(
      shape = guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(
          color = colors[names(colors)],
          size = 4
        )
      )
    ) +
    theme_void(base_size = 14) +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.text = element_text(size = 14),
      legend.key.size = grid::unit(0.45, "cm"),
      legend.spacing.x = grid::unit(0.35, "cm")
    )
  
  common_legend <- cowplot::get_legend(legend_source)
  common_legend_plot <- cowplot::ggdraw(common_legend)
  
  save_png(
    plot = common_legend_plot,
    filename = output_file,
    width = 6.7,
    height = 0.7
  )
  
  return(common_legend_plot)
}

pca_common_legend <- plot_pca_common_legend(
  metadata = metadata_av,
  colors = chipseq_colors,
  output_file = file.path(outdir, "PCA_common_legend.png")
)

# ------------------------------------------------------------------------------
# 8. PCA from multiBamSummary
# ------------------------------------------------------------------------------

bam_matrix <- read_deeptools_matrix(
  file = multi_bam_file,
  keep_canonical = TRUE
)

bam_prepared <- prepare_matrix_for_pca(
  signal_matrix = bam_matrix,
  samples = samples_av,
  metadata = metadata_av
)

pca_bam <- run_pca(
  signal_matrix_log = bam_prepared$matrix
)

p_bam <- plot_pca(
  pca = pca_bam,
  metadata = bam_prepared$metadata,
  title = "PCA Non-normalized multiBamSummary",
  output_file = file.path(outdir, "PCA_multiBamSummary_Control_vs_CPT.png"),
  colors = chipseq_colors
)

p_bam

# ------------------------------------------------------------------------------
# 9. PCA from multiBigwigSummary - spike-in DESeq2 normalized signal
# ------------------------------------------------------------------------------

bw_deseq2_matrix <- read_deeptools_matrix(
  file = multi_bw_deseq2_file,
  keep_canonical = TRUE
)

# Clean bigWig column names to match BAM sample names
colnames(bw_deseq2_matrix) <- gsub(
  ".coverage.spikein_deseq2$",
  "",
  colnames(bw_deseq2_matrix)
)

bw_deseq2_prepared <- prepare_matrix_for_pca(
  signal_matrix = bw_deseq2_matrix,
  samples = samples_av,
  metadata = metadata_av
)

pca_bw_deseq2 <- run_pca(
  signal_matrix_log = bw_deseq2_prepared$matrix
)

p_bw_deseq2 <- plot_pca(
  pca = pca_bw_deseq2,
  metadata = bw_deseq2_prepared$metadata,
  title = "PCA multiBigwigSummary DESeq2 normalized",
  output_file = file.path(outdir,"PCA_multiBigwigSummary_spikein_deseq2_Control_vs_CPT.png"),
  colors = chipseq_colors
)

p_bw_deseq2

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# 10. PCA from multiBigwigSummary - spike-in median normalized signal
# ------------------------------------------------------------------------------

bw_median_matrix <- read_deeptools_matrix(
  file = multi_bw_median_file,
  keep_canonical = TRUE
)

# Clean bigWig column names to match BAM sample names
colnames(bw_median_matrix) <- gsub(
  "\\.coverage\\.spikein_median$",
  "",
  colnames(bw_median_matrix)
)

bw_median_prepared <- prepare_matrix_for_pca(
  signal_matrix = bw_median_matrix,
  samples = samples_av,
  metadata = metadata_av
)

pca_bw_median <- run_pca(
  signal_matrix_log = bw_median_prepared$matrix
)

p_bw_median <- plot_pca(
  pca = pca_bw_median,
  metadata = bw_median_prepared$metadata,
  title = "PCA multiBigwigSummary median normalized",
  output_file = file.path(outdir, "PCA_multiBigwigSummary_spikein_median_Control_vs_CPT.png"),
  colors = chipseq_colors
)

p_bw_median


# ------------------------------------------------------------------------------
# 11. PCA from DESeq2 spike-in normalized count matrix
# ------------------------------------------------------------------------------

deseq2_spikein_matrix <- read_deseq2_normalized_counts(
  file = deseq2_spikein_counts_file
)

deseq2_spikein_prepared <- prepare_matrix_for_pca(
  signal_matrix = deseq2_spikein_matrix,
  samples = samples_deseq2,
  metadata = metadata_deseq2
)

pca_deseq2_spikein <- run_pca(
  signal_matrix_log = deseq2_spikein_prepared$matrix
)

p_deseq2_spikein <- plot_pca(
  pca = pca_deseq2_spikein,
  metadata = deseq2_spikein_prepared$metadata,
  title = "PCA DESeq2 spike-in normalized counts",
  output_file = file.path(
    outdir,
    "PCA_DESeq2_spikein_normalized_counts_Control_vs_CPT.png"
  ),
  colors = chipseq_colors
)

p_deseq2_spikein
