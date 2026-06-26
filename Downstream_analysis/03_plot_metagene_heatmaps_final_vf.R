################################################################################
# Plot RNAPII metagene heatmaps and assess replicate concordance
################################################################################
#
# This script generates the RNAPII ChIP-seq heatmaps used in the TFM figures
# from the metagene RDS objects created by 01_build_metagene_objects.R.
#
# Signal matrices are not recalculated. The script uses:
#   result$matrices        Gene-by-bin signal matrices for each sample
#   result$sample_metadata Sample conditions and replicate information
#   result$bin_metadata    Genomic or scaled positions represented by each bin
#
# For each metagene profile, biological replicate concordance is assessed
# separately for DMSO and CPT samples using Spearman correlation. Each gene is
# summarized by its mean normalized RNAPII signal across all bins of the
# corresponding profile before calculating the correlation between replicates.
#
# Required input objects:
#   Metageneplot/R_objects/
#     metagene_TSS_protein_coding_min3kb.rds
#     metagene_TSS_TES_protein_coding_min3kb.rds
#
# Output structure:
#   Metageneplot/Heatmaps_plots/
#     TSS/
#     TSS_TES/
#     RNAPII_metagene_replicate_Spearman_correlations.tsv
#
# Each heatmap is exported as:
#   1) A PDF with vector text, titles and legends.
#   2) A high-resolution PNG at 600 dpi.
#
# The correlation table contains the number of genes, Spearman's rho and
# p-value for each condition and metagene profile.
#
################################################################################

# ==============================================================================
# 1. Packages
# ==============================================================================

suppressPackageStartupMessages({
  library(EnrichedHeatmap)
  library(ComplexHeatmap)
  library(circlize)
  library(matrixStats)
  library(grid)
})

# ==============================================================================
# 2. Input and output configuration
# ==============================================================================

rds_dir <- "../Script_Limpios/Metageneplot/R_objects"
outdir <- "../Script_Limpios/Metageneplot/Heatmaps_plots"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Final RDS objects generated from protein-coding genes at least 3 kb long.
rds_info <- data.frame(
  file = c(
    "metagene_TSS_protein_coding_min3kb.rds",
    "metagene_TSS_TES_protein_coding_min3kb.rds"
  ),
  profile_folder = c(
    "TSS",
    "TSS_TES"
  ),
  profile_label = c(
    "RNAPII occupancy around the TSS",
    "RNAPII occupancy across gene bodies"
  ),
  stringsAsFactors = FALSE
)

# Output dimensions and resolution.
# The heatmap body is rasterized at high quality to limit PDF file size,
# whereas titles, annotations and legends remain vector-based.
plot_width <- 5
plot_height <- 7
plot_dpi <- 1680
pdf_raster_quality <- 4

# ==============================================================================
# 3. Sample labels and colour palettes
# ==============================================================================

# Sample labels avoid the R1/R2 notation used for paired-end sequencing reads.
sample_label_map <- c(
  "DMSO_R1" = "DMSO Rep1",
  "CPT_R1"  = "CPT Rep1",
  "DMSO_R2" = "DMSO Rep2",
  "CPT_R2"  = "CPT Rep2"
)

col_fun_control <- circlize::colorRamp2(
  c(0, 10, 25, 50),
  c("white", "#D9F0D3", "#74C476", "#00441B")
)

col_fun_cpt <- circlize::colorRamp2(
  c(0, 10, 25, 50),
  c("white", "#E7D4E8", "#AF8DC3", "#762A83")
)

# Colours used for the average signal profile above each heatmap.
profile_line_colours <- c(
  "Control" = "#007A3D",
  "CPT" = "#6A00A8"
)

# ==============================================================================
# 4. Helper functions
# ==============================================================================

check_input_files <- function(files, label) {
  missing_files <- files[!file.exists(files)]
  
  if (length(missing_files) > 0) {
    stop(
      "Missing ", label, " file(s):\n",
      paste(missing_files, collapse = "\n"),
      call. = FALSE
    )
  }
}

get_sample_label <- function(sample_name) {
  sample_label <- sample_label_map[[sample_name]]
  
  if (is.null(sample_label)) {
    sample_label <- sample_name
  }
  
  sample_label
}

get_sample_condition <- function(sample_name, sample_metadata) {
  if (
    !is.null(sample_metadata) &&
    all(c("sample", "condition") %in% colnames(sample_metadata))
  ) {
    matched_condition <- sample_metadata$condition[
      sample_metadata$sample == sample_name
    ]
    
    if (
      length(matched_condition) == 1 &&
      !is.na(matched_condition)
    ) {
      return(as.character(matched_condition))
    }
  }
  
  # Fallback used only for older objects without complete sample metadata.
  if (grepl("^CPT", sample_name)) {
    return("CPT")
  }
  
  "Control"
}

get_condition_colours <- function(condition) {
  if (condition == "CPT") {
    return(
      list(
        heatmap = col_fun_cpt,
        profile = profile_line_colours[["CPT"]]
      )
    )
  }
  
  list(
    heatmap = col_fun_control,
    profile = profile_line_colours[["Control"]]
  )
}

get_row_order <- function(
    mat,
    bin_metadata = NULL,
    row_order_method = "max",
    early_fraction = 0.20
) {
  if (row_order_method == "max") {
    # TSS heatmaps are ordered by the maximum signal within the TSS window.
    row_score <- matrixStats::rowMaxs(
      mat,
      na.rm = TRUE
    )
    
  } else if (row_order_method == "body_mean") {
    # Gene-body heatmaps are ordered by mean signal after excluding the
    # promoter-proximal fraction of the scaled gene body.
    if (
      is.null(bin_metadata) ||
      !"scaled_position" %in% colnames(bin_metadata)
    ) {
      stop(
        "row_order_method = 'body_mean' requires ",
        "bin_metadata$scaled_position.",
        call. = FALSE
      )
    }
    
    body_columns <- which(
      bin_metadata$scaled_position >= early_fraction
    )
    
    if (length(body_columns) == 0) {
      stop(
        "No bins were selected for row ordering.",
        call. = FALSE
      )
    }
    
    row_score <- rowMeans(
      mat[, body_columns, drop = FALSE],
      na.rm = TRUE
    )
    
  } else {
    stop(
      "Unknown row_order_method: ",
      row_order_method,
      call. = FALSE
    )
  }
  
  order(row_score, decreasing = TRUE)
}

build_metagene_heatmap <- function(
    mat,
    title,
    col_fun,
    profile_line_colour,
    bin_metadata = NULL,
    row_order_method = "max",
    early_fraction = 0.20
) {
  mat <- as.matrix(mat)
  
  row_order <- get_row_order(
    mat = mat,
    bin_metadata = bin_metadata,
    row_order_method = row_order_method,
    early_fraction = early_fraction
  )
  
  mat <- mat[row_order, , drop = FALSE]
  
  # EnrichedHeatmap does not allow column names.
  # Bin information remains stored separately in bin_metadata.
  dimnames(mat) <- list(
    rownames(mat),
    NULL
  )
  
  EnrichedHeatmap::EnrichedHeatmap(
    mat,
    name = "RNAPII signal",
    col = col_fun,
    column_title = title,
    column_title_gp = grid::gpar(
      fontsize = 18,
      fontface = "bold",
      lineheight = 1.1
    ),
    top_annotation = ComplexHeatmap::HeatmapAnnotation(
      enriched = EnrichedHeatmap::anno_enriched(
        gp = grid::gpar(
          col = profile_line_colour,
          lwd = 1.4
        ),
        axis_param = list(
          labels = FALSE,
          gp = grid::gpar(fontsize = 0),
          side = "left"
        )
      ),
      annotation_height = grid::unit(1.4, "cm")
    ),
    show_row_names = FALSE,
    use_raster = TRUE,
    raster_quality = pdf_raster_quality,
    heatmap_legend_param = list(
      title_gp = grid::gpar(
        fontsize = 11,
        fontface = "bold"
      ),
      labels_gp = grid::gpar(fontsize = 9)
    )
  )
}

save_heatmap_files <- function(
    heatmap_object,
    output_base,
    width = plot_width,
    height = plot_height,
    dpi = plot_dpi
) {
  pdf_file <- paste0(output_base, ".pdf")
  png_file <- paste0(output_base, ".png")
  
  grDevices::pdf(
    file = pdf_file,
    width = width,
    height = height,
    useDingbats = FALSE,
    bg = "white"
  )
  
  ComplexHeatmap::draw(
    heatmap_object,
    show_heatmap_legend = FALSE,
    show_annotation_legend = FALSE
  )
  
  grDevices::dev.off()
  
  grDevices::png(
    filename = png_file,
    width = width,
    height = height,
    units = "in",
    res = dpi,
    bg = "white"
  )
  
  ComplexHeatmap::draw(
    heatmap_object,
    show_heatmap_legend = FALSE,
    show_annotation_legend = FALSE
  )
  
  grDevices::dev.off()
  
  message("  Saved vector heatmap without legend: ", pdf_file)
  message("  Saved high-resolution heatmap without legend: ", png_file)
}

save_common_heatmap_legend <- function(
    col_fun,
    output_base,
    title = "RNAPII signal",
    width = 2.4,
    height = 3,
    dpi = plot_dpi
) {
  pdf_file <- paste0(output_base, ".pdf")
  png_file <- paste0(output_base, ".png")
  
  legend_object <- ComplexHeatmap::Legend(
    title = title,
    col_fun = col_fun,
    at = c(0, 25, 50),
    labels = c("0", "25", "50"),
    title_gp = grid::gpar(
      fontsize = 18,
      fontface = "bold"
    ),
    labels_gp = grid::gpar(
      fontsize = 16
    ),
    legend_height = grid::unit(2.3, "cm")
  )
  
  grDevices::pdf(
    file = pdf_file,
    width = width,
    height = height,
    useDingbats = FALSE,
    bg = "white"
  )
  grid::grid.newpage()
  ComplexHeatmap::draw(legend_object)
  grDevices::dev.off()
  
  grDevices::png(
    filename = png_file,
    width = width,
    height = height,
    units = "in",
    res = dpi,
    bg = "white"
  )
  grid::grid.newpage()
  ComplexHeatmap::draw(legend_object)
  grDevices::dev.off()
  
  message("  Saved common heatmap legend PDF: ", pdf_file)
  message("  Saved common heatmap legend PNG: ", png_file)
}

save_combined_heatmap_legend <- function(
    output_base,
    width = 4.2,
    height = 3.0,
    dpi = plot_dpi
) {
  pdf_file <- paste0(output_base, ".pdf")
  png_file <- paste0(output_base, ".png")
  
  legend_dmso <- ComplexHeatmap::Legend(
    title = "DMSO",
    col_fun = col_fun_control,
    at = c(0, 25, 50),
    labels = c("0", "25", "50"),
    title_gp = grid::gpar(
      fontsize = 18,
      fontface = "bold"
    ),
    labels_gp = grid::gpar(
      fontsize = 16
    ),
    legend_height = grid::unit(2.3, "cm")
  )
  
  legend_cpt <- ComplexHeatmap::Legend(
    title = "CPT",
    col_fun = col_fun_cpt,
    at = c(0, 25, 50),
    labels = c("0", "25", "50"),
    title_gp = grid::gpar(
      fontsize = 18,
      fontface = "bold"
    ),
    labels_gp = grid::gpar(
      fontsize = 16
    ),
    legend_height = grid::unit(2.3, "cm")
  )
  
  legend_pack <- ComplexHeatmap::packLegend(
    legend_dmso,
    legend_cpt,
    direction = "horizontal",
    gap = grid::unit(0.7, "cm")
  )
  
  grDevices::pdf(
    file = pdf_file,
    width = width,
    height = height,
    useDingbats = FALSE,
    bg = "white"
  )
  grid::grid.newpage()
  ComplexHeatmap::draw(legend_pack)
  grDevices::dev.off()
  
  grDevices::png(
    filename = png_file,
    width = width,
    height = height,
    units = "in",
    res = dpi,
    bg = "white"
  )
  grid::grid.newpage()
  ComplexHeatmap::draw(legend_pack)
  grDevices::dev.off()
  
  message("  Saved combined heatmap legend PDF: ", pdf_file)
  message("  Saved combined heatmap legend PNG: ", png_file)
}

calculate_spearman <- function(mat_rep1, mat_rep2) {
  mat_rep1 <- as.matrix(mat_rep1)
  mat_rep2 <- as.matrix(mat_rep2)
  
  # Match genes between replicates.
  common_genes <- intersect(
    rownames(mat_rep1),
    rownames(mat_rep2)
  )
  
  if (length(common_genes) == 0) {
    stop(
      "No common genes were found between replicate matrices.",
      call. = FALSE
    )
  }
  
  mat_rep1 <- mat_rep1[common_genes, , drop = FALSE]
  mat_rep2 <- mat_rep2[common_genes, , drop = FALSE]
  
  # Summarize each gene by its mean signal across all bins.
  signal_rep1 <- rowMeans(
    mat_rep1,
    na.rm = TRUE
  )
  
  signal_rep2 <- rowMeans(
    mat_rep2,
    na.rm = TRUE
  )
  
  valid_genes <- is.finite(signal_rep1) &
    is.finite(signal_rep2)
  
  signal_rep1 <- signal_rep1[valid_genes]
  signal_rep2 <- signal_rep2[valid_genes]
  
  correlation_test <- stats::cor.test(
    signal_rep1,
    signal_rep2,
    method = "spearman",
    exact = FALSE
  )
  
  list(
    n_genes = length(signal_rep1),
    rho = unname(correlation_test$estimate),
    p_value = correlation_test$p.value
  )
}

# ==============================================================================
# 5. Calculate replicate correlations and generate heatmaps
# ==============================================================================

rds_paths <- file.path(rds_dir, rds_info$file)
check_input_files(rds_paths, label = "metagene RDS")

message("============================================================================")
message("Plotting final RNAPII metagene heatmaps")

# Store replicate correlation results for all metagene profiles.
replicate_correlation_results <- list()

for (i in seq_len(nrow(rds_info))) {
  rds_path <- file.path(
    rds_dir,
    rds_info$file[i]
  )
  
  message("Processing: ", rds_path)
  
  metagene_object <- readRDS(rds_path)
  
  if (is.null(metagene_object$matrices)) {
    stop(
      "The RDS object does not contain result$matrices: ",
      rds_path,
      call. = FALSE
    )
  }
  # --------------------------------------------------------------------------
  # Calculate Spearman correlations between biological replicates
  # --------------------------------------------------------------------------
  
  required_samples <- c(
    "DMSO_R1",
    "DMSO_R2",
    "CPT_R1",
    "CPT_R2"
  )
  
  missing_samples <- setdiff(
    required_samples,
    names(metagene_object$matrices)
  )
  
  if (length(missing_samples) > 0) {
    stop(
      "Missing sample matrices in ",
      rds_info$profile_folder[i],
      ": ",
      paste(missing_samples, collapse = ", "),
      call. = FALSE
    )
  }
  
  dmso_correlation <- calculate_spearman(
    mat_rep1 = metagene_object$matrices[["DMSO_R1"]],
    mat_rep2 = metagene_object$matrices[["DMSO_R2"]]
  )
  
  cpt_correlation <- calculate_spearman(
    mat_rep1 = metagene_object$matrices[["CPT_R1"]],
    mat_rep2 = metagene_object$matrices[["CPT_R2"]]
  )
  
  replicate_correlation_results[[rds_info$profile_folder[i]]] <- data.frame(
    profile = rep(
      rds_info$profile_folder[i],
      2
    ),
    condition = c(
      "DMSO",
      "CPT"
    ),
    replicate_1 = c(
      "DMSO Rep1",
      "CPT Rep1"
    ),
    replicate_2 = c(
      "DMSO Rep2",
      "CPT Rep2"
    ),
    n_genes = c(
      dmso_correlation$n_genes,
      cpt_correlation$n_genes
    ),
    spearman_rho = c(
      dmso_correlation$rho,
      cpt_correlation$rho
    ),
    p_value = c(
      dmso_correlation$p_value,
      cpt_correlation$p_value
    ),
    stringsAsFactors = FALSE
  )
  
  message("  Replicate Spearman correlations:")
  
  print(
    replicate_correlation_results[[
      rds_info$profile_folder[i]
    ]]
  )
  
  object_outdir <- file.path(
    outdir,
    rds_info$profile_folder[i]
  )
  
  dir.create(
    object_outdir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  sample_names <- names(metagene_object$matrices)
  
  for (sample_name in sample_names) {
    message("  Sample: ", sample_name)
    
    mat <- metagene_object$matrices[[sample_name]]
    sample_label <- get_sample_label(sample_name)
    message("    Genes in heatmap: ", nrow(mat))
    condition <- get_sample_condition(
      sample_name = sample_name,
      sample_metadata = metagene_object$sample_metadata
    )
    
    colours <- get_condition_colours(condition)
    
    plot_title <- sample_label
    
    safe_sample_label <- gsub(
      " ",
      "_",
      sample_label
    )
    
    row_order_method <- if (
      rds_info$profile_folder[i] == "TSS_TES"
    ) {
      "body_mean"
    } else {
      "max"
    }
    
    heatmap_object <- build_metagene_heatmap(
      mat = mat,
      title = plot_title,
      col_fun = colours$heatmap,
      profile_line_colour = colours$profile,
      bin_metadata = metagene_object$bin_metadata,
      row_order_method = row_order_method,
      early_fraction = 0.20
    )
    
    output_base <- file.path(
      object_outdir,
      paste0(
        "heatmap_",
        rds_info$profile_folder[i],
        "_protein_coding_min3kb_",
        safe_sample_label
      )
    )
    
    save_heatmap_files(
      heatmap_object = heatmap_object,
      output_base = output_base
    )
  }
}
# ==============================================================================
# 6. Export replicate correlation results
# ==============================================================================

replicate_correlation_results <- do.call(
  rbind,
  replicate_correlation_results
)

rownames(replicate_correlation_results) <- NULL

correlation_output_file <- file.path(
  outdir,
  "RNAPII_metagene_replicate_Spearman_correlations.tsv"
)

write.table(
  replicate_correlation_results,
  file = correlation_output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ==============================================================================
# 7. Export common heatmap legends
# ==============================================================================

save_common_heatmap_legend(
  col_fun = col_fun_control,
  output_base = file.path(
    outdir,
    "common_heatmap_legend_DMSO"
  ),
  title = "RNAPII signal"
)

save_common_heatmap_legend(
  col_fun = col_fun_cpt,
  output_base = file.path(
    outdir,
    "common_heatmap_legend_CPT"
  ),
  title = "RNAPII signal"
)

save_combined_heatmap_legend(
  output_base = file.path(
    outdir,
    "common_heatmap_legend_DMSO_CPT_combined"
  )
)

message("============================================================================")
message("Replicate correlation results saved in: ", correlation_output_file)

message("============================================================================")
message("Done. Heatmaps saved in: ", outdir)

message("============================================================================")
message("Done. Heatmaps legends saved in: ", outdir)