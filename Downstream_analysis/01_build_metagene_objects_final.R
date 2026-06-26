################################################################################
# Build RNA Pol II metagene RDS objects
################################################################################
#
# This script generates the R objects used by the downstream metagene profile and
# heatmap scripts. The final analysis uses a single gene annotation containing
# protein-coding genes annotated as "gene" features.
#
# For each profile type, the script stores:
#   result$profile         Mean signal per bin and sample
#   result$matrices        Named list of genes x bins matrices, one per sample
#   result$gene_metadata   Gene annotation matching the matrix rows
#   result$genes           GRanges object used to build the matrices
#   result$bin_metadata    Bin coordinates and labels for matrix columns
#   result$sample_metadata Sample table used in the analysis
#   result$settings        Parameters and input files used to generate the object
#
# Main outputs:
#   Metageneplot/R_objects/metagene_TSS_protein_coding_min3kb.rds
#   Metageneplot/R_objects/metagene_TSS_TES_protein_coding_min3kb.rds
#   Metageneplot/object_metadata/*_gene_metadata.tsv
#
################################################################################

# ==============================================================================
# 1. Packages
# ==============================================================================

suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(EnrichedHeatmap)
  library(dplyr)
})

# ==============================================================================
# 2. Input and output configuration
# ==============================================================================

outdir_rds <- "../Script_Limpios/Metageneplot/R_objects"
outdir_tables <- "../Script_Limpios/Metageneplot/object_metadata"

dir.create(outdir_rds, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir_tables, recursive = TRUE, showWarnings = FALSE)

# Final annotation used for metagene profiles and heatmaps.
# The annotation contains protein-coding genes located on standard chromosomes,
# with defined strand and a minimum length of 3 kb.
annotation_id <- "protein_coding_min3kb"
annotation_label <- "Protein-coding genes ≥3 kb"
gtf_file <- "../Script_Limpios/gtf_processing/genes_protein_coding_min3kb.gtf"


# Spike-in-normalized BigWig files used to calculate signal matrices.
sample_metadata <- data.frame(
  sample = c("DMSO_R1", "CPT_R1", "DMSO_R2", "CPT_R2"),
  condition = c("DMSO", "CPT", "DMSO", "CPT"),
  replicate = c("Rep1", "Rep1", "Rep2", "Rep2"),
  sample_label = c("DMSO Rep1", "CPT Rep1", "DMSO Rep2", "CPT Rep2"),
  bigwig_file = c(
    "Data_final_CPT/bigwig_spikein_deseq2/R1_RNAPII_AV_DM_S1.coverage.spikein_deseq2.bw",
    "Data_final_CPT/bigwig_spikein_deseq2/R1_RNAPII_AV_CPT_S4.coverage.spikein_deseq2.bw",
    "Data_final_CPT/bigwig_spikein_deseq2/R2_RNAPII_AV_DM_S7.coverage.spikein_deseq2.bw",
    "Data_final_CPT/bigwig_spikein_deseq2/R2_RNAPII_AV_CPT_S10.coverage.spikein_deseq2.bw"
  ),
  stringsAsFactors = FALSE
)
# Metagene parameters.
tss_upstream <- 3000
tss_downstream <- 3000
tss_bin_size <- 50
gene_body_bins <- 100
mean_mode <- "w0"
background <- 0

# ==============================================================================
# 3. Helper functions
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

harmonize_seqlevels_to_signal <- function(target, signal, sample_name = NULL) {
  common_chr <- intersect(seqlevels(target), seqlevels(signal))

  if (length(common_chr) == 0) {
    seqlevelsStyle(target) <- seqlevelsStyle(signal)[1]
    common_chr <- intersect(seqlevels(target), seqlevels(signal))
  }

  if (length(common_chr) == 0) {
    stop(
      "No common chromosomes between annotation and bigWig",
      if (!is.null(sample_name)) paste0(" for sample: ", sample_name) else "",
      ". Check chromosome names, for example '1' vs 'chr1'.",
      call. = FALSE
    )
  }

  keepSeqlevels(target, common_chr, pruning.mode = "coarse")
}

get_metadata_column <- function(gr, candidate_columns) {
  meta <- as.data.frame(mcols(gr))
  selected <- candidate_columns[candidate_columns %in% colnames(meta)]

  if (length(selected) == 0) {
    return(rep(NA_character_, length(gr)))
  }

  as.character(meta[[selected[1]]])
}

make_safe_row_ids <- function(gene_id) {
  # Ensembl gene IDs are used as matrix row identifiers to simplify integration
  # with count matrices, DESeq2 results and other gene-level analyses.
  row_ids <- ifelse(
    !is.na(gene_id) & gene_id != "",
    gene_id,
    paste0("gene_", seq_along(gene_id))
  )

  if (anyDuplicated(row_ids)) {
    warning(
      "Duplicated gene_id values were found. make.unique() will be applied ",
      "only to matrix row names; the original gene_id is kept in gene_metadata."
    )
  }

  make.unique(row_ids, sep = "__dup")
}

make_gene_metadata <- function(genes, profile_type) {
  gene_id <- get_metadata_column(genes, c("gene_id", "geneID", "gene", "ID"))
  gene_name <- get_metadata_column(genes, c("gene_name", "gene_symbol", "symbol", "Name"))
  gene_biotype <- get_metadata_column(genes, c("gene_biotype", "gene_type", "biotype"))

  gene_id_versionless <- sub("\\.[0-9]+$", "", gene_id)

  tss_position <- ifelse(
    as.character(strand(genes)) == "-",
    end(genes),
    start(genes)
  )

  data.frame(
    gene_index = seq_along(genes),
    row_id = make_safe_row_ids(gene_id),
    gene_id = gene_id,
    gene_id_versionless = gene_id_versionless,
    gene_name = gene_name,
    gene_biotype = gene_biotype,
    seqnames = as.character(seqnames(genes)),
    start = start(genes),
    end = end(genes),
    strand = as.character(strand(genes)),
    gene_length_bp = width(genes),
    tss_position = tss_position,
    profile_type = profile_type,
    stringsAsFactors = FALSE
  )
}

make_tss_bin_metadata <- function(n_bins, upstream, downstream, bin_size) {
  upstream_bins <- upstream / bin_size
  downstream_bins <- downstream / bin_size

  if (upstream_bins != floor(upstream_bins) || downstream_bins != floor(downstream_bins)) {
    stop("upstream and downstream must be exact multiples of bin_size.", call. = FALSE)
  }

  upstream_bins <- as.integer(upstream_bins)
  downstream_bins <- as.integer(downstream_bins)
  target_bins <- n_bins - upstream_bins - downstream_bins

  if (target_bins < 0) {
    stop(
      "The number of matrix columns is smaller than upstream + downstream bins. ",
      "Check normalizeToMatrix() output.",
      call. = FALSE
    )
  }

  upstream_x <- seq(
    from = -upstream + bin_size / 2,
    to = -bin_size / 2,
    by = bin_size
  )

  downstream_x <- seq(
    from = bin_size / 2,
    to = downstream - bin_size / 2,
    by = bin_size
  )

  target_x <- rep(0, target_bins)
  x_values <- c(upstream_x, target_x, downstream_x)

  if (length(x_values) != n_bins) {
    stop(
      "TSS bin metadata does not match the number of matrix columns. ",
      "Expected ", n_bins, " bins but generated ", length(x_values), ".",
      call. = FALSE
    )
  }

  region_type <- c(
    rep("upstream", upstream_bins),
    rep("target_TSS", target_bins),
    rep("downstream", downstream_bins)
  )

  data.frame(
    bin_index = seq_len(n_bins),
    x = x_values,
    region_type = region_type,
    bin_start_relative_bp = ifelse(region_type == "target_TSS", 0, x_values - bin_size / 2),
    bin_end_relative_bp = ifelse(region_type == "target_TSS", 1, x_values + bin_size / 2),
    bin_label = paste0("bin_", seq_len(n_bins)),
    stringsAsFactors = FALSE
  )
}

make_gene_body_bin_metadata <- function(n_bins) {
  data.frame(
    bin_index = seq_len(n_bins),
    x = seq(0, 1, length.out = n_bins),
    scaled_position = seq(0, 1, length.out = n_bins),
    bin_label = paste0("bin_", seq_len(n_bins)),
    stringsAsFactors = FALSE
  )
}

prepare_reference_genes <- function(genes, reference_bigwig, profile_type) {
  # The first BigWig is used to define the chromosome set and gene order stored
  # in the object. The remaining samples are checked against this same order.
  signal <- rtracklayer::import(reference_bigwig, format = "BigWig")

  genes_ref <- harmonize_seqlevels_to_signal(
    target = genes,
    signal = signal,
    sample_name = basename(reference_bigwig)
  )

  gene_metadata <- make_gene_metadata(
    genes = genes_ref,
    profile_type = profile_type
  )

  list(
    genes = genes_ref,
    gene_metadata = gene_metadata
  )
}

save_object_and_metadata <- function(result, output_file) {
  saveRDS(result, file = output_file)

  metadata_file <- file.path(
    outdir_tables,
    paste0(tools::file_path_sans_ext(basename(output_file)), "_gene_metadata.tsv")
  )

  write.table(
    result$gene_metadata,
    file = metadata_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  message("  Saved RDS: ", output_file)
  message("  Saved gene metadata table: ", metadata_file)
}

# ==============================================================================
# 4. Metagene object builders
# ==============================================================================

compute_tss_metagene_object <- function(
    genes,
    sample_metadata,
    annotation_id,
    annotation_label,
    upstream = 3000,
    downstream = 3000,
    bin_size = 50,
    mean_mode = "w0",
    background = 0
) {
  message("Building TSS object")

  reference <- prepare_reference_genes(
    genes = genes,
    reference_bigwig = sample_metadata$bigwig_file[1],
    profile_type = "TSS_window"
  )

  genes_ref <- reference$genes
  gene_metadata <- reference$gene_metadata
  row_ids <- gene_metadata$row_id

  # Convert gene ranges into strand-aware 1 bp TSS regions.
  tss_ref <- promoters(genes_ref, upstream = 0, downstream = 1)

  all_profiles <- list()
  all_matrices <- list()
  bin_metadata <- NULL

  for (i in seq_len(nrow(sample_metadata))) {
    sample_name <- sample_metadata$sample[i]
    bigwig_file <- sample_metadata$bigwig_file[i]

    message("  Sample: ", sample_name)

    signal <- rtracklayer::import(bigwig_file, format = "BigWig")

    target <- harmonize_seqlevels_to_signal(
      target = tss_ref,
      signal = signal,
      sample_name = sample_name
    )

    if (length(target) != length(row_ids)) {
      stop(
        "The number of TSS regions does not match the number of gene identifiers ",
        "for sample ", sample_name, ". Check chromosome filtering.",
        call. = FALSE
      )
    }

    mat <- EnrichedHeatmap::normalizeToMatrix(
      signal = signal,
      target = target,
      value_column = "score",
      extend = c(upstream, downstream),
      w = bin_size,
      mean_mode = mean_mode,
      background = background
    )

    mat <- as.matrix(mat)

    if (nrow(mat) != length(row_ids)) {
      stop(
        "The number of matrix rows does not match the number of gene identifiers ",
        "for sample ", sample_name,
        call. = FALSE
      )
    }

    rownames(mat) <- row_ids

    if (is.null(bin_metadata)) {
      bin_metadata <- make_tss_bin_metadata(
        n_bins = ncol(mat),
        upstream = upstream,
        downstream = downstream,
        bin_size = bin_size
      )
    }

    if (ncol(mat) != nrow(bin_metadata)) {
      stop(
        "The number of matrix columns does not match bin_metadata for sample ",
        sample_name,
        call. = FALSE
      )
    }

    colnames(mat) <- bin_metadata$bin_label

    all_profiles[[sample_name]] <- data.frame(
      sample = sample_name,
      x = bin_metadata$x,
      signal = colMeans(mat, na.rm = TRUE),
      stringsAsFactors = FALSE
    )

    all_matrices[[sample_name]] <- mat
  }

  list(
    profile = dplyr::bind_rows(all_profiles),
    matrices = all_matrices,
    gene_metadata = gene_metadata,
    genes = genes_ref,
    bin_metadata = bin_metadata,
    sample_metadata = sample_metadata,
    settings = list(
      object_version = "metagene_object_v3",
      profile_type = "TSS_window",
      annotation_id = annotation_id,
      annotation_label = annotation_label,
      upstream = upstream,
      downstream = downstream,
      bin_size = bin_size,
      mean_mode = mean_mode,
      background = background,
      n_genes = length(genes_ref),
      sample_names = sample_metadata$sample,
      bigwig_files = sample_metadata$bigwig_file
    )
  )
}

compute_gene_body_metagene_object <- function(
    genes,
    sample_metadata,
    annotation_id,
    annotation_label,
    n_bins = 100,
    mean_mode = "w0",
    background = 0
) {
  message("Building TSS-TES object")

  reference <- prepare_reference_genes(
    genes = genes,
    reference_bigwig = sample_metadata$bigwig_file[1],
    profile_type = "scaled_TSS_to_TES"
  )

  genes_ref <- reference$genes
  gene_metadata <- reference$gene_metadata
  row_ids <- gene_metadata$row_id

  all_profiles <- list()
  all_matrices <- list()
  bin_metadata <- NULL

  for (i in seq_len(nrow(sample_metadata))) {
    sample_name <- sample_metadata$sample[i]
    bigwig_file <- sample_metadata$bigwig_file[i]

    message("  Sample: ", sample_name)

    signal <- rtracklayer::import(bigwig_file, format = "BigWig")

    target <- harmonize_seqlevels_to_signal(
      target = genes_ref,
      signal = signal,
      sample_name = sample_name
    )

    if (length(target) != length(row_ids)) {
      stop(
        "The number of gene regions does not match the number of gene identifiers ",
        "for sample ", sample_name, ". Check chromosome filtering.",
        call. = FALSE
      )
    }

    mat <- EnrichedHeatmap::normalizeToMatrix(
      signal = signal,
      target = target,
      value_column = "score",
      extend = c(0, 0),
      k = n_bins,
      mean_mode = mean_mode,
      background = background
    )

    mat <- as.matrix(mat)

    if (nrow(mat) != length(row_ids)) {
      stop(
        "The number of matrix rows does not match the number of gene identifiers ",
        "for sample ", sample_name,
        call. = FALSE
      )
    }

    rownames(mat) <- row_ids

    if (is.null(bin_metadata)) {
      bin_metadata <- make_gene_body_bin_metadata(n_bins = ncol(mat))
    }

    if (ncol(mat) != nrow(bin_metadata)) {
      stop(
        "The number of matrix columns does not match bin_metadata for sample ",
        sample_name,
        call. = FALSE
      )
    }

    colnames(mat) <- bin_metadata$bin_label

    all_profiles[[sample_name]] <- data.frame(
      sample = sample_name,
      x = bin_metadata$x,
      signal = colMeans(mat, na.rm = TRUE),
      stringsAsFactors = FALSE
    )

    all_matrices[[sample_name]] <- mat
  }

  list(
    profile = dplyr::bind_rows(all_profiles),
    matrices = all_matrices,
    gene_metadata = gene_metadata,
    genes = genes_ref,
    bin_metadata = bin_metadata,
    sample_metadata = sample_metadata,
    settings = list(
      object_version = "metagene_object_v3",
      profile_type = "scaled_TSS_to_TES",
      annotation_id = annotation_id,
      annotation_label = annotation_label,
      n_bins = n_bins,
      mean_mode = mean_mode,
      background = background,
      n_genes = length(genes_ref),
      sample_names = sample_metadata$sample,
      bigwig_files = sample_metadata$bigwig_file
    )
  )
}

# ==============================================================================
# 5. Run object generation
# ==============================================================================

check_input_files(gtf_file, label = "GTF")
check_input_files(sample_metadata$bigwig_file, label = "bigWig")

message("============================================================================")
message("Loading final gene annotation")
message("File: ", gtf_file)

genes <- rtracklayer::import(gtf_file)
message("Genes loaded: ", length(genes))

tss_object <- compute_tss_metagene_object(
  genes = genes,
  sample_metadata = sample_metadata,
  annotation_id = annotation_id,
  annotation_label = annotation_label,
  upstream = tss_upstream,
  downstream = tss_downstream,
  bin_size = tss_bin_size,
  mean_mode = mean_mode,
  background = background
)

save_object_and_metadata(
  result = tss_object,
  output_file = file.path(outdir_rds, paste0("metagene_TSS_", annotation_id, ".rds"))
)

gene_body_object <- compute_gene_body_metagene_object(
  genes = genes,
  sample_metadata = sample_metadata,
  annotation_id = annotation_id,
  annotation_label = annotation_label,
  n_bins = gene_body_bins,
  mean_mode = mean_mode,
  background = background
)

save_object_and_metadata(
  result = gene_body_object,
  output_file = file.path(outdir_rds, paste0("metagene_TSS_TES_", annotation_id, ".rds"))
)

message("============================================================================")
message("Done. Metagene RDS objects saved in: ", outdir_rds)
message("Gene metadata tables saved in: ", outdir_tables)
