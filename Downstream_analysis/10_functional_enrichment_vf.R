# ==============================================================================
# FUNCTIONAL ENRICHMENT ANALYSIS OF CONCORDANT UP GENES
# ==============================================================================
# This script analyses genes classified as Up_in_CPT in both RNA-seq and
# RNAPII ChIP-seq.
#
# Target gene set:
#   Concordant UP genes from the final RNA-seq/RNAPII ChIP-seq correspondence
#   table.
#
# Background universe:
#   All genes in the common active-gene universe used for the correspondence
#   analysis.
#
# Main outputs:
#   - GO BP, GO MF, GO CC and KEGG ORA result tables
#   - GO BP result after semantic simplification (cutoff = 0.70)
#   - GO DNA damage-response and DNA-repair branch tables
#   - Descriptive overlap with the GO DNA-repair branch
#   - Main figures: KEGG dotplot, DNA-damage-response gene barplot,
#     DNA-repair gene dumbbell plot and GO BP term-similarity tree
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(GO.db)
library(GOSemSim)
library(cowplot)

# ------------------------------------------------------------------------------
# 2. Input and output paths
# ------------------------------------------------------------------------------

correspondence_file <- file.path(
  "../Script_Limpios",
  "RNAseq_ChIPseq_correspondence",
  "correspondence_analysis",
  "tables",
  "RNAseq_RNAPII_ChIPseq_correspondence_common_active_genes.rds"
)

enrichment_dir <- file.path(
  "../Script_Limpios",
  "RNAseq_ChIPseq_correspondence",
  "functional_enrichment",
  "concordant_UP_ORA"
)

table_dir <- file.path(enrichment_dir, "tables")
plot_dir <- file.path(enrichment_dir, "plots")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(correspondence_file)) {
  stop(
    "The correspondence table was not found:\n",
    correspondence_file,
    call. = FALSE
  )
}
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
  ggplot2::ggsave(
    filename = file.path(plot_dir, paste0(filename, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    device = grDevices::cairo_pdf,
    bg = "white"
  )
  
  ggplot2::ggsave(
    filename = file.path(plot_dir, paste0(filename, ".png")),
    plot = plot,
    width = width,
    height = height,
    dpi = safe_png_dpi(width, height),
    bg = "white"
  )
}

write_tsv_base <- function(x, filename) {
  write.table(
    x,
    file = file.path(table_dir, filename),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

wrap_enrichment_descriptions <- function(enrichment_object, width = 42L) {
  if (!is.null(enrichment_object) && nrow(as.data.frame(enrichment_object)) > 0L) {
    enrichment_object@result$Description <- stringr::str_wrap(
      enrichment_object@result$Description,
      width = width
    )
  }

  enrichment_object
}

rank_enrichment_terms <- function(enrichment_df, n_terms = 15L) {
  enrichment_df %>%
    dplyr::arrange(p.adjust, pvalue, dplyr::desc(Count)) %>%
    dplyr::slice_head(n = min(n_terms, nrow(enrichment_df)))
}

build_GO_gene_table_from_terms <- function(selected_terms) {
  if (is.null(selected_terms) || nrow(selected_terms) == 0L) {
    return(
      tibble::tibble(
        gene_symbol = character(),
        matching_GO_terms = character(),
        matching_GO_IDs = character(),
        n_matching_terms = integer(),
        minimum_p_adjust = numeric()
      )
    )
  }

  selected_terms %>%
    dplyr::select(ID, Description, p.adjust, geneID) %>%
    tidyr::separate_rows(geneID, sep = "/") %>%
    dplyr::rename(gene_symbol = geneID) %>%
    dplyr::filter(!is.na(gene_symbol), gene_symbol != "") %>%
    dplyr::group_by(gene_symbol) %>%
    dplyr::summarise(
      matching_GO_terms = paste(sort(unique(Description)), collapse = "; "),
      matching_GO_IDs = paste(sort(unique(ID)), collapse = "; "),
      n_matching_terms = dplyr::n_distinct(ID),
      minimum_p_adjust = min(p.adjust, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(minimum_p_adjust, gene_symbol)
}

plot_kegg_dotplot <- function(enrichment_object, subtitle, filename) {
  enrichment_df <- as.data.frame(enrichment_object)

  if (nrow(enrichment_df) == 0L) {
    message("No significant KEGG pathways available.")
    return(NULL)
  }

  plot_object <- wrap_enrichment_descriptions(enrichment_object, width = 20L)

  p <- enrichplot::dotplot(
    plot_object,
    showCategory = nrow(enrichment_df),
    x = "GeneRatio"
  ) +
    ggplot2::labs(
      title = "KEGG pathway enrichment",
      subtitle = subtitle,
      x = "Gene ratio",
      y = NULL
    ) +
    ggplot2::theme_classic(base_size = 18) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 14),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10),
      axis.text.y = ggplot2::element_text(size = 18),
      plot.margin = ggplot2::margin(t = 10, r = 15, b = 10, l = 15)
    )

  print(p)
  save_plot(
    p,
    filename,
    width = 8,
    height = max(4.5, nrow(enrichment_df) * 0.9)
  )

  p
}

# ------------------------------------------------------------------------------
# 4. Load and validate the correspondence table
# ------------------------------------------------------------------------------

correspondence_table <- readRDS(correspondence_file)

required_columns <- c(
  "gene_id",
  "RNAseq_log2FC",
  "RNAseq_padj",
  "RNAseq_significance",
  "RNAPII_ChIPseq_log2FC",
  "RNAPII_ChIPseq_padj",
  "RNAPII_ChIPseq_significance",
  "correspondence_class"
)

missing_columns <- setdiff(required_columns, colnames(correspondence_table))

if (length(missing_columns) > 0L) {
  stop(
    "The correspondence table is missing required columns:\n",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

correspondence_table <- correspondence_table %>%
  dplyr::mutate(gene_id = clean_gene_id(gene_id))

if (anyDuplicated(correspondence_table$gene_id)) {
  stop("Duplicated gene IDs were found in the correspondence table.", call. = FALSE)
}

# ------------------------------------------------------------------------------
# 5. Define the target gene set and background universe
# ------------------------------------------------------------------------------

background_genes <- correspondence_table %>%
  dplyr::pull(gene_id) %>%
  unique()

concordant_UP_genes <- correspondence_table %>%
  dplyr::filter(
    RNAseq_significance == "Up_in_CPT",
    RNAPII_ChIPseq_significance == "Up_in_CPT"
  ) %>%
  dplyr::arrange(dplyr::desc(RNAseq_log2FC)) %>%
  dplyr::select(
    gene_id,
    RNAseq_log2FC,
    RNAseq_padj,
    RNAPII_ChIPseq_log2FC,
    RNAPII_ChIPseq_padj,
    RNAseq_significance,
    RNAPII_ChIPseq_significance,
    correspondence_class
  )

target_genes <- unique(concordant_UP_genes$gene_id)

if (length(target_genes) == 0L) {
  stop("No concordant UP genes were found.", call. = FALSE)
}

if (!all(target_genes %in% background_genes)) {
  stop("Some target genes are absent from the background universe.", call. = FALSE)
}

message("Common active-gene universe: ", length(background_genes))
message("Concordant UP target genes: ", length(target_genes))

# ------------------------------------------------------------------------------
# 6. Save analysis input lists
# ------------------------------------------------------------------------------

write_tsv_base(
  concordant_UP_genes,
  "concordant_UP_RNAseq_RNAPII_ChIPseq_genes.tsv"
)

write_tsv_base(
  data.frame(gene_id = target_genes),
  "concordant_UP_gene_list.tsv"
)

write_tsv_base(
  data.frame(gene_id = background_genes),
  "common_active_gene_universe.tsv"
)

# ------------------------------------------------------------------------------
# 7. Map Ensembl IDs to Entrez IDs and gene symbols
# ------------------------------------------------------------------------------

all_ids <- unique(c(background_genes, target_genes))

id_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = all_ids,
  keytype = "ENSEMBL",
  columns = c("ENTREZID", "SYMBOL")
) %>%
  tibble::as_tibble() %>%
  dplyr::distinct()

background_id_map <- id_map %>%
  dplyr::filter(ENSEMBL %in% background_genes)

target_id_map <- id_map %>%
  dplyr::filter(ENSEMBL %in% target_genes)

mapped_background_ensembl <- background_id_map %>%
  dplyr::filter(!is.na(ENSEMBL)) %>%
  dplyr::pull(ENSEMBL) %>%
  unique()

mapped_target_ensembl <- target_id_map %>%
  dplyr::filter(!is.na(ENSEMBL)) %>%
  dplyr::pull(ENSEMBL) %>%
  unique()

background_entrez <- background_id_map %>%
  dplyr::filter(!is.na(ENTREZID)) %>%
  dplyr::pull(ENTREZID) %>%
  unique()

target_entrez <- target_id_map %>%
  dplyr::filter(!is.na(ENTREZID)) %>%
  dplyr::pull(ENTREZID) %>%
  unique()

unmapped_background <- setdiff(background_genes, mapped_background_ensembl)
unmapped_target <- setdiff(target_genes, mapped_target_ensembl)

mapping_summary <- tibble::tibble(
  gene_set = c("Concordant UP", "Background universe"),
  input_ensembl = c(length(target_genes), length(background_genes)),
  mapped_ensembl = c(
    length(mapped_target_ensembl),
    length(mapped_background_ensembl)
  ),
  mapped_entrez = c(length(target_entrez), length(background_entrez)),
  unmapped_ensembl = c(length(unmapped_target), length(unmapped_background))
)

print(mapping_summary)

write_tsv_base(target_id_map, "concordant_UP_gene_ID_mapping.tsv")
write_tsv_base(background_id_map, "background_gene_ID_mapping.tsv")
write_tsv_base(mapping_summary, "gene_ID_mapping_summary.tsv")
write_tsv_base(
  data.frame(gene_id = unmapped_target),
  "concordant_UP_unmapped_gene_IDs.tsv"
)
write_tsv_base(
  data.frame(gene_id = unmapped_background),
  "background_unmapped_gene_IDs.tsv"
)

# ------------------------------------------------------------------------------
# 8. GO over-representation analysis
# ------------------------------------------------------------------------------

run_go_ora <- function(ontology) {
  clusterProfiler::enrichGO(
    gene = target_genes,
    universe = background_genes,
    OrgDb = org.Hs.eg.db,
    keyType = "ENSEMBL",
    ont = ontology,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
}

ego_bp <- run_go_ora("BP")
ego_mf <- run_go_ora("MF")
ego_cc <- run_go_ora("CC")

ego_bp_df <- as.data.frame(ego_bp)
ego_mf_df <- as.data.frame(ego_mf)
ego_cc_df <- as.data.frame(ego_cc)

write_tsv_base(ego_bp_df, "GO_BP_concordant_UP_ORA.tsv")
write_tsv_base(ego_mf_df, "GO_MF_concordant_UP_ORA.tsv")
write_tsv_base(ego_cc_df, "GO_CC_concordant_UP_ORA.tsv")

message("Significant GO BP terms: ", nrow(ego_bp_df))
message("Significant GO MF terms: ", nrow(ego_mf_df))
message("Significant GO CC terms: ", nrow(ego_cc_df))

# ------------------------------------------------------------------------------
# 9. Semantic simplification of GO BP terms
# ------------------------------------------------------------------------------
# GO BP terms with Wang semantic similarity above 0.70 are treated as
# redundant. Within each redundant group, the term with the lowest adjusted
# p-value is retained. The complete unsimplified ORA result is not modified.

ego_bp_simplified_070 <- clusterProfiler::simplify(
  ego_bp,
  cutoff = 0.70,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang"
)

ego_bp_simplified_070_df <- as.data.frame(ego_bp_simplified_070)

write_tsv_base(
  ego_bp_simplified_070_df,
  "GO_BP_concordant_UP_ORA_simplified_cutoff_0.70.tsv"
)

message(
  "GO BP terms retained after semantic simplification (cutoff 0.70): ",
  nrow(ego_bp_simplified_070_df)
)

# ------------------------------------------------------------------------------
# 10. KEGG over-representation analysis
# ------------------------------------------------------------------------------

if (length(target_entrez) > 0L && length(background_entrez) > 0L) {
  ekegg <- clusterProfiler::enrichKEGG(
    gene = target_entrez,
    universe = background_entrez,
    organism = "hsa",
    keyType = "kegg",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )

  ekegg_df <- as.data.frame(ekegg)

  if (nrow(ekegg_df) > 0L) {
    ekegg <- clusterProfiler::setReadable(
      ekegg,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID"
    )
    ekegg_df <- as.data.frame(ekegg)
  }
} else {
  ekegg <- NULL
  ekegg_df <- data.frame()
}

write_tsv_base(ekegg_df, "KEGG_concordant_UP_ORA.tsv")
message("Significant KEGG pathways: ", nrow(ekegg_df))

# ------------------------------------------------------------------------------
# 11. Select the GO BP terms used in the main treeplot
# ------------------------------------------------------------------------------
# The treeplot uses the 15 terms with the lowest adjusted p-values. Raw p-value
# and Count are used only to resolve ties. The complete GO BP result is retained
# in the corresponding output table.

top15_go_bp <- rank_enrichment_terms(
  ego_bp_df,
  n_terms = 15L
)

write_tsv_base(
  top15_go_bp,
  "GO_BP_top15_terms_shown_without_simplification.tsv"
)

# ------------------------------------------------------------------------------
# 12. Extract significant GO BP terms from predefined ontology branches
# ------------------------------------------------------------------------------
# Branch membership is defined from the GO hierarchy rather than by keyword
# matching. These are directed summaries of the significant GO BP ORA result,
# not additional enrichment tests.

DNA_damage_response_root_GO <- "GO:0006974"
DNA_repair_root_GO <- "GO:0006281"

DNA_damage_response_GO_ids <- unique(
  c(
    DNA_damage_response_root_GO,
    GO.db::GOBPOFFSPRING[[DNA_damage_response_root_GO]]
  )
)

DNA_repair_GO_ids <- unique(
  c(
    DNA_repair_root_GO,
    GO.db::GOBPOFFSPRING[[DNA_repair_root_GO]]
  )
)

DNA_damage_response_significant_terms <- ego_bp_df %>%
  dplyr::filter(ID %in% DNA_damage_response_GO_ids) %>%
  dplyr::arrange(p.adjust, pvalue, dplyr::desc(Count))

DNA_repair_significant_terms <- ego_bp_df %>%
  dplyr::filter(ID %in% DNA_repair_GO_ids) %>%
  dplyr::arrange(p.adjust, pvalue, dplyr::desc(Count))

write_tsv_base(
  DNA_damage_response_significant_terms,
  "GO_BP_significant_DNA_damage_response_branch_terms.tsv"
)
write_tsv_base(
  DNA_repair_significant_terms,
  "GO_BP_significant_DNA_repair_branch_terms.tsv"
)

DNA_damage_response_genes <- build_GO_gene_table_from_terms(
  DNA_damage_response_significant_terms
)

DNA_repair_genes_from_significant_terms <- build_GO_gene_table_from_terms(
  DNA_repair_significant_terms
)

write_tsv_base(
  DNA_damage_response_genes,
  "GO_BP_significant_DNA_damage_response_branch_genes.tsv"
)
write_tsv_base(
  DNA_repair_genes_from_significant_terms,
  "GO_BP_significant_DNA_repair_branch_genes.tsv"
)

message(
  "Significant GO BP terms in the DNA damage-response branch: ",
  nrow(DNA_damage_response_significant_terms)
)
message(
  "Significant GO BP terms in the DNA-repair branch: ",
  nrow(DNA_repair_significant_terms)
)
message(
  "Genes represented in significant DNA damage-response GO BP terms: ",
  nrow(DNA_damage_response_genes)
)
message(
  "Genes represented in significant DNA-repair GO BP terms: ",
  nrow(DNA_repair_genes_from_significant_terms)
)

# ------------------------------------------------------------------------------
# 13. Descriptive overlap with the complete GO DNA-repair branch
# ------------------------------------------------------------------------------
# This analysis identifies concordant UP genes annotated to GO:0006281 or any
# descendant Biological Process term. It is an annotation overlap and does not
# test functional enrichment.

target_GO_annotations <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = target_genes,
  keytype = "ENSEMBL",
  columns = c("SYMBOL", "GOALL", "ONTOLOGYALL", "EVIDENCEALL")
) %>%
  tibble::as_tibble() %>%
  dplyr::filter(!is.na(GOALL), ONTOLOGYALL == "BP") %>%
  dplyr::distinct()

DNA_repair_GO_overlap_long <- target_GO_annotations %>%
  dplyr::filter(GOALL %in% DNA_repair_GO_ids) %>%
  dplyr::rename(
    gene_id = ENSEMBL,
    gene_symbol = SYMBOL,
    GO_ID = GOALL,
    GO_evidence = EVIDENCEALL
  ) %>%
  dplyr::mutate(GO_term = AnnotationDbi::Term(GO_ID)) %>%
  dplyr::left_join(
    concordant_UP_genes %>%
      dplyr::select(
        gene_id,
        RNAseq_log2FC,
        RNAseq_padj,
        RNAPII_ChIPseq_log2FC,
        RNAPII_ChIPseq_padj
      ),
    by = "gene_id"
  ) %>%
  dplyr::filter(!is.na(gene_symbol), gene_symbol != "") %>%
  dplyr::distinct(
    gene_id,
    gene_symbol,
    GO_ID,
    GO_term,
    GO_evidence,
    .keep_all = TRUE
  ) %>%
  dplyr::arrange(gene_symbol, GO_term)

DNA_repair_GO_gene_summary <- DNA_repair_GO_overlap_long %>%
  dplyr::group_by(gene_id, gene_symbol) %>%
  dplyr::summarise(
    n_DNA_repair_GO_terms = dplyr::n_distinct(GO_ID),
    DNA_repair_GO_IDs = paste(sort(unique(GO_ID)), collapse = "; "),
    DNA_repair_GO_terms = paste(sort(unique(GO_term)), collapse = "; "),
    GO_evidence_codes = paste(sort(unique(GO_evidence)), collapse = "; "),
    RNAseq_log2FC = dplyr::first(RNAseq_log2FC),
    RNAseq_padj = dplyr::first(RNAseq_padj),
    RNAPII_ChIPseq_log2FC = dplyr::first(RNAPII_ChIPseq_log2FC),
    RNAPII_ChIPseq_padj = dplyr::first(RNAPII_ChIPseq_padj),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(n_DNA_repair_GO_terms), gene_symbol)

DNA_repair_GO_overlap_summary <- tibble::tibble(
  target_gene_set = "Concordant UP genes",
  total_target_genes = length(target_genes),
  DNA_repair_annotated_genes = dplyr::n_distinct(
    DNA_repair_GO_gene_summary$gene_id
  ),
  fraction_DNA_repair_annotated = dplyr::n_distinct(
    DNA_repair_GO_gene_summary$gene_id
  ) / length(target_genes),
  DNA_repair_GO_terms_represented = dplyr::n_distinct(
    DNA_repair_GO_overlap_long$GO_ID
  )
)

write_tsv_base(
  DNA_repair_GO_overlap_long,
  "GO_DNA_repair_concordant_UP_term_gene_overlap.tsv"
)
write_tsv_base(
  DNA_repair_GO_gene_summary,
  "GO_DNA_repair_concordant_UP_gene_summary.tsv"
)
write_tsv_base(
  DNA_repair_GO_overlap_summary,
  "GO_DNA_repair_concordant_UP_overlap_summary.tsv"
)

print(DNA_repair_GO_overlap_summary)

# ------------------------------------------------------------------------------
# 14. Save a concise analysis summary
# ------------------------------------------------------------------------------

enrichment_summary <- tibble::tibble(
  analysis = c(
    "GO BP",
    "GO BP simplified cutoff 0.70",
    "GO MF",
    "GO CC",
    "KEGG",
    "Significant DNA damage-response GO BP branch",
    "Significant DNA-repair GO BP branch"
  ),
  target_input_genes = c(
    length(target_genes),
    length(target_genes),
    length(target_genes),
    length(target_genes),
    length(target_entrez),
    length(target_genes),
    length(target_genes)
  ),
  background_input_genes = c(
    length(background_genes),
    length(background_genes),
    length(background_genes),
    length(background_genes),
    length(background_entrez),
    length(background_genes),
    length(background_genes)
  ),
  retained_or_significant_terms = c(
    nrow(ego_bp_df),
    nrow(ego_bp_simplified_070_df),
    nrow(ego_mf_df),
    nrow(ego_cc_df),
    nrow(ekegg_df),
    nrow(DNA_damage_response_significant_terms),
    nrow(DNA_repair_significant_terms)
  )
)

write_tsv_base(enrichment_summary, "ORA_analysis_summary.tsv")
print(enrichment_summary)

# ==============================================================================
# MAIN FIGURES
# ==============================================================================

# ------------------------------------------------------------------------------
# 15. KEGG pathway-enrichment dotplot
# ------------------------------------------------------------------------------

if (!is.null(ekegg) && nrow(ekegg_df) > 0L) {
  p_kegg_dotplot <- plot_kegg_dotplot(
    ekegg,
    subtitle = "",
    filename = "KEGG_concordant_UP_dotplot"
  )
} else {
  p_kegg_dotplot <- NULL
}

# ------------------------------------------------------------------------------
# 16. Genes represented in significant DNA damage-response GO BP terms
# ------------------------------------------------------------------------------
# The bar length shows how many significant GO terms contain each gene. The
# colour reports the smallest adjusted p-value among those terms. This count
# reflects GO annotation structure and should not be interpreted as gene-level
# biological importance.

if (nrow(DNA_damage_response_genes) > 0L) {
  
  DNA_damage_response_plot_data <- DNA_damage_response_genes %>%
    dplyr::mutate(
      minus_log10_minimum_p_adjust = -log10(minimum_p_adjust)
    ) %>%
    dplyr::arrange(
      dplyr::desc(n_matching_terms),
      dplyr::desc(minus_log10_minimum_p_adjust),
      gene_symbol
    ) %>%
    dplyr::slice_head(
      n = min(15L, nrow(DNA_damage_response_genes))
    ) %>%
    dplyr::mutate(
      gene_symbol = factor(
        gene_symbol,
        levels = rev(gene_symbol)
      )
    )
  
  p_DNA_damage_response_genes <- ggplot2::ggplot(
    DNA_damage_response_plot_data,
    ggplot2::aes(
      x = gene_symbol,
      y = n_matching_terms,
      fill = minus_log10_minimum_p_adjust
    )
  ) +
    ggplot2::geom_col(
      width = 0.72
    ) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(
      breaks = scales::pretty_breaks()
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Number of associated GO terms",
      fill = expression(-log[10]("adjusted P-value"))
    ) +
    ggplot2::theme_classic(
      base_size = 24
    ) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(
        size = 24,
        colour = "black"
      ),
      axis.text.x = ggplot2::element_text(
        size = 20,
        colour = "black"
      ),
      axis.title.x = ggplot2::element_text(
        size = 22
      ),
      legend.position = "none",
      plot.margin = ggplot2::margin(
        t = 10,
        r = 15,
        b = 10,
        l = 15
      )
    )
  
  print(
    p_DNA_damage_response_genes
  )
  
  save_plot(
    p_DNA_damage_response_genes,
    "GO_BP_significant_DNA_damage_response_genes_barplot",
    width = 9,
    height = 7
  )
  
  DNA_damage_response_legend_source <- ggplot2::ggplot(
    DNA_damage_response_plot_data,
    ggplot2::aes(
      x = gene_symbol,
      y = n_matching_terms,
      fill = minus_log10_minimum_p_adjust
    )
  ) +
    ggplot2::geom_col(
      width = 0.72
    ) +
    ggplot2::scale_fill_gradient(
      name = expression(-log[10]("adjusted P-value")),
      low = "#132E43",
      high = "#56B1F7"
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = grid::unit(0.45, "cm"),
        barheight = grid::unit(3.5, "cm")
      )
    ) +
    ggplot2::theme_void(
      base_size = 18
    ) +
    ggplot2::theme(
      legend.position = "right",
      legend.direction = "vertical",
      legend.title = ggplot2::element_text(
        size = 18
      ),
      legend.text = ggplot2::element_text(
        size = 16
      )
    )
  
  DNA_damage_response_legend <- cowplot::get_legend(
    DNA_damage_response_legend_source
  )
  
  p_DNA_damage_response_genes_legend <- cowplot::ggdraw(
    DNA_damage_response_legend
  )
  
  save_plot(
    p_DNA_damage_response_genes_legend,
    "GO_BP_significant_DNA_damage_response_genes_barplot_legend",
    width = 2.8,
    height = 3
  )
  
} else {
  
  p_DNA_damage_response_genes <- NULL
}

# ------------------------------------------------------------------------------
# 17. Concordant UP genes annotated in the GO DNA-repair branch
# ------------------------------------------------------------------------------
# The dumbbell plot compares RNA-seq and RNAPII ChIP-seq log2 fold changes for
# the genes in the descriptive DNA-repair overlap. Connected points represent
# the two measurements for each gene; this is not an enrichment test. A
# complementary table stores the plotted values and GO annotations.

if (nrow(DNA_repair_GO_gene_summary) > 0L) {
  
  DNA_repair_dumbbell_data <- DNA_repair_GO_gene_summary %>%
    dplyr::mutate(
      log2FC_difference = RNAPII_ChIPseq_log2FC - RNAseq_log2FC,
      absolute_log2FC_difference = abs(log2FC_difference)
    ) %>%
    dplyr::arrange(
      RNAPII_ChIPseq_log2FC,
      RNAseq_log2FC,
      gene_symbol
    ) %>%
    dplyr::mutate(
      gene_symbol = factor(
        gene_symbol,
        levels = gene_symbol
      )
    )
  
  DNA_repair_dumbbell_summary_table <- DNA_repair_dumbbell_data %>%
    dplyr::transmute(
      gene_id,
      gene_symbol = as.character(gene_symbol),
      RNAseq_log2FC,
      RNAseq_padj,
      RNAPII_ChIPseq_log2FC,
      RNAPII_ChIPseq_padj,
      log2FC_difference,
      absolute_log2FC_difference,
      n_DNA_repair_GO_terms,
      DNA_repair_GO_IDs,
      DNA_repair_GO_terms,
      GO_evidence_codes
    ) %>%
    dplyr::arrange(
      dplyr::desc(RNAPII_ChIPseq_log2FC),
      dplyr::desc(RNAseq_log2FC),
      gene_symbol
    )
  
  write_tsv_base(
    DNA_repair_dumbbell_summary_table,
    "GO_DNA_repair_concordant_UP_dumbbell_summary.tsv"
  )
  
  p_GO_DNA_repair_dumbbell <- ggplot2::ggplot(
    DNA_repair_dumbbell_data,
    ggplot2::aes(
      y = gene_symbol
    )
  ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = RNAseq_log2FC,
        xend = RNAPII_ChIPseq_log2FC,
        yend = gene_symbol
      ),
      linewidth = 0.9
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        x = RNAseq_log2FC,
        shape = "RNA-seq"
      ),
      size = 4.5
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        x = RNAPII_ChIPseq_log2FC,
        shape = "RNAPII ChIP-seq"
      ),
      size = 4.5
    ) +
    ggplot2::scale_shape_manual(
      name = NULL,
      values = c(
        "RNA-seq" = 16,
        "RNAPII ChIP-seq" = 17
      )
    ) +
    ggplot2::geom_vline(
      xintercept = 1,
      linetype = "dashed",
      linewidth = 0.7
    ) +
    ggplot2::labs(
      x = "log2 fold change",
      y = NULL
    ) +
    ggplot2::theme_classic(
      base_size = 20
    ) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(
        size = 18,
        colour = "black"
      ),
      axis.text.x = ggplot2::element_text(
        size = 16,
        colour = "black"
      ),
      axis.title.x = ggplot2::element_text(
        size = 18
      ),
      legend.position = "top",
      legend.text = ggplot2::element_text(
        size = 17
      ),
      legend.key.size = grid::unit(
        0.6,
        "cm"
      ),
      plot.margin = ggplot2::margin(
        t = 10,
        r = 15,
        b = 10,
        l = 15
      )
    )
  
  print(
    p_GO_DNA_repair_dumbbell
  )
  
  save_plot(
    p_GO_DNA_repair_dumbbell,
    "GO_DNA_repair_concordant_UP_dumbbell_plot",
    width = 8,
    height = 5.5
  )
  
  ggplot2::ggsave(
    filename = file.path(
      plot_dir,
      "GO_DNA_repair_concordant_UP_dumbbell_plot.jpg"
    ),
    plot = p_GO_DNA_repair_dumbbell,
    width = 8,
    height = 5.5,
    dpi = safe_png_dpi(8, 5.5),
    bg = "white"
  )
  
} else {
  
  DNA_repair_dumbbell_data <- NULL
  DNA_repair_dumbbell_summary_table <- NULL
  p_GO_DNA_repair_dumbbell <- NULL
}

# ------------------------------------------------------------------------------
# 18. GO BP treeplot using the same top 15 non-simplified terms
# ------------------------------------------------------------------------------

if (nrow(top15_go_bp) >= 2L) {
  
  # ---------------------------------------------------------------------------
  # 1. Prepare semantic information required for Wang similarity
  # ---------------------------------------------------------------------------
  
  GO_BP_semData <- GOSemSim::godata(
    annoDb = org.Hs.eg.db,
    ont = "BP",
    computeIC = FALSE
  )
  
  
  # ---------------------------------------------------------------------------
  # 2. Keep exactly the same 15 terms shown in the non-simplified GO plots
  # ---------------------------------------------------------------------------
  
  tree_object_top15 <- ego_bp
  
  tree_object_top15@result <- tree_object_top15@result %>%
    dplyr::filter(
      ID %in% top15_go_bp$ID
    ) %>%
    dplyr::arrange(
      p.adjust,
      pvalue,
      dplyr::desc(Count)
    )
  
  tree_object_top15@termsim <- matrix(
    numeric(0),
    nrow = 0,
    ncol = 0
  )
  
  
  # ---------------------------------------------------------------------------
  # 3. Recalculate semantic similarity
  # ---------------------------------------------------------------------------
  
  tree_similarity_top15 <- enrichplot::pairwise_termsim(
    tree_object_top15,
    method = "Wang",
    semData = GO_BP_semData
  )
  
  
  # ---------------------------------------------------------------------------
  # 4. Treeplot with original legends
  # ---------------------------------------------------------------------------
  # This object is used only to extract the original legends.
  
  p_go_bp_treeplot_with_legends <- enrichplot::treeplot(
    tree_similarity_top15,
    showCategory = nrow(tree_object_top15@result),
    color = "p.adjust",
    nCluster = 4,
    cluster_method = "ward.D",
    
    # Wrap long terms earlier and increase readability
    label_format = 25,
    fontsize_tiplab = 4.4,
    tiplab_offset = 0.40,
    
    fontsize_cladelab = 0.01,
    cladelab_offset = 0,
    extend = 0,
    
    hilight = TRUE,
    align = "both",
    
    # More horizontal room for GO term labels
    hexpand = 0.38
  ) +
    ggplot2::labs(
      title = "Semantic similarity tree of the 15 most significant enriched GO Biological Process terms"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5,
        size = 19,
        margin = ggplot2::margin(
          b = 12
        )
      ),
      legend.position = "right",
      plot.margin = ggplot2::margin(
        t = 10,
        r = 55,
        b = 5,
        l = 15
      )
    ) +
    ggplot2::coord_cartesian(
      clip = "off"
    )
  
  
  # ---------------------------------------------------------------------------
  # 5. Visible treeplot without automatic legends
  # ---------------------------------------------------------------------------
  
  p_go_bp_treeplot_only <- p_go_bp_treeplot_with_legends +
    ggplot2::theme(
      legend.position = "none",
      plot.margin = ggplot2::margin(
        t = 10,
        r = 55,
        b = 5,
        l = 15
      )
    )
  
  
  # ---------------------------------------------------------------------------
  # 6. Extract original legends from treeplot
  # ---------------------------------------------------------------------------
  
  original_tree_legend_box <- cowplot::get_legend(
    p_go_bp_treeplot_with_legends
  )
  
  if (is.null(original_tree_legend_box)) {
    stop(
      "The treeplot legends could not be extracted.",
      call. = FALSE
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # 7. Helper function to identify legends by title
  # ---------------------------------------------------------------------------
  
  grob_contains_label <- function(grob_object, target_labels) {
    
    labels_found <- character()
    
    collect_grob_labels <- function(x) {
      
      if (inherits(x, "text")) {
        labels_found <<- c(
          labels_found,
          as.character(x$label)
        )
      }
      
      if (!is.null(x$grobs)) {
        lapply(
          x$grobs,
          collect_grob_labels
        )
      }
      
      if (!is.null(x$children)) {
        lapply(
          x$children,
          collect_grob_labels
        )
      }
      
      invisible(NULL)
    }
    
    collect_grob_labels(grob_object)
    
    any(target_labels %in% labels_found)
  }
  
  
  # ---------------------------------------------------------------------------
  # 8. Arrange Count and p.adjust legends horizontally
  # ---------------------------------------------------------------------------
  # The two original guides are placed side by side:
  #   Count | p.adjust
  #
  # They are not reconstructed, so their scales remain exactly those generated
  # by treeplot().
  
  # ---------------------------------------------------------------------------
  # Arrange Count and p.adjust legends horizontally with manual vertical alignment
  # ---------------------------------------------------------------------------
  
  # Retain only the original Count and p.adjust guides. The automatic
  # cluster-group guide is excluded because cluster names are displayed below.
  selected_tree_guides <- original_tree_legend_box$grobs[
    vapply(
      original_tree_legend_box$grobs,
      grob_contains_label,
      logical(1),
      target_labels = c("Count", "p.adjust")
    )
  ]

  if (length(selected_tree_guides) == 0L) {
    stop(
      "The Count and p.adjust guides could not be extracted.",
      call. = FALSE
    )
  }

  # Identify the two guides explicitly by their titles
  count_guide <- selected_tree_guides[
    vapply(
      selected_tree_guides,
      grob_contains_label,
      logical(1),
      target_labels = "Count"
    )
  ][[1]]
  
  padjust_guide <- selected_tree_guides[
    vapply(
      selected_tree_guides,
      grob_contains_label,
      logical(1),
      target_labels = "p.adjust"
    )
  ][[1]]
  
  
  # Convert each original guide into an independent plot
  p_count_guide <- cowplot::ggdraw() +
    cowplot::draw_grob(
      count_guide,
      x = 0,
      y = 0,
      width = 1,
      height = 1
    )
  
  p_padjust_guide <- cowplot::ggdraw() +
    cowplot::draw_grob(
      padjust_guide,
      x = 0,
      y = 0,
      width = 1,
      height = 1
    )
  
  
  # Place both guides manually inside the same legend panel
  tree_statistics_legend <- cowplot::ggdraw() +
    
    # Count legend
    cowplot::draw_plot(
      p_count_guide,
      x = 0.00,
      y = 0.1,
      width = 0.46,
      height = 0.965
    ) +
    
    # p.adjust legend
    cowplot::draw_plot(
      p_padjust_guide,
      x = 0.46,
      y = 0.00,
      width = 0.54,
      height = 1.00
    )
  
  
  # ---------------------------------------------------------------------------
  # 9. Place the statistical legends inside the upper-left empty area
  # ---------------------------------------------------------------------------
  #
  # The legend panel is made sufficiently wide and high to avoid clipping the
  # Count circles or the p.adjust colour bar.
  
  p_tree_with_internal_legend <- cowplot::ggdraw() +
    
    cowplot::draw_plot(
      p_go_bp_treeplot_only,
      x = 0,
      y = 0,
      width = 1,
      height = 1
    ) +
    
    cowplot::draw_plot(
      tree_statistics_legend,
      x = 0.035,
      y = 0.74,
      width = 0.15,
      height = 0.1
    )
  
  
  # ---------------------------------------------------------------------------
  # 10. Manual cluster legend below
  # ---------------------------------------------------------------------------
  
  manual_cluster_legend_data <- tibble::tibble(
    cluster = factor(
      c(
        "General stress and stimulus responses",
        "Regulation of programmed cell death",
        "Responses to bacterial stimuli",
        "p53-mediated DNA damage response and apoptosis"
      ),
      levels = c(
        "General stress and stimulus responses",
        "Regulation of programmed cell death",
        "Responses to bacterial stimuli",
        "p53-mediated DNA damage response and apoptosis"
      )
    ),
    x = c(1, 2, 3, 4),
    y = 1
  )
  
  p_manual_cluster_legend <- ggplot2::ggplot(
    manual_cluster_legend_data,
    ggplot2::aes(
      x = x,
      y = y,
      fill = cluster
    )
  ) +
    ggplot2::geom_tile(
      width = 0.20,
      height = 0.20
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = stringr::str_wrap(
          as.character(cluster),
          width = 28
        )
      ),
      nudge_y = -0.28,
      size = 6,
      fontface = "bold",
      lineheight = 0.9
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "General stress and stimulus responses" = "#00BFC4",
        "Regulation of programmed cell death" = "#C77CFF",
        "Responses to bacterial stimuli" = "#F8766D",
        "p53-mediated DNA damage response and apoptosis" = "#7CAE00"
      )
    ) +
    ggplot2::coord_cartesian(
      xlim = c(0.5, 4.5),
      ylim = c(0.45, 1.15),
      clip = "off"
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "none",
      plot.margin = ggplot2::margin(
        t = 0,
        r = 10,
        b = 10,
        l = 10
      )
    )
  
  
  # ---------------------------------------------------------------------------
  # 11. Combine top treeplot and bottom cluster legend
  # ---------------------------------------------------------------------------
  
  p_go_bp_treeplot_final <- cowplot::plot_grid(
    p_tree_with_internal_legend,
    p_manual_cluster_legend,
    ncol = 1,
    rel_heights = c(0.86, 0.14)
  )
  
  
  # ---------------------------------------------------------------------------
  # 12. Display and save
  # ---------------------------------------------------------------------------
  
  print(p_go_bp_treeplot_final)
  
  ggplot2::ggsave(
    filename = file.path(
      plot_dir,
      "GO_BP_top15_non_simplified_treeplot_final.pdf"
    ),
    plot = p_go_bp_treeplot_final,
    width = 15.5,
    height = 9,
    device = grDevices::cairo_pdf,
    bg = "white"
  )
  
  ggplot2::ggsave(
    filename = file.path(
      plot_dir,
      "GO_BP_top15_non_simplified_treeplot_final.png"
    ),
    plot = p_go_bp_treeplot_final,
    width = 15.5,
    height = 9,
    dpi = 600,
    bg = "white"
  )
  
  message(
    "Final GO BP treeplot saved in:\n",
    normalizePath(plot_dir)
  )
  
} else {
  
  message(
    "Fewer than two GO BP terms are available; ",
    "the treeplot was not generated."
  )
}