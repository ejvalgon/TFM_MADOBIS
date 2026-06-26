# TFM_MADOBIS

Computational workflows and downstream analysis scripts associated with the Master's Thesis:

**Camptothecin-induced topoisomerase I poisoning drives gene-length-dependent RNA polymerase II redistribution and stress response**

This repository contains the code used to process RNAPII ChIP-seq and bulk RNA-seq data and to perform the downstream analyses reported in the thesis. The project focuses on the genome-wide response of RNA polymerase II (RNAPII) to short-term camptothecin (CPT) treatment, using spike-in-normalized ChIP-seq data and integrating these results with an independent CPT RNA-seq dataset.

The repository is organized into three main components:

1.  a modular Snakemake workflow for RNAPII ChIP-seq processing;
2.  a modular Snakemake workflow for bulk RNA-seq processing;
3.  ordered R scripts for downstream genomic, statistical and functional analyses.

Raw sequencing data and large processed files are not distributed in this repository. Configuration files contain placeholder paths that must be adapted before running the workflows.

------------------------------------------------------------------------

## Repository structure

``` text
TFM_MADOBIS/
├── CHIPseq_pipeline/
│   ├── Snakefile_ChipSeq
│   ├── config/
│   │   └── config_CHIPseq_module_comentado.yaml
│   ├── envs/
│   │   ├── *_environment.yaml
│   │   ├── *_packages.txt
│   │   ├── *_packages.json
│   │   ├── export_summary.tsv
│   │   └── workflow_manager_versions.tsv
│   ├── logs/
│   ├── sbatch/
│   │   ├── controller.sh
│   │   └── star_index_generator.sh
│   ├── scripts/
│   │   └── spikein_scale_factors.R
│   ├── slurm_config/
│   │   ├── config.yaml
│   │   ├── parseJobID.sh
│   │   └── slurm_jobscript.sh
│   └── workflow/
│       ├── main_workflow.smk
│       ├── main_qc.smk
│       └── peak_calling.smk
│
├── RNAseq_pipeline/
│   ├── Snakefile_RNAseq
│   ├── config/
│   │   └── config_RNAseq.yaml
│   ├── envs/
│   │   ├── *_environment.yaml
│   │   ├── *_packages.txt
│   │   ├── *_packages.json
│   │   └── workflow_manager_versions.tsv
│   ├── logs/
│   ├── sbatch/
│   │   ├── controller.sh
│   │   └── RNAseq_star_index_generator.sh
│   ├── slurm_config/
│   │   ├── config.yaml
│   │   ├── parseJobID.sh
│   │   └── slurm_jobscript.sh
│   └── workflow/
│       ├── workflow_RNAseq.smk
│       └── qc_RNAseq.smk
│
└── Downstream_analysis/
    ├── 00_prepare_gene_annotations.R
    ├── 01_build_metagene_objects_final.R
    ├── 02_plot_metagene_profiles_final_vf.R
    ├── 03_plot_metagene_heatmaps_final_vf.R
    ├── 04_count_RNAPII_ChIPseq_internal_gene_bodies_vf.R
    ├── 05_DESeq2_ChIPseq_vf.R
    ├── 06_DESeq2_RNAseq_vf.R
    ├── 07_define_common_gene_universe_vf.R
    ├── 08_analyze_RNAseq_RNAPII_ChIPseq_correspondence_vf.R
    ├── 09_analyze_RNAPII_gene_length_vf.R
    ├── 10_functional_enrichment_vf.R
    ├── 11_DESeq2_Compare_spikein_vs_NOspikein_vf.R
    ├── 12_PCA.R
    ├── 13_titration_assay.R
    └── 14_Compare_spikein_vs_NOspikein.R
```

------------------------------------------------------------------------

## Project overview

The analysis was designed to quantify how CPT-induced TOP1 poisoning affects chromatin-associated RNAPII and how these changes relate to gene length and RNA abundance.

The RNAPII ChIP-seq workflow processes paired-end sequencing data generated from RPE-1 cells treated with DMSO or 10 µM CPT for 1 hour. A yeast chromatin spike-in was used to support quantitative normalization between samples. The RNA-seq workflow processes an independent bulk RNA-seq dataset from HeLa cells treated with DMSO or 10 µM CPT for 24 hours. The downstream R analyses integrate both datasets at the gene level.

The main analytical steps are:

-   generation of processed BAM files and genome-wide bigWig tracks;
-   spike-in-based normalization of RNAPII ChIP-seq signal;
-   quantification of RNAPII occupancy in internal gene-body regions;
-   differential RNAPII occupancy and differential RNA abundance analysis with DESeq2;
-   definition of a common active-gene universe for cross-platform comparison;
-   analysis of the relationship between RNAPII response and gene length;
-   functional enrichment analysis of genes concordantly increased in RNAPII occupancy and RNA abundance.

------------------------------------------------------------------------

## RNAPII ChIP-seq workflow

The ChIP-seq pipeline is implemented in `CHIPseq_pipeline/` and is controlled by:

``` text
CHIPseq_pipeline/Snakefile_ChipSeq
CHIPseq_pipeline/config/config_CHIPseq_module_comentado.yaml
```

The main Snakefile imports three workflow modules:

``` text
workflow/main_workflow.smk
workflow/main_qc.smk
workflow/peak_calling.smk
```

### Main processing steps

The core workflow performs:

1.  adapter and quality trimming with Trim Galore;
2.  read alignment with STAR or Bowtie2, depending on the configuration;
3.  BAM sorting, indexing and duplicate removal;
4.  optional splitting of reads into main-genome and spike-in BAM files;
5.  spike-in read counting with `samtools idxstats`;
6.  calculation of spike-in scaling factors with `scripts/spikein_scale_factors.R`;
7.  generation of normalized bigWig tracks with deepTools `bamCoverage`;
8.  generation of genome-wide summary matrices with `multiBamSummary` and `multiBigwigSummary`.

When `spikeIN: "yes"` is used, reads are aligned to a combined main-organism plus spike-in reference. The workflow then separates BAM files using the configured genome prefixes. In the provided configuration, `main_prefix` is set to `hs_` and `spikein_prefix` to `sc_`.

### Quality-control outputs

The QC module generates read-, alignment-, library- and signal-level quality-control files, including:

-   FastQC reports before and after trimming;
-   SAMtools `flagstat`, `idxstats` and `stats`;
-   Picard alignment and insert-size metrics;
-   phantompeakqualtools cross-correlation metrics;
-   blacklist-fraction estimates;
-   deepTools fingerprint, correlation and PCA plots;
-   an integrated MultiQC report.

### Peak-calling module

A peak-calling module is included in `workflow/peak_calling.smk`. It supports MACS3 peak calling, consensus and union peak generation, peak-level read counting, FRiP calculation and peak-centered signal heatmaps.

In the provided ChIP-seq configuration, this module is disabled by default:

``` yaml
calling_peaks: "no"
```

This reflects the RNAPII occupancy-focused analysis used in the thesis, where the main downstream quantification was performed over annotated internal gene-body regions rather than peak calls.

------------------------------------------------------------------------

## RNA-seq workflow

The RNA-seq pipeline is implemented in `RNAseq_pipeline/` and is controlled by:

``` text
RNAseq_pipeline/Snakefile_RNAseq
RNAseq_pipeline/config/config_RNAseq.yaml
```

The main Snakefile imports two modules:

``` text
workflow/workflow_RNAseq.smk
workflow/qc_RNAseq.smk
```

### Main processing steps

The RNA-seq workflow performs:

1.  adapter and quality trimming with Trim Galore;
2.  splice-aware alignment with STAR;
3.  BAM read-group annotation, sorting and indexing;
4.  gene-level quantification with featureCounts;
5.  exon-level and gene-body count generation;
6.  CPM-normalized bigWig generation with deepTools `bamCoverage`;
7.  genome-wide summary matrix generation with `multiBamSummary` and `multiBigwigSummary`.

The provided configuration is set up for paired-end, unstranded RNA-seq data:

``` yaml
layout: "PAIRED"
featurecounts_strand: 0
bamcoverage_normalization: "CPM"
```

The STAR maximum intron length is explicitly configurable. In the provided configuration it is set to:

``` yaml
star_align_intron_max: 1500000
```

### Quality-control outputs

The RNA-seq QC module generates:

-   FastQC reports before and after trimming;
-   SAMtools alignment statistics;
-   featureCounts summary files;
-   deepTools correlation and PCA plots from BAM and bigWig summary matrices;
-   an integrated MultiQC report.

------------------------------------------------------------------------

## Downstream R analyses

The `Downstream_analysis/` directory contains ordered R scripts. They are intended to be run sequentially after the corresponding Snakemake workflows have generated the required BAM files, bigWig tracks, count matrices and spike-in scaling factors.

The scripts use relative input and output paths matching the local analysis structure used for the thesis. If the repository is cloned into a different layout, the file paths defined at the beginning of each script should be edited before execution.

### Script summary

| Script | Purpose |
|------------------------------------|------------------------------------|
| `00_prepare_gene_annotations.R` | Prepares protein-coding gene annotations and strand-aware internal gene-body regions. |
| `01_build_metagene_objects_final.R` | Builds RNAPII metagene RDS objects from spike-in-normalized bigWig files. |
| `02_plot_metagene_profiles_final_vf.R` | Generates metagene profile plots from precomputed RDS objects. |
| `03_plot_metagene_heatmaps_final_vf.R` | Generates RNAPII heatmaps and replicate-concordance statistics. |
| `04_count_RNAPII_ChIPseq_internal_gene_bodies_vf.R` | Counts paired-end RNAPII ChIP-seq fragments in internal gene-body regions. |
| `05_DESeq2_ChIPseq_vf.R` | Performs differential RNAPII occupancy analysis using spike-in-derived size factors. |
| `06_DESeq2_RNAseq_vf.R` | Performs RNA-seq differential-expression analysis with DESeq2. |
| `07_define_common_gene_universe_vf.R` | Defines the common active-gene universe used for RNA-seq and RNAPII ChIP-seq integration. |
| `08_analyze_RNAseq_RNAPII_ChIPseq_correspondence_vf.R` | Compares RNAPII occupancy changes with RNA abundance changes. |
| `09_analyze_RNAPII_gene_length_vf.R` | Tests the relationship between RNAPII occupancy changes and internal gene-body length. |
| `10_functional_enrichment_vf.R` | Performs GO and KEGG over-representation analysis on concordantly increased genes. |
| `11_DESeq2_Compare_spikein_vs_NOspikein_vf.R` | Compares RNAPII differential-occupancy classifications with and without spike-in size factors. |
| `12_PCA.R` | Generates PCA plots from deepTools matrices and DESeq2-normalized ChIP-seq counts. |
| `13_titration_assay.R` | Tests the linear relationship between spike-in chromatin proportion and recovered spike-in reads. |
| `14_Compare_spikein_vs_NOspikein.R` | Additional comparison of spike-in and non-spike-in RNAPII ChIP-seq normalization strategies. |

### Main analysis thresholds

The downstream analysis uses the following main thresholds and definitions:

-   protein-coding genes on standard chromosomes;
-   minimum gene length of 3 kb for RNAPII metagene and gene-body analyses;
-   internal gene-body regions generated by removing 2.5 kb from the TSS side and 0.5 kb from the TES side;
-   low-count filtering before DESeq2 analysis: at least 100 total counts or fragments across all samples;
-   differential categories defined with adjusted *P* value ≤ 0.05 and absolute log2 fold change \> 1;
-   active RNA-seq genes defined by mean DMSO expression ≥ 1 RPKM;
-   common active-gene universe defined as the intersection of filtered RNA-seq and RNAPII ChIP-seq genes after applying matched basal-expression and basal-occupancy filters.

------------------------------------------------------------------------

## Input data expected by the workflows

The workflows expect raw FASTQ files and reference files to be provided by the user. The configuration files should be edited before execution.

### ChIP-seq inputs

The ChIP-seq configuration expects:

-   raw FASTQ files;
-   a STAR or Bowtie2 index for the main genome;
-   if spike-in normalization is enabled, a combined main-genome plus spike-in reference index;
-   main-genome and spike-in contig prefixes for BAM splitting;
-   a reference FASTA matching the selected alignment index;
-   a blacklist BED file for the main organism;
-   optional metadata for the peak-calling module.

### RNA-seq inputs

The RNA-seq configuration expects:

-   raw FASTQ files;
-   a STAR genome index generated for RNA-seq;
-   a gene annotation GTF file for featureCounts;
-   a reference FASTA file for documentation or extension of the workflow.

FASTQ sample discovery supports paired-end files ending in either:

``` text
_R1_001.fastq.gz / _R2_001.fastq.gz
_R1.fastq.gz     / _R2.fastq.gz
```

The extension can be modified in the corresponding configuration file.

------------------------------------------------------------------------

## Running the workflows

Before running either workflow, edit the relevant YAML configuration file to point to the local raw FASTQ directory, output directory, log directory and reference files.

### Dry run

``` bash
cd CHIPseq_pipeline
snakemake -s Snakefile_ChipSeq --use-conda --cores 1 --dry-run
```

``` bash
cd RNAseq_pipeline
snakemake -s Snakefile_RNAseq --use-conda --cores 1 --dry-run
```

### Local execution

``` bash
cd CHIPseq_pipeline
snakemake -s Snakefile_ChipSeq --use-conda --cores <N>
```

``` bash
cd RNAseq_pipeline
snakemake -s Snakefile_RNAseq --use-conda --cores <N>
```

### SLURM execution

The repository includes SLURM-related files in:

``` text
CHIPseq_pipeline/sbatch/
CHIPseq_pipeline/slurm_config/
RNAseq_pipeline/sbatch/
RNAseq_pipeline/slurm_config/
```

These files were used to submit and manage workflow execution on a high-performance computing environment. Cluster parameters should be checked and adapted to the local scheduler before use.

------------------------------------------------------------------------

## Software and environments

The workflows were developed using Snakemake and Conda-managed software environments. The `envs/` directories contain environment YAML files and exported package lists for reproducibility.

Main tools used by the workflows include:

-   Snakemake;
-   Conda;
-   Trim Galore and Cutadapt;
-   FastQC;
-   STAR;
-   Bowtie2, supported by the ChIP-seq workflow when selected;
-   SAMtools;
-   Picard;
-   deepTools;
-   Subread featureCounts;
-   MACS3, for the optional ChIP-seq peak-calling module;
-   phantompeakqualtools;
-   MultiQC.

The downstream analysis was performed in R using packages including:

-   DESeq2;
-   rtracklayer;
-   GenomicRanges;
-   GenomeInfoDb;
-   EnrichedHeatmap;
-   ComplexHeatmap;
-   Rsubread;
-   ggplot2;
-   dplyr;
-   tidyr;
-   cowplot;
-   clusterProfiler;
-   enrichplot;
-   org.Hs.eg.db;
-   AnnotationDbi;
-   GO.db;
-   GOSemSim.

Exact package versions are documented in in the thesis supplementary material.

------------------------------------------------------------------------

## Notes on reproducibility

This repository is intended to document the computational analysis performed for the thesis and to provide reusable workflows for related RNAPII ChIP-seq and RNA-seq analyses.

Several paths in the configuration files and R scripts are intentionally project-specific or represented as placeholders. To reproduce the analysis in a different environment, update:

-   raw FASTQ paths;
-   output and log directories;
-   genome indexes;
-   FASTA and GTF annotation files;
-   blacklist files;
-   spike-in reference settings;
-   paths to count matrices, bigWig files and DESeq2 outputs used by the downstream R scripts.

Because raw sequencing data and large intermediate files are not included, the repository cannot be run end-to-end immediately after cloning without providing the required input data and references.

------------------------------------------------------------------------

## Data availability

The primary RNAPII ChIP-seq dataset was generated for the thesis project and is not included in this repository. The code, workflow definitions and exported computational environments are provided to support transparency and reproducibility of the analyses.

------------------------------------------------------------------------

## Citation

If this repository is used or adapted, please cite the associated Master's Thesis:

Valenzuela-González, E. J. **Camptothecin-induced topoisomerase I poisoning drives gene-length-dependent RNA polymerase II redistribution and stress response.** Master's Thesis, Máster Universitario en Análisis de Datos Ómicos y Biología de Sistemas, Universidad de Sevilla and Universidad Internacional de Andalucía.

Repository: `https://github.com/ejvalgon/TFM_MADOBIS`
