# ==============================================================================
# Main ChIP-seq quality control workflow
# ==============================================================================
# This module collects quality control metrics at different stages of the ChIP-seq
# processing workflow, from raw reads to normalized genome-wide signal tracks.
#
# QC modules:
#   1. Read-level QC
#   2. BAM-level QC
#   3. Alignment, library and enrichment QC
#   4. BAM-level reproducibility QC
#   5. Normalized signal reproducibility QC
#   6. Integrated MultiQC report
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Read-level QC
# ------------------------------------------------------------------------------
# FastQC is run before and after adapter trimming to assess raw sequencing quality
# and the effect of the trimming step.

# ------------------------------------------------------------------------------
# 1.1 Raw FASTQ quality control
# ------------------------------------------------------------------------------
# FastQC is used to assess the quality of raw sequencing reads before trimming.

rule fastqc_raw:
    input:
        get_raw_reads
    output:
        directory(os.path.join(QC_ROOT, "fastqc", "raw", "{sample}"))
    params:
        job_name = f"{config['GSE']}_{{sample}}_FASTQC_RAW"
    resources:
        time_min = int(custom_param("fastqc_time", 60)),
        cpus = int(custom_param("fastqc_cpus", 2)),
        mem_gb = int(custom_param("fastqc_mem_gb", 8))
    log:
        os.path.join(LOG_ROOT, "fastqc", "raw", "FastQC_raw.{sample}.log")
    benchmark:
        os.path.join(LOG_ROOT, "fastqc", "raw", "benchmark", "FastQC_raw.{sample}.benchmark")
    conda:
        os.path.join(ENV_DIR, "fastqc.yaml")
    shell:
        """
        mkdir -p {QC_ROOT}/fastqc/raw/{wildcards.sample} \
                 {LOG_ROOT}/fastqc/raw/benchmark

        fastqc \
          --threads {resources.cpus} \
          --outdir {QC_ROOT}/fastqc/raw/{wildcards.sample} \
          {input} \
          &> {log}
        """

# ------------------------------------------------------------------------------
# 1.2 Trimmed FASTQ quality control
# ------------------------------------------------------------------------------
# FastQC is used to assess read quality after adapter and quality trimming.

rule fastqc_trimmed:
    input:
        rules.trimgalore_paired.output if LAYOUT == "PAIRED" else rules.trimgalore_single.output
    output:
        directory(os.path.join(QC_ROOT, "fastqc", "trimmed", "{sample}"))
    params:
        job_name = f"{config['GSE']}_{{sample}}_FASTQC_TRIMMED"
    resources:
        time_min = int(custom_param("fastqc_time", 60)),
        cpus = int(custom_param("fastqc_cpus", 2)),
        mem_gb = int(custom_param("fastqc_mem_gb", 8))
    log:
        os.path.join(LOG_ROOT, "fastqc", "trimmed", "FastQC_trimmed.{sample}.log")
    benchmark:
        os.path.join(LOG_ROOT, "fastqc", "trimmed", "benchmark", "FastQC_trimmed.{sample}.benchmark")
    conda:
        os.path.join(ENV_DIR, "fastqc.yaml")
    shell:
        """
        mkdir -p {QC_ROOT}/fastqc/trimmed/{wildcards.sample} \
                 {LOG_ROOT}/fastqc/trimmed/benchmark

        fastqc \
          --threads {resources.cpus} \
          --outdir {QC_ROOT}/fastqc/trimmed/{wildcards.sample} \
          {input} \
          &> {log}
        """

# ------------------------------------------------------------------------------
# 3. BAM alignment quality control
# ------------------------------------------------------------------------------
# samtools is used to generate basic alignment statistics from final BAM files.

rule samtools_alignment_qc:
    input:
        bam = rules.process_bam.output.bam,
        bai = rules.process_bam.output.bai
    output:
        flagstat = os.path.join(QC_ROOT, "alignment", "samtools", "{sample}.flagstat.txt"),
        idxstats = os.path.join(QC_ROOT, "alignment", "samtools", "{sample}.idxstats.txt"),
        stats = os.path.join(QC_ROOT, "alignment", "samtools", "{sample}.stats.txt")
    params:
        job_name = f"{config['GSE']}_{{sample}}_SAMTOOLS_QC"
    resources:
        time_min = int(custom_param("samtools_qc_time", 30)),
        cpus = int(custom_param("samtools_qc_cpus", 2)),
        mem_gb = int(custom_param("samtools_qc_mem_gb", 8))
    log:
        os.path.join(LOG_ROOT, "samtools_qc", "samtools_qc.{sample}.log")
    benchmark:
        os.path.join(LOG_ROOT, "samtools_qc", "benchmark", "samtools_qc.{sample}.benchmark")
    conda:
        os.path.join(ENV_DIR, "samtools.yaml")
    shell:
        """
        mkdir -p {QC_ROOT}/alignment/samtools \
                 {LOG_ROOT}/samtools_qc/benchmark

        samtools flagstat \
          -@ {resources.cpus} \
          {input.bam} \
          > {output.flagstat} \
          2> {log}

        samtools idxstats \
          {input.bam} \
          > {output.idxstats} \
          2>> {log}

        samtools stats \
          -@ {resources.cpus} \
          {input.bam} \
          > {output.stats} \
          2>> {log}
        """
# ------------------------------------------------------------------------------
# 4. BAM-level reproducibility QC
# ------------------------------------------------------------------------------
# Correlation and PCA analyses are performed on genome-wide read-count matrices
# computed from the final BAM files. This provides an overview of sample similarity
# before coverage-track normalization.

rule plot_correlation_multibam:
    input:
        npz = os.path.join(
            OUT_ROOT,
            "multiBamSummary",
            f"{config['GSE']}_multiBamSummary.npz"
        )
    output:
        heatmap = os.path.join(
            QC_ROOT,
            "reproducibility",
            "correlation",
            f"{config['GSE']}.multiBamSummary.spearman_heatmap.pdf"
        ),
        matrix = os.path.join(
            QC_ROOT,
            "reproducibility",
            "correlation",
            f"{config['GSE']}.multiBamSummary.spearman_matrix.tab"
        )
    params:
        job_name = f"{config['GSE']}_PLOT_CORR_MULTIBAM"
    resources:
        time_min = int(custom_param("plot_correlation_multibam_time", 30)),
        cpus = int(custom_param("plot_correlation_multibam_cpus", 1)),
        mem_gb = int(custom_param("plot_correlation_multibam_mem_gb", 16))
    log:
        os.path.join(
            LOG_ROOT,
            "Quality_Control",
            f"plotCorrelation.multiBamSummary.{config['GSE']}.log"
        )
    conda:
        os.path.join(ENV_DIR, "deeptools.yaml")
    shell:
        """
        mkdir -p {QC_ROOT}/reproducibility/correlation \
                 {LOG_ROOT}/Quality_Control

        plotCorrelation \
          -in {input.npz} \
          --corMethod spearman \
          --skipZeros \
          --whatToPlot heatmap \
          --plotNumbers \
          -o {output.heatmap} \
          --outFileCorMatrix {output.matrix} \
          &> {log}
        """


rule plot_pca_multibam:
    input:
        npz = os.path.join(
            OUT_ROOT,
            "multiBamSummary",
            f"{config['GSE']}_multiBamSummary.npz"
        )
    output:
        pdf = os.path.join(
            QC_ROOT,
            "reproducibility",
            "PCA",
            f"{config['GSE']}.multiBamSummary.PCA.pdf"
        ),
        table = os.path.join(
            QC_ROOT,
            "reproducibility",
            "PCA",
            f"{config['GSE']}.multiBamSummary.PCA.tab"
        )
    params:
        job_name = f"{config['GSE']}_PLOT_PCA_MULTIBAM"
    resources:
        time_min = int(custom_param("plot_pca_multibam_time", 30)),
        cpus = int(custom_param("plot_pca_multibam_cpus", 1)),
        mem_gb = int(custom_param("plot_pca_multibam_mem_gb", 16))
    log:
        os.path.join(
            LOG_ROOT,
            "Quality_Control",
            f"plotPCA.multiBamSummary.{config['GSE']}.log"
        )
    conda:
        os.path.join(ENV_DIR, "deeptools.yaml")
    shell:
        """
        mkdir -p {QC_ROOT}/reproducibility/PCA \
                 {LOG_ROOT}/Quality_Control

        plotPCA \
          -in {input.npz} \
          -o {output.pdf} \
          --outFileNameData {output.table} \
          &> {log}
        """

# ------------------------------------------------------------------------------
# 5. Normalized signal reproducibility QC
# ------------------------------------------------------------------------------
# Correlation and PCA analyses are performed on normalized bigWig signal matrices.
# This evaluates sample similarity after applying the selected coverage
# normalization method.

rule plot_correlation_multibigwig:
    input:
        npz = os.path.join(
            OUT_ROOT,
            "multiBigwigSummary",
            f"{config['GSE']}_multiBigwigSummary_{BW_NORMALIZATION}.npz"
        )
    output:
        heatmap = os.path.join(
            QC_ROOT,
            "reproducibility",
            "normalized_signal",
            BW_NORMALIZATION,
            "correlation",
            f"{config['GSE']}.{BW_NORMALIZATION}.spearman_heatmap.pdf"
        ),
        matrix = os.path.join(
            QC_ROOT,
            "reproducibility",
            "normalized_signal",
            BW_NORMALIZATION,
            "correlation",
            f"{config['GSE']}.{BW_NORMALIZATION}.spearman_matrix.tab"
        )
    params:
        job_name = f"{config['GSE']}_PLOT_CORR_MULTIBW_{BW_NORMALIZATION}"
    resources:
        time_min = int(custom_param("plot_correlation_multibw_time", 30)),
        cpus = int(custom_param("plot_correlation_multibw_cpus", 1)),
        mem_gb = int(custom_param("plot_correlation_multibw_mem_gb", 16))
    log:
        os.path.join(
            LOG_ROOT,
            "Quality_Control",
            f"plotCorrelation.multiBigwigSummary.{BW_NORMALIZATION}.{config['GSE']}.log"
        )
    conda:
        os.path.join(ENV_DIR, "deeptools.yaml")
    shell:
        """
        mkdir -p {QC_ROOT}/reproducibility/normalized_signal/{BW_NORMALIZATION}/correlation \
                 {LOG_ROOT}/Quality_Control

        plotCorrelation \
          -in {input.npz} \
          --corMethod spearman \
          --skipZeros \
          --whatToPlot heatmap \
          --plotNumbers \
          -o {output.heatmap} \
          --outFileCorMatrix {output.matrix} \
          &> {log}
        """


rule plot_pca_multibigwig:
    input:
        npz = os.path.join(
            OUT_ROOT,
            "multiBigwigSummary",
            f"{config['GSE']}_multiBigwigSummary_{BW_NORMALIZATION}.npz"
        )
    output:
        pdf = os.path.join(
            QC_ROOT,
            "reproducibility",
            "normalized_signal",
            BW_NORMALIZATION,
            "PCA",
            f"{config['GSE']}.{BW_NORMALIZATION}.PCA.pdf"
        ),
        table = os.path.join(
            QC_ROOT,
            "reproducibility",
            "normalized_signal",
            BW_NORMALIZATION,
            "PCA",
            f"{config['GSE']}.{BW_NORMALIZATION}.PCA.tab"
        )
    params:
        job_name = f"{config['GSE']}_PLOT_PCA_MULTIBW_{BW_NORMALIZATION}"
    resources:
        time_min = int(custom_param("plot_pca_multibw_time", 30)),
        cpus = int(custom_param("plot_pca_multibw_cpus", 1)),
        mem_gb = int(custom_param("plot_pca_multibw_mem_gb", 16))
    log:
        os.path.join(
            LOG_ROOT,
            "Quality_Control",
            f"plotPCA.multiBigwigSummary.{BW_NORMALIZATION}.{config['GSE']}.log"
        )
    conda:
        os.path.join(ENV_DIR, "deeptools.yaml")
    shell:
        """
        mkdir -p {QC_ROOT}/reproducibility/normalized_signal/{BW_NORMALIZATION}/PCA \
                 {LOG_ROOT}/Quality_Control

        plotPCA \
          -in {input.npz} \
          -o {output.pdf} \
          --outFileNameData {output.table} \
          &> {log}
        """

# ------------------------------------------------------------------------------
# 6. MultiQC report
# ------------------------------------------------------------------------------
# MultiQC aggregates quality-control outputs from FastQC, Trim Galore, STAR,
# samtools, featureCounts and deepTools into a single HTML report.

rule multiqc:
    input:
        raw_fastqc = expand(
            os.path.join(QC_ROOT, "fastqc", "raw", "{sample}"),
            sample = SAMPLES
        ),
        trimmed_fastqc = expand(
            os.path.join(QC_ROOT, "fastqc", "trimmed", "{sample}"),
            sample = SAMPLES
        ),
        star_logs = expand(
            os.path.join(OUT_ROOT, ALIGN_DIR, "{sample}.Log.final.out"),
            sample = SAMPLES
        ),
        flagstat = expand(
            os.path.join(QC_ROOT, "alignment", "samtools", "{sample}.flagstat.txt"),
            sample = SAMPLES
        ),
        idxstats = expand(
            os.path.join(QC_ROOT, "alignment", "samtools", "{sample}.idxstats.txt"),
            sample = SAMPLES
        ),
        stats = expand(
            os.path.join(QC_ROOT, "alignment", "samtools", "{sample}.stats.txt"),
            sample = SAMPLES
        ),
        exon_summary = os.path.join(
            OUT_ROOT,
            "counts",
            "exon",
            "featureCounts_exon.txt.summary"
        ),
        gene_body_summary = os.path.join(
            OUT_ROOT,
            "counts",
            "gene_body",
            "featureCounts_gene_body.txt.summary"
        ),
        multibam_npz = os.path.join(
            OUT_ROOT,
            "multiBamSummary",
            f"{config['GSE']}_multiBamSummary.npz"
        ),
        multibigwig_npz = os.path.join(
            OUT_ROOT,
            "multiBigwigSummary",
            f"{config['GSE']}_multiBigwigSummary_{BW_NORMALIZATION}.npz"
        )
    output:
        html = os.path.join(OUT_ROOT, "multiqc", f"{config['GSE']}_multiqc_report.html")
    params:
        job_name = f"{config['GSE']}_MULTIQC"
    resources:
        time_min = int(custom_param("multiqc_time", 30)),
        cpus = int(custom_param("multiqc_cpus", 1)),
        mem_gb = int(custom_param("multiqc_mem_gb", 8))
    log:
        os.path.join(LOG_ROOT, "multiqc", f"MultiQC.{config['GSE']}.log")
    benchmark:
        os.path.join(LOG_ROOT, "multiqc", "benchmark", f"MultiQC.{config['GSE']}.benchmark")
    conda:
        os.path.join(ENV_DIR, "multiqc.yaml")
    shell:
        """
        mkdir -p {OUT_ROOT}/multiqc \
                 {LOG_ROOT}/multiqc/benchmark

        multiqc \
          {OUT_ROOT} \
          {QC_ROOT} \
          {LOG_ROOT} \
          --outdir {OUT_ROOT}/multiqc \
          --filename {config[GSE]}_multiqc_report.html \
          &> {log}
        """