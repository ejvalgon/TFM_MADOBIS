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

rule fastqc_raw:
    input:
        get_raw_reads
    output:
        fastqc_dir = directory(os.path.join(QC_ROOT, "fastq", "raw_fastqc", "{sample}"))
    params:
        job_name = f"{config['GSE']}_{{sample}}_FASTQC_RAW"
    resources:
        time_min = int(custom_param("fastqc_raw_time", 30)),
        cpus = int(custom_param("fastqc_raw_cpus", 4)),
        mem_gb = int(custom_param("fastqc_raw_mem_gb", 16))
    log:
        os.path.join(LOG_ROOT, "fastq", "fastqc_raw.{sample}.log")
    benchmark:
        os.path.join(LOG_ROOT, "fastq", "benchmark", "fastqc_raw.{sample}.benchmark")
    conda:
        os.path.join(ENV_DIR, "trimgalore_samtools.yaml")
    shell:
        """
        mkdir -p {output.fastqc_dir} \
                 {LOG_ROOT}/fastq/benchmark

        fastqc \
          -t {resources.cpus} \
          -o {output.fastqc_dir} \
          {input} \
          &> {log}
        """

rule fastqc_trimmed:
    input:
        r1 = os.path.join(OUT_ROOT, "fastq_trimmed", "{sample}_val_1.fq.gz") if LAYOUT == "PAIRED" else os.path.join(OUT_ROOT, "fastq_trimmed", "{sample}_trimmed.fq.gz"),
        r2 = os.path.join(OUT_ROOT, "fastq_trimmed", "{sample}_val_2.fq.gz") if LAYOUT == "PAIRED" else []
    output:
        fastqc_dir = directory(os.path.join(QC_ROOT, "fastq", "trimmed_fastqc", "{sample}"))
    params:
        job_name = f"{config['GSE']}_{{sample}}_FASTQC_TRIMMED"
    resources:
        time_min = int(custom_param("fastqc_trimmed_time", 30)),
        cpus = int(custom_param("fastqc_trimmed_cpus", 4)),
        mem_gb = int(custom_param("fastqc_trimmed_mem_gb", 16))
    log:
        os.path.join(LOG_ROOT, "fastq", "fastqc_trimmed.{sample}.log")
    benchmark:
        os.path.join(LOG_ROOT, "fastq", "benchmark", "fastqc_trimmed.{sample}.benchmark")
    conda:
        os.path.join(ENV_DIR, "trimgalore_samtools.yaml")
    shell:
        """
        mkdir -p {output.fastqc_dir} \
                 {LOG_ROOT}/fastq/benchmark

        fastqc \
          -t {resources.cpus} \
          -o {output.fastqc_dir} \
          {input} \
          &> {log}
        """

# ------------------------------------------------------------------------------
# 2. BAM-level QC
# ------------------------------------------------------------------------------
# Basic alignment statistics are computed from the final processed BAM files using
# samtools. These metrics summarize mapping rates, read-level flags and the
# distribution of aligned reads across reference sequences.

rule samtools_qc:
    input:
        bam = os.path.join(OUT_ROOT, ALIGN_DIR, "{sample}.bam"),
        bai = os.path.join(OUT_ROOT, ALIGN_DIR, "{sample}.bam.bai")
    output:
        flagstat = os.path.join(QC_ROOT, "alignment", "samtools_flagstat", "{sample}.flagstat.txt"),
        stats    = os.path.join(QC_ROOT, "alignment", "samtools_stats", "{sample}.stats.txt"),
        idxstats = os.path.join(QC_ROOT, "alignment", "samtools_idxstats", "{sample}.idxstats.txt")
    params:
        job_name = f"{config['GSE']}_{{sample}}_SAMTOOLS_QC"
    resources:
        time_min = int(custom_param("samtools_qc_time", 20)),
        cpus = int(custom_param("samtools_qc_cpus", 2)),
        mem_gb = int(custom_param("samtools_qc_mem_gb", 16))
    log:
        os.path.join(LOG_ROOT, "Quality_Control", "samtools_qc.{sample}.log")
    benchmark:
        os.path.join(LOG_ROOT, "Quality_Control", "benchmark", "samtools_qc.{sample}.benchmark")
    conda:
        os.path.join(ENV_DIR, "samtools_picard.yaml")
    shell:
        """
        mkdir -p {QC_ROOT}/alignment/samtools_flagstat \
                 {QC_ROOT}/alignment/samtools_stats \
                 {QC_ROOT}/alignment/samtools_idxstats \
                 {LOG_ROOT}/Quality_Control/benchmark

        samtools flagstat -@ {resources.cpus} {input.bam} > {output.flagstat} 2> {log}
        samtools stats -@ {resources.cpus} {input.bam} > {output.stats} 2>> {log}
        samtools idxstats {input.bam} > {output.idxstats} 2>> {log}
        """

# ------------------------------------------------------------------------------
# 3. Alignment, library and enrichment QC
# ------------------------------------------------------------------------------
# Picard and deepTools-based metrics are used to evaluate alignment quality,
# library properties, ChIP enrichment and potential technical artifacts.
#
# For spike-in datasets, alignment metrics are computed both on the combined BAM
# and on the main-organism BAM after host/spike-in separation.

if spikeIN == "yes":

    rule picard_alignment_metrics_combined:
        input:
            bam = os.path.join(OUT_ROOT, ALIGN_DIR, "{sample}.bam")
        output:
            metrics = os.path.join(
                QC_ROOT,
                "alignment",
                "picard_alignment",
                "combined",
                "{sample}.combined.alignment_summary_metrics.txt"
            )
        params:
            job_name = f"{config['GSE']}_{{sample}}_PICARD_ALIGN_COMBINED",
            reference = config[ORGANISM]["reference_fasta_spike_in"]
        resources:
            time_min = int(custom_param("picard_alignment_time", 30)),
            cpus = int(custom_param("picard_alignment_cpus", 1)),
            mem_gb = int(custom_param("picard_alignment_mem_gb", 16))
        log:
            os.path.join(LOG_ROOT, "Quality_Control", "picard_alignment.combined.{sample}.log")
        benchmark:
            os.path.join(LOG_ROOT, "Quality_Control", "benchmark", "picard_alignment.combined.{sample}.benchmark")
        conda:
            os.path.join(ENV_DIR, "samtools_picard.yaml")
        shell:
            """
            mkdir -p {QC_ROOT}/alignment/picard_alignment/combined \
                     {LOG_ROOT}/Quality_Control/benchmark

            picard CollectAlignmentSummaryMetrics \
              R={params.reference} \
              I={input.bam} \
              O={output.metrics} \
              &> {log}
            """


    rule picard_alignment_metrics_main:
        input:
            bam = os.path.join(OUT_ROOT, "split_bam", "main", "{sample}.bam")
        output:
            metrics = os.path.join(
                QC_ROOT,
                "alignment",
                "picard_alignment",
                "main",
                "{sample}.main.alignment_summary_metrics.txt"
            )
        params:
            job_name = f"{config['GSE']}_{{sample}}_PICARD_ALIGN_MAIN",
            reference = config[ORGANISM]["reference_fasta"]
        resources:
            time_min = int(custom_param("picard_alignment_time", 30)),
            cpus = int(custom_param("picard_alignment_cpus", 1)),
            mem_gb = int(custom_param("picard_alignment_mem_gb", 16))
        log:
            os.path.join(LOG_ROOT, "Quality_Control", "picard_alignment.main.{sample}.log")
        benchmark:
            os.path.join(LOG_ROOT, "Quality_Control", "benchmark", "picard_alignment.main.{sample}.benchmark")
        conda:
            os.path.join(ENV_DIR, "samtools_picard.yaml")
        shell:
            """
            mkdir -p {QC_ROOT}/alignment/picard_alignment/main \
                     {LOG_ROOT}/Quality_Control/benchmark

            picard CollectAlignmentSummaryMetrics \
              R={params.reference} \
              I={input.bam} \
              O={output.metrics} \
              &> {log}
            """

else:

    rule picard_alignment_metrics:
        input:
            bam = os.path.join(OUT_ROOT, ALIGN_DIR, "{sample}.bam")
        output:
            metrics = os.path.join(
                QC_ROOT,
                "alignment",
                "picard_alignment",
                "{sample}.alignment_summary_metrics.txt"
            )
        params:
            job_name = f"{config['GSE']}_{{sample}}_PICARD_ALIGN",
            reference = config[ORGANISM]["reference_fasta"]
        resources:
            time_min = int(custom_param("picard_alignment_time", 30)),
            cpus = int(custom_param("picard_alignment_cpus", 1)),
            mem_gb = int(custom_param("picard_alignment_mem_gb", 16))
        log:
            os.path.join(LOG_ROOT, "Quality_Control", "picard_alignment.{sample}.log")
        benchmark:
            os.path.join(LOG_ROOT, "Quality_Control", "benchmark", "picard_alignment.{sample}.benchmark")
        conda:
            os.path.join(ENV_DIR, "samtools_picard.yaml")
        shell:
            """
            mkdir -p {QC_ROOT}/alignment/picard_alignment \
                     {LOG_ROOT}/Quality_Control/benchmark

            picard CollectAlignmentSummaryMetrics \
              R={params.reference} \
              I={input.bam} \
              O={output.metrics} \
              &> {log}
            """

if LAYOUT == "PAIRED":
    rule picard_insert_size_metrics:
        input:
            bam = os.path.join(OUT_ROOT, ALIGN_DIR, "{sample}.bam")
        output:
            metrics = os.path.join(QC_ROOT, "insert_size", "{sample}.insert_size_metrics.txt"),
            pdf = os.path.join(QC_ROOT, "insert_size", "{sample}.insert_size_histogram.pdf")
        params:
            job_name = f"{config['GSE']}_{{sample}}_INSERTSIZE"
        resources:
            time_min = int(custom_param("insert_size_time", 30)),
            cpus = int(custom_param("insert_size_cpus", 1)),
            mem_gb = int(custom_param("insert_size_mem_gb", 16))
        log:
            os.path.join(LOG_ROOT, "Quality_Control", "insert_size.{sample}.log")
        benchmark:
            os.path.join(LOG_ROOT, "Quality_Control", "benchmark", "insert_size.{sample}.benchmark")
        conda:
            os.path.join(ENV_DIR, "samtools_picard.yaml")
        shell:
            """
            mkdir -p {QC_ROOT}/insert_size \
                     {LOG_ROOT}/Quality_Control/benchmark

            picard CollectInsertSizeMetrics \
              I={input.bam} \
              O={output.metrics} \
              H={output.pdf} \
              M=0.5 \
              &> {log}
            """

# Cross-correlation metrics are used to estimate ChIP enrichment and fragment
# length consistency using phantompeakqualtools/SPP.
# Note: phantompeakqualtools is tag/read-based and was designed for single-end
# ChIP-seq alignments. For paired-end BAMs, it can be used as a complementary QC,
# but it does not use paired-end fragment information explicitly.
rule phantompeakqualtools:
    input:
        bam = os.path.join(MAIN_BAM_DIR, "{sample}.bam")
    output:
        spp = os.path.join(QC_ROOT, "enrichment", "cross_correlation", "{sample}.spp.out"),
        pdf = os.path.join(QC_ROOT, "enrichment", "cross_correlation", "{sample}.cross_correlation.pdf")
    params:
        job_name = f"{config['GSE']}_{{sample}}_PHANTOMPEAK"
    resources:
        time_min = int(custom_param("phantompeakqualtools_time", 60)),
        cpus = int(custom_param("phantompeakqualtools_cpus", 2)),
        mem_gb = int(custom_param("phantompeakqualtools_mem_gb", 16))
    log:
        os.path.join(LOG_ROOT, "Quality_Control", "phantompeakqualtools.{sample}.log")
    benchmark:
        os.path.join(LOG_ROOT, "Quality_Control", "benchmark", "phantompeakqualtools.{sample}.benchmark")
    conda:
        os.path.join(ENV_DIR, "phantompeakqualtools.yaml")
    shell:
        """
        mkdir -p {QC_ROOT}/enrichment/cross_correlation \
                 {LOG_ROOT}/Quality_Control/benchmark

        run_spp.R \
          -c={input.bam} \
          -p={resources.cpus} \
          -savp={output.pdf} \
          -out={output.spp} \
          &> {log}
        """

# The fraction of reads falling within blacklisted regions is used as an
# indicator of potential technical artifacts.
rule blacklist_fraction:
    input:
        bam = os.path.join(MAIN_BAM_DIR, "{sample}.bam"),
        bai = os.path.join(MAIN_BAM_DIR, "{sample}.bam.bai")
    output:
        tsv = os.path.join(QC_ROOT, "blacklist", "{sample}.blacklist_fraction.tsv")
    params:
        job_name = f"{config['GSE']}_{{sample}}_BLACKLIST",
        blacklist = config[ORGANISM]["blacklist_bed"]
    resources:
        time_min = int(custom_param("blacklist_fraction_time", 20)),
        cpus = int(custom_param("blacklist_fraction_cpus", 2)),
        mem_gb = int(custom_param("blacklist_fraction_mem_gb", 16))
    log:
        os.path.join(LOG_ROOT, "Quality_Control", "blacklist_fraction.{sample}.log")
    benchmark:
        os.path.join(LOG_ROOT, "Quality_Control", "benchmark", "blacklist_fraction.{sample}.benchmark")
    conda:
        os.path.join(ENV_DIR, "samtools_picard.yaml")
    shell:
        """
        mkdir -p {QC_ROOT}/blacklist \
                 {LOG_ROOT}/Quality_Control/benchmark

        total=$(samtools view -@ {resources.cpus} -c -F 4 {input.bam} 2> {log})
        blacklisted=$(samtools view -@ {resources.cpus} -c -F 4 -L {params.blacklist} {input.bam} 2>> {log})

        awk -v s="{wildcards.sample}" -v t="$total" -v b="$blacklisted" \
        'BEGIN {{
            frac = (t > 0) ? b/t : 0;
            print "sample\\ttotal_mapped_reads\\tblacklisted_reads\\tblacklist_fraction";
            print s"\\t"t"\\t"b"\\t"frac;
        }}' > {output.tsv}
        """

# Fingerprint plots provide a genome-wide assessment of signal enrichment across
# samples.
rule plot_fingerprint:
    input:
        bams = expand(os.path.join(MAIN_BAM_DIR, "{sample}.bam"), sample=SAMPLES),
        bais = expand(os.path.join(MAIN_BAM_DIR, "{sample}.bam.bai"), sample=SAMPLES)
    output:
        pdf = os.path.join(QC_ROOT, "enrichment", "fingerprint", f"{config['GSE']}.plotFingerprint.pdf"),
        raw = os.path.join(QC_ROOT, "enrichment", "fingerprint", f"{config['GSE']}.plotFingerprint.raw_counts.tsv")
    params:
        job_name = f"{config['GSE']}_PLOT_FINGERPRINT",
        blacklist = config[ORGANISM]["blacklist_bed"]
    resources:
        time_min = int(custom_param("plot_fingerprint_time", 120)),
        cpus = int(custom_param("plot_fingerprint_cpus", 8)),
        mem_gb = int(custom_param("plot_fingerprint_mem_gb", 32))
    log:
        os.path.join(LOG_ROOT, "Quality_Control", f"plotFingerprint.{config['GSE']}.log")
    benchmark:
        os.path.join(LOG_ROOT, "Quality_Control", "benchmark", f"plotFingerprint.{config['GSE']}.benchmark")
    conda:
        os.path.join(ENV_DIR, "deeptools.yaml")
    shell:
        """
        mkdir -p {QC_ROOT}/enrichment/fingerprint \
                 {LOG_ROOT}/Quality_Control/benchmark

        plotFingerprint \
          -b {input.bams} \
          --smartLabels \
          --skipZeros \
          --numberOfSamples 500000 \
          --binSize 500 \
          --blackListFileName {params.blacklist} \
          -p {resources.cpus} \
          --plotFile {output.pdf} \
          --outRawCounts {output.raw} \
          &> {log}
        """


# ------------------------------------------------------------------------------
# 4. BAM-level reproducibility QC
# ------------------------------------------------------------------------------
# Correlation and PCA analyses are performed on genome-wide read-count matrices
# computed from the final BAM files. This provides an overview of sample similarity
# before coverage-track normalization.

rule plot_correlation_multibam:
    input:
        npz = os.path.join(OUT_ROOT, "multiBamSummary", f"{config['GSE']}_multiBamSummary.npz")
    output:
        heatmap = os.path.join(QC_ROOT, "reproducibility", "correlation", f"{config['GSE']}.multiBamSummary.spearman_heatmap.pdf"),
        matrix = os.path.join(QC_ROOT, "reproducibility", "correlation", f"{config['GSE']}.multiBamSummary.spearman_matrix.tab")
    params:
        job_name = f"{config['GSE']}_PLOT_CORR_MULTIBAM"
    resources:
        time_min = int(custom_param("plot_correlation_multibam_time", 30)),
        cpus = int(custom_param("plot_correlation_multibam_cpus", 1)),
        mem_gb = int(custom_param("plot_correlation_multibam_mem_gb", 16))
    log:
        os.path.join(LOG_ROOT, "Quality_Control", f"plotCorrelation.multiBamSummary.{config['GSE']}.log")
    conda:
        os.path.join(ENV_DIR, "deeptools.yaml")
    shell:
        """
        mkdir -p {QC_ROOT}/reproducibility/correlation

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
        npz = os.path.join(OUT_ROOT, "multiBamSummary", f"{config['GSE']}_multiBamSummary.npz")
    output:
        pdf = os.path.join(QC_ROOT, "reproducibility", "PCA", f"{config['GSE']}.multiBamSummary.PCA.pdf"),
        table = os.path.join(QC_ROOT, "reproducibility", "PCA", f"{config['GSE']}.multiBamSummary.PCA.tab")
    params:
        job_name = f"{config['GSE']}_PLOT_PCA_MULTIBAM"
    resources:
        time_min = int(custom_param("plot_pca_multibam_time", 30)),
        cpus =  int(custom_param("plot_pca_multibam_cpus", 1)),
        mem_gb = int(custom_param("plot_pca_multibam_mem_gb", 16))
    log:
        os.path.join(LOG_ROOT, "Quality_Control", f"plotPCA.multiBamSummary.{config['GSE']}.log")
    conda:
        os.path.join(ENV_DIR, "deeptools.yaml")
    shell:
        """
        mkdir -p {QC_ROOT}/reproducibility/PCA

        plotPCA \
          -in {input.npz} \
          -o {output.pdf} \
          --outFileNameData {output.table} \
          &> {log}
        """

# ------------------------------------------------------------------------------
# 5. Normalized signal reproducibility QC
# ------------------------------------------------------------------------------
# Correlation and PCA analyses are also performed on normalized bigWig signal
# matrices. This evaluates sample similarity after applying either spike-in or CPM
# normalization, depending on the experimental design.

if spikeIN == "yes":

    # For spike-in datasets, reproducibility is assessed for both median-based and
    # DESeq2-based normalized signal tracks.
    rule plot_correlation_multibigwig:
        input:
            median = os.path.join(OUT_ROOT, "multiBigwigSummary", f"{config['GSE']}_multiBigwigSummary.spikein_median.npz"),
            deseq2 = os.path.join(OUT_ROOT, "multiBigwigSummary", f"{config['GSE']}_multiBigwigSummary.spikein_deseq2.npz")
        output:
            median_heatmap = os.path.join(QC_ROOT, "reproducibility", "normalized_signal", "spikein_median", "correlation", f"{config['GSE']}.spikein_median.spearman_heatmap.pdf"),
            median_matrix  = os.path.join(QC_ROOT, "reproducibility", "normalized_signal", "spikein_median", "correlation", f"{config['GSE']}.spikein_median.spearman_matrix.tab"),
            deseq2_heatmap = os.path.join(QC_ROOT, "reproducibility", "normalized_signal", "spikein_deseq2", "correlation", f"{config['GSE']}.spikein_deseq2.spearman_heatmap.pdf"),
            deseq2_matrix  = os.path.join(QC_ROOT, "reproducibility", "normalized_signal", "spikein_deseq2", "correlation", f"{config['GSE']}.spikein_deseq2.spearman_matrix.tab")
        params:
            job_name = f"{config['GSE']}_PLOT_CORR_MULTIBW"
        resources:
            time_min = int(custom_param("plot_correlation_multibw_time", 30)),
            cpus = int(custom_param("plot_correlation_multibw_mem_cpus", 1)),
            mem_gb = int(custom_param("plot_correlation_multibw_mem_gb", 16))
        log:
            os.path.join(LOG_ROOT, "Quality_Control", f"plotCorrelation.multiBigwigSummary.{config['GSE']}.log")
        conda:
            os.path.join(ENV_DIR, "deeptools.yaml")
        shell:
            """
            mkdir -p {QC_ROOT}/reproducibility/normalized_signal/spikein_median/correlation \
                     {QC_ROOT}/reproducibility/normalized_signal/spikein_deseq2/correlation \
                     {LOG_ROOT}/Quality_Control

            plotCorrelation -in {input.median} \
              --corMethod spearman --skipZeros --whatToPlot heatmap --plotNumbers \
              -o {output.median_heatmap} \
              --outFileCorMatrix {output.median_matrix} \
              &> {log}

            plotCorrelation -in {input.deseq2} \
              --corMethod spearman --skipZeros --whatToPlot heatmap --plotNumbers \
              -o {output.deseq2_heatmap} \
              --outFileCorMatrix {output.deseq2_matrix} \
              &>> {log}
            """
    rule plot_pca_multibigwig:
        input:
            median = os.path.join(OUT_ROOT, "multiBigwigSummary", f"{config['GSE']}_multiBigwigSummary.spikein_median.npz"),
            deseq2 = os.path.join(OUT_ROOT, "multiBigwigSummary", f"{config['GSE']}_multiBigwigSummary.spikein_deseq2.npz")
        output:
            median_pdf   = os.path.join(QC_ROOT, "reproducibility", "normalized_signal", "spikein_median", "PCA", f"{config['GSE']}.spikein_median.PCA.pdf"),
            median_table = os.path.join(QC_ROOT, "reproducibility", "normalized_signal", "spikein_median", "PCA", f"{config['GSE']}.spikein_median.PCA.tab"),
            deseq2_pdf   = os.path.join(QC_ROOT, "reproducibility", "normalized_signal", "spikein_deseq2", "PCA", f"{config['GSE']}.spikein_deseq2.PCA.pdf"),
            deseq2_table = os.path.join(QC_ROOT, "reproducibility", "normalized_signal", "spikein_deseq2", "PCA", f"{config['GSE']}.spikein_deseq2.PCA.tab")
        params:
            job_name = f"{config['GSE']}_PLOT_PCA_MULTIBW"
        resources:
            time_min = int(custom_param("plot_pca_multibw_time", 30)),
            cpus = int(custom_param("plot_pca_multibw_mem_cpus", 1)),
            mem_gb = int(custom_param("plot_pca_multibw_mem_gb", 16))
        log:
            os.path.join(LOG_ROOT, "Quality_Control", f"plotPCA.multiBigwigSummary.{config['GSE']}.log")
        conda:
            os.path.join(ENV_DIR, "deeptools.yaml")
        shell:
            """
            mkdir -p {QC_ROOT}/reproducibility/normalized_signal/spikein_median/PCA \
                     {QC_ROOT}/reproducibility/normalized_signal/spikein_deseq2/PCA \
                     {LOG_ROOT}/Quality_Control

            plotPCA -in {input.median} \
              -o {output.median_pdf} \
              --outFileNameData {output.median_table} \
              &> {log}

            plotPCA -in {input.deseq2} \
              -o {output.deseq2_pdf} \
              --outFileNameData {output.deseq2_table} \
              &>> {log}
            """
else:

    # For datasets without spike-in, reproducibility is assessed using CPM-normalized
    # signal tracks.
    rule plot_correlation_multibigwig:
        input:
            npz = os.path.join(OUT_ROOT, "multiBigwigSummary", f"{config['GSE']}_multiBigwigSummary.npz")
        output:
            heatmap = os.path.join(QC_ROOT, "reproducibility", "normalized_signal", "CPM", "correlation", f"{config['GSE']}.CPM.spearman_heatmap.pdf"),
            matrix  = os.path.join(QC_ROOT, "reproducibility", "normalized_signal", "CPM", "correlation", f"{config['GSE']}.CPM.spearman_matrix.tab")
        params:
            job_name = f"{config['GSE']}_PLOT_CORR_MULTIBW"
        resources:
            time_min = int(custom_param("plot_correlation_multibw_time", 30)),
            cpus = int(custom_param("plot_correlation_multibw_mem_cpus", 1)),
            mem_gb = int(custom_param("plot_correlation_multibw_mem_gb", 16))
        log:
            os.path.join(LOG_ROOT, "Quality_Control", f"plotCorrelation.multiBigwigSummary.CPM.{config['GSE']}.log")
        conda:
            os.path.join(ENV_DIR, "deeptools.yaml")
        shell:
            """
            mkdir -p {QC_ROOT}/reproducibility/normalized_signal/CPM/correlation \
                     {LOG_ROOT}/Quality_Control

            plotCorrelation -in {input.npz} \
              --corMethod spearman --skipZeros --whatToPlot heatmap --plotNumbers \
              -o {output.heatmap} \
              --outFileCorMatrix {output.matrix} \
              &> {log}
            """
    rule plot_pca_multibigwig:
        input:
            npz = os.path.join(OUT_ROOT, "multiBigwigSummary", f"{config['GSE']}_multiBigwigSummary.npz")
        output:
            pdf   = os.path.join(QC_ROOT, "reproducibility", "normalized_signal", "CPM", "PCA", f"{config['GSE']}.CPM.PCA.pdf"),
            table = os.path.join(QC_ROOT, "reproducibility", "normalized_signal", "CPM", "PCA", f"{config['GSE']}.CPM.PCA.tab")
        params:
            job_name = f"{config['GSE']}_PLOT_PCA_MULTIBW"
        resources:
            time_min = int(custom_param("plot_pca_multibw_time", 30)),
            cpus = int(custom_param("plot_pca_multibw_mem_cpus", 1)),
            mem_gb = int(custom_param("plot_pca_multibw_mem_gb", 16))
        log:
            os.path.join(LOG_ROOT, "Quality_Control", f"plotPCA.multiBigwigSummary.CPM.{config['GSE']}.log")
        conda:
            os.path.join(ENV_DIR, "deeptools.yaml")
        shell:
            """
            mkdir -p {QC_ROOT}/reproducibility/normalized_signal/CPM/PCA \
                     {LOG_ROOT}/Quality_Control

            plotPCA -in {input.npz} \
              -o {output.pdf} \
              --outFileNameData {output.table} \
              &> {log}
            """

# ------------------------------------------------------------------------------
# 6. Integrated MultiQC report
# ------------------------------------------------------------------------------
# MultiQC aggregates QC metrics, logs and summary files generated across the whole
# workflow into a single interactive report.

rule multiqc:
    input:
        MULTIQC_INPUTS
    output:
        html = os.path.join(QC_ROOT, "multiqc", f"{config['GSE']}_multiqc_report.html"),
        data = directory(os.path.join(QC_ROOT, "multiqc", f"{config['GSE']}_multiqc_report_data"))
    params:
        job_name = f"{config['GSE']}_MULTIQC",
        search_dirs = f"{QC_ROOT} {LOG_ROOT}",
        outdir = os.path.join(QC_ROOT, "multiqc"),
        name = f"{config['GSE']}_multiqc_report"
    resources:
        time_min = int(custom_param("multiqc_time", 30)),
        cpus = int(custom_param("multiqc_cpus", 2)),
        mem_gb = int(custom_param("multiqc_mem_gb", 16))
    log:
        os.path.join(LOG_ROOT, "Quality_Control", f"multiqc.{config['GSE']}.log")
    benchmark:
        os.path.join(LOG_ROOT, "Quality_Control", "benchmark", f"multiqc.{config['GSE']}.benchmark")
    conda:
        os.path.join(ENV_DIR, "multiqc.yaml")
    shell:
        """
        mkdir -p {params.outdir} \
                 {LOG_ROOT}/Quality_Control/benchmark

        multiqc {params.search_dirs} \
          --outdir {params.outdir} \
          --filename {params.name}.html \
          --force \
          &> {log}
        """