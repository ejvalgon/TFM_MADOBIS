# ==============================================================================
# Main bulk RNA-seq processing workflow
# ==============================================================================
# This module performs the core processing steps required to generate processed
# BAM files, gene-level count matrices and genome-wide coverage tracks from raw
# sequencing reads.
#
# Main steps:
#   1. Adapter and quality trimming
#   2. Splice-aware read alignment with STAR
#   3. BAM read-group annotation, sorting and indexing
#   4. Gene-level quantification with featureCounts
#   5. bigWig coverage track generation
#   6. Summary matrices for downstream reproducibility QC and exploratory analysis
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Adapter and quality trimming
# ------------------------------------------------------------------------------
# Trim Galore is used to remove sequencing adapters and low-quality bases.
# Separate rules are defined for paired-end and single-end libraries.

if LAYOUT == "PAIRED":

    rule trimgalore_paired:
        input:
            get_raw_reads
        output:
            r1 = os.path.join(OUT_ROOT, "fastq_trimmed", "{sample}_val_1.fq.gz"),
            r2 = os.path.join(OUT_ROOT, "fastq_trimmed", "{sample}_val_2.fq.gz")
        params:
            job_name = f"{config['GSE']}_{{sample}}_TRIMPAIRED"
        resources:
            time_min = int(custom_param("trim_time", 120)),
            cpus = int(custom_param("trim_cpus", 4)),
            mem_gb = int(custom_param("trim_mem_gb", 16))
        log:
            os.path.join(LOG_ROOT, "trimgalore", "Trim_galore.{sample}.log")
        benchmark:
            os.path.join(LOG_ROOT, "trimgalore", "benchmark", "Trim_galore.{sample}.benchmark")
        conda:
            os.path.join(ENV_DIR, "trimgalore_samtools.yaml")
        shell:
            """
            mkdir -p {OUT_ROOT}/fastq_trimmed \
                     {LOG_ROOT}/trimgalore/benchmark

            trim_galore \
              -j {resources.cpus} \
              --stringency 2 \
              -q 30 \
              --paired \
              --basename {wildcards.sample} \
              -o {OUT_ROOT}/fastq_trimmed \
              {input} \
              --gzip \
              &> {log}
            """
else:

    rule trimgalore_single:
        input:
            get_raw_reads
        output:
            fq = os.path.join(OUT_ROOT, "fastq_trimmed", "{sample}_trimmed.fq.gz")
        params:
            job_name = f"{config['GSE']}_{{sample}}_TRIMSINGLE"
        resources:
            time_min = int(custom_param("trim_time", 120)),
            cpus = int(custom_param("trim_cpus", 2)),
            mem_gb = int(custom_param("trim_mem_gb", 16))
        log:
            os.path.join(LOG_ROOT, "trimgalore", "Trim_galore.{sample}.log")
        benchmark:
            os.path.join(LOG_ROOT, "trimgalore", "benchmark", "Trim_galore.{sample}.benchmark")
        conda:
            os.path.join(ENV_DIR, "trimgalore_samtools.yaml")
        shell:
            """
            mkdir -p {OUT_ROOT}/fastq_trimmed \
                     {LOG_ROOT}/trimgalore/benchmark

            trim_galore \
              -j {resources.cpus} \
              --stringency 3 \
              -q 30 \
              --basename {wildcards.sample} \
              -o {OUT_ROOT}/fastq_trimmed \
              {input} \
              --gzip \
              &> {log}
            """


# ------------------------------------------------------------------------------
# 2. Read alignment with STAR
# ------------------------------------------------------------------------------
# Reads are aligned to the reference genome using STAR in RNA-seq mode.
# The output BAM is kept unsorted and will be processed in the next rule.

rule star_align:
    input:
        rules.trimgalore_paired.output if LAYOUT == "PAIRED" else rules.trimgalore_single.output
    output:
        bam = os.path.join(OUT_ROOT, ALIGN_DIR, "{sample}.Aligned.out.bam"),
        log_final = os.path.join(OUT_ROOT, ALIGN_DIR, "{sample}.Log.final.out"),
        sj = os.path.join(OUT_ROOT, ALIGN_DIR, "{sample}.SJ.out.tab"),
        gene_counts = os.path.join(OUT_ROOT, ALIGN_DIR, "{sample}.ReadsPerGene.out.tab")
    params:
        job_name = f"{config['GSE']}_{{sample}}_STAR",
        genome_dir = config[ORGANISM]["star_genome"],
        prefix = os.path.join(OUT_ROOT, ALIGN_DIR, "{sample}."),
        align_intron_max = int(custom_param("star_align_intron_max", 0))
    resources:
        time_min = int(custom_param("star_time", 180)),
        cpus = int(custom_param("star_cpus", 16)),
        mem_gb = int(custom_param("star_mem_gb", 128))
    log:
        os.path.join(LOG_ROOT, ALIGN_DIR, "STAR_mapping.{sample}.log")
    benchmark:
        os.path.join(LOG_ROOT, ALIGN_DIR, "benchmark", "STAR_mapping.{sample}.benchmark")
    conda:
        os.path.join(ENV_DIR, "STAR_samtools.yaml")
    shell:
        """
        mkdir -p {OUT_ROOT}/{ALIGN_DIR} \
                 {LOG_ROOT}/{ALIGN_DIR}/benchmark

        STAR \
          --runThreadN {resources.cpus} \
          --runMode alignReads \
          --genomeDir {params.genome_dir} \
          --readFilesIn {input} \
          --readFilesCommand zcat \
          --outFileNamePrefix {params.prefix} \
          --outSAMtype BAM Unsorted \
          --quantMode GeneCounts \
          --alignIntronMax {params.align_intron_max} \
          &> {log}
        """

# ------------------------------------------------------------------------------
# 3. BAM processing
# ------------------------------------------------------------------------------
# Aligned reads are converted into final analysis-ready BAM files.
# This step adds read-group information, removes PCR duplicates, sorts the BAM by
# genomic coordinates and creates the corresponding BAM index.

rule process_bam:
    input:
        bam = rules.star_align.output.bam
    output:
        bam_grouped = temp(os.path.join(OUT_ROOT, ALIGN_DIR, "{sample}.grouped.bam")),
        bam = os.path.join(OUT_ROOT, ALIGN_DIR, "{sample}.bam"),
        bai = os.path.join(OUT_ROOT, ALIGN_DIR, "{sample}.bam.bai")
    params:
        job_name = f"{config['GSE']}_{{sample}}_PROCESSBAM"
    resources:
        time_min = int(custom_param("process_bam_time", 60)),
        cpus = int(custom_param("process_bam_cpus", 4)),
        mem_gb = int(custom_param("process_bam_mem_gb", 16))
    log:
        os.path.join(LOG_ROOT, ALIGN_DIR, "process_bam.{sample}.log")
    benchmark:
        os.path.join(LOG_ROOT, ALIGN_DIR, "benchmark", "process_bam.{sample}.benchmark")
    conda:
        os.path.join(ENV_DIR, "samtools.yaml")
    shell:
        """
        mkdir -p {OUT_ROOT}/{ALIGN_DIR} \
                 {LOG_ROOT}/{ALIGN_DIR}/benchmark

        samtools addreplacerg \
          -r "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:ILLUMINA\\tLB:{wildcards.sample}" \
          -o {output.bam_grouped} \
          {input.bam} \
          &> {log}

        samtools sort \
          -@ {resources.cpus} \
          -o {output.bam} \
          {output.bam_grouped} \
          &>> {log}

        samtools index \
          -@ {resources.cpus} \
          {output.bam} \
          &>> {log}
        """


# ------------------------------------------------------------------------------
# 4. Gene-level quantification
# ------------------------------------------------------------------------------
# Gene-level raw counts are generated from the final coordinate-sorted BAM files
# using featureCounts and the reference gene annotation.

rule featurecounts:
    input:
        bam = expand(
            os.path.join(OUT_ROOT, ALIGN_DIR, "{sample}.bam"),
            sample = SAMPLES
        )
    output:
        counts = os.path.join(OUT_ROOT, "counts", "featureCounts.txt"),
        summary = os.path.join(OUT_ROOT, "counts", "featureCounts.txt.summary")
    params:
        job_name = f"{config['GSE']}_FEATURECOUNTS",
        gtf = config[ORGANISM]["gtf"],
        paired_flag = "-p --countReadPairs" if LAYOUT == "PAIRED" else "",
        strandedness = int(custom_param("featurecounts_strand", 0))
    resources:
        time_min = int(custom_param("featurecounts_time", 120)),
        cpus = int(custom_param("featurecounts_cpus", 8)),
        mem_gb = int(custom_param("featurecounts_mem_gb", 32))
    log:
        os.path.join(LOG_ROOT, "featurecounts", "featureCounts.log")
    benchmark:
        os.path.join(LOG_ROOT, "featurecounts", "benchmark", "featureCounts.benchmark")
    conda:
        os.path.join(ENV_DIR, "subread.yaml")
    shell:
        """
        mkdir -p {OUT_ROOT}/counts \
                 {LOG_ROOT}/featurecounts/benchmark

        featureCounts \
          -T {resources.cpus} \
          -a {params.gtf} \
          -o {output.counts} \
          -t exon \
          -g gene_id \
          -s {params.strandedness} \
          {params.paired_flag} \
          {input.bam} \
          &> {log}
        """


# ------------------------------------------------------------------------------
# 4. Gene-level quantification with featureCounts
# ------------------------------------------------------------------------------
# Raw count matrices are generated from the final coordinate-sorted BAM files.
#
# Two complementary count matrices are produced:
#   1. Exonic counts: counts reads overlapping annotated exons.
#   2. Gene body counts: counts reads overlapping complete gene regions,
#      including exons and introns.
#
# Exonic counts are closer to conventional mature RNA-seq quantification, whereas
# gene body counts may be more informative for chromatin-associated / nascent RNA.

rule featurecounts_exon:
    input:
        bam = expand(
            os.path.join(OUT_ROOT, ALIGN_DIR, "{sample}.bam"),
            sample = SAMPLES
        )
    output:
        counts = os.path.join(OUT_ROOT, "counts", "exon", "featureCounts_exon.txt"),
        summary = os.path.join(OUT_ROOT, "counts", "exon", "featureCounts_exon.txt.summary")
    params:
        job_name = f"{config['GSE']}_FEATURECOUNTS_EXON",
        gtf = config[ORGANISM]["gtf"],
        paired_flag = "-p --countReadPairs" if LAYOUT == "PAIRED" else "",
        strandedness = int(custom_param("featurecounts_strand", 0))
    resources:
        time_min = int(custom_param("featurecounts_time", 120)),
        cpus = int(custom_param("featurecounts_cpus", 8)),
        mem_gb = int(custom_param("featurecounts_mem_gb", 32))
    log:
        os.path.join(LOG_ROOT, "featurecounts", "featureCounts_exon.log")
    benchmark:
        os.path.join(LOG_ROOT, "featurecounts", "benchmark", "featureCounts_exon.benchmark")
    conda:
        os.path.join(ENV_DIR, "subread.yaml")
    shell:
        """
        mkdir -p {OUT_ROOT}/counts/exon \
                 {LOG_ROOT}/featurecounts/benchmark

        featureCounts \
          -T {resources.cpus} \
          -a {params.gtf} \
          -o {output.counts} \
          -t exon \
          -g gene_id \
          -s {params.strandedness} \
          {params.paired_flag} \
          {input.bam} \
          &> {log}
        """


rule featurecounts_gene_body:
    input:
        bam = expand(
            os.path.join(OUT_ROOT, ALIGN_DIR, "{sample}.bam"),
            sample = SAMPLES
        )
    output:
        counts = os.path.join(OUT_ROOT, "counts", "gene_body", "featureCounts_gene_body.txt"),
        summary = os.path.join(OUT_ROOT, "counts", "gene_body", "featureCounts_gene_body.txt.summary")
    params:
        job_name = f"{config['GSE']}_FEATURECOUNTS_GENEBODY",
        gtf = config[ORGANISM]["gtf"],
        paired_flag = "-p --countReadPairs" if LAYOUT == "PAIRED" else "",
        strandedness = int(custom_param("featurecounts_strand", 0))
    resources:
        time_min = int(custom_param("featurecounts_time", 120)),
        cpus = int(custom_param("featurecounts_cpus", 8)),
        mem_gb = int(custom_param("featurecounts_mem_gb", 32))
    log:
        os.path.join(LOG_ROOT, "featurecounts", "featureCounts_gene_body.log")
    benchmark:
        os.path.join(LOG_ROOT, "featurecounts", "benchmark", "featureCounts_gene_body.benchmark")
    conda:
        os.path.join(ENV_DIR, "subread.yaml")
    shell:
        """
        mkdir -p {OUT_ROOT}/counts/gene_body \
                 {LOG_ROOT}/featurecounts/benchmark

        featureCounts \
          -T {resources.cpus} \
          -a {params.gtf} \
          -o {output.counts} \
          -t gene \
          -g gene_id \
          -s {params.strandedness} \
          {params.paired_flag} \
          {input.bam} \
          &> {log}
        """


# ------------------------------------------------------------------------------
# 5. Genome-wide coverage tracks
# ------------------------------------------------------------------------------
# BigWig files are generated from final coordinate-sorted BAM files using
# deepTools bamCoverage. These tracks are intended for genome browser
# visualization and downstream metagene/profile analyses.

rule bamcoverage_bigwig:
    input:
        bam = rules.process_bam.output.bam,
        bai = rules.process_bam.output.bai
    output:
        bw = os.path.join(OUT_ROOT, "bigwig", "{sample}." + BW_NORMALIZATION + ".bw")
    params:
        job_name = f"{config['GSE']}_{{sample}}_BAMCOVERAGE",
        bin_size = int(custom_param("bamcoverage_bin_size", 1)),
        normalization = BW_NORMALIZATION
    resources:
        time_min = int(custom_param("bamcoverage_time", 120)),
        cpus = int(custom_param("bamcoverage_cpus", 4)),
        mem_gb = int(custom_param("bamcoverage_mem_gb", 16))
    log:
        os.path.join(LOG_ROOT, "bamcoverage", "bamCoverage.{sample}." + BW_NORMALIZATION + ".log")
    benchmark:
        os.path.join(LOG_ROOT, "bamcoverage", "benchmark", "bamCoverage.{sample}." + BW_NORMALIZATION + ".benchmark")
    conda:
        os.path.join(ENV_DIR, "deeptools.yaml")
    shell:
        """
        mkdir -p {OUT_ROOT}/bigwig \
                 {LOG_ROOT}/bamcoverage/benchmark

        bamCoverage \
          -b {input.bam} \
          -o {output.bw} \
          --binSize {params.bin_size} \
          --normalizeUsing {params.normalization} \
          --numberOfProcessors {resources.cpus} \
          &> {log}
        """


# ------------------------------------------------------------------------------
# 5. Genome-wide BAM summary matrix
# ------------------------------------------------------------------------------
# multiBamSummary computes read-count matrices across genomic bins from the final
# BAM files. These matrices are used later for sample correlation and PCA analyses.


rule multiBamSummary:
    input:
        expand(
            os.path.join(OUT_ROOT, ALIGN_DIR, "{sample}.bam"),
            sample = SAMPLES
        )
    output:
        npz = os.path.join(OUT_ROOT, "multiBamSummary", f"{config['GSE']}_multiBamSummary.npz"),
        raw_counts = os.path.join(OUT_ROOT, "multiBamSummary", f"{config['GSE']}_multiBamSummary.txt")
    params:
        job_name = f"{config['GSE']}_multiBamSummary",
        bin_size = int(custom_param("multiBamSummary_binSize", 100000))
    resources:
        time_min = int(custom_param("multiBamSummary_time", 60)),
        cpus = int(custom_param("multiBamSummary_cpus", 20)),
        mem_gb = int(custom_param("multiBamSummary_mem_gb", 32))
    log:
        os.path.join(LOG_ROOT, "multiBamSummary", f"multiBamSummary.{config['GSE']}.log")
    benchmark:
        os.path.join(LOG_ROOT, "multiBamSummary", "benchmark", f"multiBamSummary.{config['GSE']}.benchmark")
    conda:
        os.path.join(ENV_DIR, "deeptools.yaml")
    shell:
        """
        mkdir -p {OUT_ROOT}/multiBamSummary \
                 {LOG_ROOT}/multiBamSummary/benchmark

        multiBamSummary bins \
          --bamfiles {input} \
          -o {output.npz} \
          --smartLabels \
          --binSize {params.bin_size} \
          -p {resources.cpus} \
          --verbose \
          --outRawCounts {output.raw_counts} \
          &> {log}
        """

# ------------------------------------------------------------------------------
# 9. Genome-wide bigWig summary matrix
# ------------------------------------------------------------------------------
# multiBigwigSummary summarizes normalized coverage tracks across genomic bins.
# The resulting matrices are used for downstream correlation and PCA plots.

rule multiBigwigSummary:
    input:
        expand(
            os.path.join(OUT_ROOT, "bigwig", "{sample}." + BW_NORMALIZATION + ".bw"),
            sample = SAMPLES
        )
    output:
        npz = os.path.join(
            OUT_ROOT,
            "multiBigwigSummary",
            f"{config['GSE']}_multiBigwigSummary_{BW_NORMALIZATION}.npz"
        ),
        txt = os.path.join(
            OUT_ROOT,
            "multiBigwigSummary",
            f"{config['GSE']}_multiBigwigSummary_{BW_NORMALIZATION}.txt"
        )
    params:
        job_name = f"{config['GSE']}_multiBigwigSummary_{BW_NORMALIZATION}",
        bin_size = int(custom_param("multiBigwigSummary_binSize", 100000)),
        normalization = BW_NORMALIZATION
    resources:
        time_min = int(custom_param("multiBigwigSummary_time", 60)),
        cpus = int(custom_param("multiBigwigSummary_cpus", 20)),
        mem_gb = int(custom_param("multiBigwigSummary_mem_gb", 32))
    log:
        os.path.join(LOG_ROOT, "multiBigwigSummary", f"multiBigwigSummary.{config['GSE']}.{BW_NORMALIZATION}.log")
    benchmark:
        os.path.join(LOG_ROOT, "multiBigwigSummary", "benchmark", f"multiBigwigSummary.{config['GSE']}.{BW_NORMALIZATION}.benchmark")
    conda:
        os.path.join(ENV_DIR, "deeptools.yaml")
    shell:
        """
        set -euo pipefail

        mkdir -p {OUT_ROOT}/multiBigwigSummary \
                 {LOG_ROOT}/multiBigwigSummary/benchmark

        multiBigwigSummary bins \
          --bwfiles {input} \
          -o {output.npz} \
          --smartLabels \
          --binSize {params.bin_size} \
          -p {resources.cpus} \
          --verbose \
          --outRawCounts {output.txt} \
          &> {log}
        """

