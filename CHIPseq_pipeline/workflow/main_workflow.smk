# ==============================================================================
# Main ChIP-seq processing workflow
# ==============================================================================
# This module performs the core processing steps required to generate analysis-ready
# BAM files and genome-wide coverage tracks from raw sequencing reads.
#
# Main steps:
#   1. Adapter and quality trimming
#   2. Read alignment
#   3. BAM sorting, duplicate removal and indexing
#   4. Spike-in/main-genome BAM splitting, if spike-in normalization is enabled
#   5. Spike-in scaling and bigWig generation
#   6. Summary matrices for downstream QC and exploratory analysis
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
# 2. Read alignment
# ------------------------------------------------------------------------------
# Reads are aligned either with STAR or Bowtie2, depending on the configuration.
# When spike-in normalization is enabled, reads are mapped against a combined
# genome containing both the main organism and the spike-in reference.


if ALIGNER == "star":

 # STAR is run in end-to-end mode with spliced alignment disabled, which is more
 # appropriate for ChIP-seq than the default RNA-seq-oriented configuration.  
    rule star_align:
        input:
            rules.trimgalore_paired.output if LAYOUT == 'PAIRED' else rules.trimgalore_single.output
        output:
            bam = os.path.join(OUT_ROOT, ALIGN_DIR, "{sample}.Aligned.out.bam")
        params:
            job_name = f"{config['GSE']}_{{sample}}_STAR",
            genome_dir = config[ORGANISM]["star_spike_in_genome" if spikeIN == "yes" else "star_genome"],
            prefix = os.path.join(OUT_ROOT, ALIGN_DIR, "{sample}."),
            paired_extra = "--alignEndsProtrude 10 ConcordantPair --alignMatesGapMax 200" if LAYOUT == "PAIRED" else ""
        resources:
            time_min = int(custom_param('star_time', 180)),
            cpus = int(custom_param('star_cpus', 16)),
            mem_gb = int(custom_param("star_mem_gb", 128))
        log:
            os.path.join(LOG_ROOT, ALIGN_DIR, "STAR_mapping.{sample}.log")
        benchmark:
            os.path.join(LOG_ROOT, ALIGN_DIR, "benchmark", "STAR_mapping.{sample}.benchmark")
        conda:
            os.path.join(ENV_DIR, "STAR_samtools.yaml")
        shell:
            """
            mkdir -p {OUT_ROOT}/{ALIGN_DIR} {LOG_ROOT}/{ALIGN_DIR}/benchmark

            STAR --runThreadN {resources.cpus} \
            --runMode alignReads \
            --alignEndsType EndToEnd \
            --alignIntronMax 1 \
            --readFilesIn {input} \
            --readFilesCommand zcat \
            --outFileNamePrefix {params.prefix} \
            --outFilterMatchNmin 3 \
            --outFilterMismatchNoverLmax 0.1 \
            --outSAMtype BAM Unsorted \
            --outMultimapperOrder Random \
            --outFilterMultimapNmax 1 \
            --outSAMmultNmax 1 \
            --genomeDir {params.genome_dir} \
            {params.paired_extra} \
            &> {log}
            """
elif ALIGNER == "bowtie2":

 # Bowtie2 provides an alternative aligner for standard ChIP-seq read mapping.
    rule bowtie2_align:
        input:
            rules.trimgalore_paired.output if LAYOUT == "PAIRED" else rules.trimgalore_single.output
        output:
            bam = os.path.join(OUT_ROOT, ALIGN_DIR, "{sample}.bowtie2.unsorted.bam")
        params:
            job_name = f"{config['GSE']}_{{sample}}_BOWTIE2",
            index = config[ORGANISM]["bowtie2_spike_in_index" if spikeIN == "yes" else "bowtie2_index"],
            reads = lambda wildcards, input: (
                f"-1 {input.r1} -2 {input.r2}" if LAYOUT == "PAIRED"
                else f"-U {input.fq}")
        resources:
            time_min = int(custom_param("bowtie2_time", 180)),
            cpus = int(custom_param("bowtie2_cpus", 16)),
            mem_gb = int(custom_param("bowtie2_mem_gb", 32))
        log:
            os.path.join(LOG_ROOT, ALIGN_DIR, "BOWTIE2_mapping.{sample}.log")
        benchmark:
            os.path.join(LOG_ROOT, ALIGN_DIR, "benchmark", "BOWTIE2_mapping.{sample}.benchmark")
        conda:
            os.path.join(ENV_DIR, "bowtie2_samtools.yaml")
        shell:
            """
            set -euo pipefail

            mkdir -p {OUT_ROOT}/{ALIGN_DIR}
            mkdir -p {LOG_ROOT}/{ALIGN_DIR}/benchmark

            bowtie2 \
              -x {params.index} \
              {params.reads} \
              --very-sensitive \
              --no-unal \
              -p {resources.cpus} \
              2> {log} \
            | samtools view -bS - > {output.bam}
            """

# ------------------------------------------------------------------------------
# 3. BAM processing
# ------------------------------------------------------------------------------
# Aligned reads are converted into final analysis-ready BAM files.
# This step adds read-group information, removes PCR duplicates, sorts the BAM by
# genomic coordinates and creates the corresponding BAM index.

rule index:
    input:
        bam = rules.star_align.output.bam if ALIGNER == "star" else rules.bowtie2_align.output.bam
    output:
        bam_grouped = temp(os.path.join(OUT_ROOT, ALIGN_DIR, "{sample}_grouped.bam")),
        bam_sorted  = temp(os.path.join(OUT_ROOT, ALIGN_DIR, "{sample}_sorted.bam")),
        marked_bam  = temp(os.path.join(OUT_ROOT, ALIGN_DIR, "{sample}_marked.bam")),
        dup_metrics = os.path.join(QC_ROOT, "duplicates", "picard_markduplicates", "{sample}.marked_dup_metrics.txt"),
        bam         = os.path.join(OUT_ROOT, ALIGN_DIR, "{sample}.bam"),
        bai         = os.path.join(OUT_ROOT, ALIGN_DIR, "{sample}.bam.bai")
    params:
        job_name = f"{config['GSE']}_{{sample}}_INDEX",
        prefix = os.path.join(OUT_ROOT, ALIGN_DIR, "{sample}")
    resources:
        time_min = int(custom_param("index_time", 30)),
        cpus = int(custom_param("index_cpus", 1)),
        mem_gb = int(custom_param("index_mem_gb", 16))
    log:
        os.path.join(LOG_ROOT, ALIGN_DIR, "index.{sample}.log")
    benchmark:
        os.path.join(LOG_ROOT, ALIGN_DIR, "benchmark", "index.{sample}.benchmark")
    conda:
        os.path.join(ENV_DIR, "samtools_picard.yaml")
    shell: 
        """
        mkdir -p {OUT_ROOT}/{ALIGN_DIR} {LOG_ROOT}/{ALIGN_DIR}/benchmark {QC_ROOT}/duplicates/picard_markduplicates

        samtools addreplacerg \
          -r "@RG\tID:{wildcards.sample}\tSM:{wildcards.sample}\tPL:ILLUMINA\tLB:{wildcards.sample}"\
          -o {output.bam_grouped} {input.bam} \
          &>> {log}

        picard SortSam \
          -I {output.bam_grouped} -O {output.bam_sorted} -SO queryname \
          &>> {log}

        picard MarkDuplicates \
          -I {output.bam_sorted} -O {output.marked_bam} \
          -M {output.dup_metrics} \
          --ASSUME_SORT_ORDER queryname --REMOVE_DUPLICATES true \
          &>> {log}

        samtools sort -@ {resources.cpus} -o {output.bam} {output.marked_bam} \
          &>> {log}

        samtools index {output.bam} \
          &>> {log}
        """

# ------------------------------------------------------------------------------
# 4. Split combined BAM into main-genome and spike-in BAMs
# ------------------------------------------------------------------------------
# When spike-in normalization is enabled, reads are first aligned to a combined
# reference genome. This rule separates the final BAM into:
#   - main-organism reads, used for downstream ChIP-seq analyses
#   - spike-in reads, used only to estimate normalization factors
#
# The main-organism prefix is removed from chromosome names to keep compatibility
# with standard genome annotations, blacklist files and downstream tools.
# The spike-in prefix is retained to clearly distinguish spike-in contigs.


if spikeIN == "yes":
    rule split_bam:
        input:
            bam = os.path.join(OUT_ROOT, ALIGN_DIR, "{sample}.bam"),
            bai = os.path.join(OUT_ROOT, ALIGN_DIR, "{sample}.bam.bai")
        output:
            main_bam = os.path.join(OUT_ROOT, "split_bam/main", "{sample}.bam"),
            main_bai = os.path.join(OUT_ROOT, "split_bam/main", "{sample}.bam.bai"),
            spikein_bam = os.path.join(OUT_ROOT, "split_bam/spikein", "{sample}.spikein.bam"),
            spikein_bai = os.path.join(OUT_ROOT, "split_bam/spikein", "{sample}.spikein.bam.bai"),
            combined_idxstats = os.path.join(OUT_ROOT, "split_bam/idxstats", "{sample}.combined.idxstats.txt"),
            main_idxstats_prerename = os.path.join(OUT_ROOT, "split_bam/idxstats", "{sample}.main.prerename.idxstats.txt"),
            main_idxstats_postrename = os.path.join(OUT_ROOT, "split_bam/idxstats", "{sample}.main.postrename.idxstats.txt"),
            spikein_idxstats = os.path.join(OUT_ROOT, "split_bam/idxstats", "{sample}.spikein.idxstats.txt")
        params:
            job_name = f"{config['GSE']}_{{sample}}_SPLITBAM",
            main_prefix = config[ORGANISM]["main_prefix"],
            spikein_prefix = config[ORGANISM]["spikein_prefix"]
        resources:
            time_min = int(custom_param("split_bam_time", 30)),
            cpus = int(custom_param("split_bam_cpus", 1)),
            mem_gb = int(custom_param("split_bam_mem_gb", 16))
        log:
            os.path.join(LOG_ROOT, "split_bam", "split_bam.{sample}.log")
        conda:
            os.path.join(ENV_DIR, "samtools_picard.yaml")
        shell:
            r"""
            set -euo pipefail

            mkdir -p {OUT_ROOT}/split_bam/main \
                     {OUT_ROOT}/split_bam/spikein \
                     {OUT_ROOT}/split_bam/idxstats \
                     {LOG_ROOT}/split_bam

            : > {log}

            samtools idxstats {input.bam} > {output.combined_idxstats} 2>> {log}

            main_refs=$(awk -v p="{params.main_prefix}" 'index($1,p)==1 {{print $1}}' {output.combined_idxstats} | paste -sd " " -)
            spikein_refs=$(awk -v p="{params.spikein_prefix}" 'index($1,p)==1 {{print $1}}' {output.combined_idxstats} | paste -sd " " -)

            [ -n "$main_refs" ] || (echo "ERROR: no main-organism refs found with prefix {params.main_prefix}" >> {log}; exit 1)
            [ -n "$spikein_refs" ] || (echo "ERROR: no spike-in refs found with prefix {params.spikein_prefix}" >> {log}; exit 1)

            echo "Main refs: $main_refs" >> {log}
            echo "Spike-in refs: $spikein_refs" >> {log}

            # Main BAM before renaming:
            # - keep only reads mapped to main-organism references
            # - keep only main-organism @SQ header lines
            samtools view -h {input.bam} $main_refs 2>> {log} \
                | awk -v mp="{params.main_prefix}" 'BEGIN {{OFS="\t"}}
                    /^@SQ/ {{
                        if ($0 ~ "SN:" mp) print
                        next
                    }}
                    /^@/ {{
                        print
                        next
                    }}
                    {{
                        if ($7 != "=" && $7 != "*" && $7 !~ "^" mp) {{
                            $7 = "*"
                            $8 = 0
                            $9 = 0
                        }}
                        print
                    }}' \
                | samtools view -b -o {output.main_bam}.tmp.bam - 2>> {log}

            samtools index {output.main_bam}.tmp.bam 2>> {log}
            samtools idxstats {output.main_bam}.tmp.bam > {output.main_idxstats_prerename} 2>> {log}

            # Main BAM after renaming:
            # - remove main_prefix from header and alignment records
            #   e.g. hs_chr1 -> chr1
            samtools view -h {output.main_bam}.tmp.bam 2>> {log} \
                | awk -v mp="{params.main_prefix}" 'BEGIN {{OFS="\t"}}
                    /^@SQ/ {{
                        sub("SN:" mp, "SN:")
                        print
                        next
                    }}
                    /^@/ {{
                        print
                        next
                    }}
                    {{
                        sub("^" mp, "", $3)
                        if ($7 != "=" && $7 != "*") sub("^" mp, "", $7)
                        print
                    }}' \
                | samtools view -b -o {output.main_bam} - 2>> {log}

            samtools index {output.main_bam} 2>> {log}
            samtools idxstats {output.main_bam} > {output.main_idxstats_postrename} 2>> {log}

            # Spike-in BAM:
            # - keep only reads mapped to spike-in references
            # - keep only spike-in @SQ header lines
            # - keep spike-in prefix
            samtools view -h {input.bam} $spikein_refs 2>> {log} \
                | awk -v sp="{params.spikein_prefix}" 'BEGIN {{OFS="\t"}}
                    /^@SQ/ {{
                        if ($0 ~ "SN:" sp) print
                        next
                    }}
                    /^@/ {{
                        print
                        next
                    }}
                    {{
                        if ($7 != "=" && $7 != "*" && $7 !~ "^" sp) {{
                            $7 = "*"
                            $8 = 0
                            $9 = 0
                        }}
                        print
                    }}' \
                | samtools view -b -o {output.spikein_bam} - 2>> {log}

            samtools index {output.spikein_bam} 2>> {log}
            samtools idxstats {output.spikein_bam} > {output.spikein_idxstats} 2>> {log}

            # Sanity checks to ensure that genome splitting and chromosome renaming worked as expected.
            if grep -q '^{params.spikein_prefix}' {output.main_idxstats_postrename}; then
                echo "ERROR: spike-in references found in main BAM" >> {log}
                exit 1
            fi

            if grep -q '^{params.main_prefix}' {output.main_idxstats_postrename}; then
                echo "ERROR: main prefix still present in main BAM after renaming" >> {log}
                exit 1
            fi

            if grep -q '^{params.main_prefix}' {output.spikein_idxstats}; then
                echo "ERROR: main-organism references found in spike-in BAM" >> {log}
                exit 1
            fi

            main_pre_count=$(awk '$1!="*" {{sum+=$3}} END {{print sum+0}}' {output.main_idxstats_prerename})
            main_post_count=$(awk '$1!="*" {{sum+=$3}} END {{print sum+0}}' {output.main_idxstats_postrename})

            if [ "$main_pre_count" != "$main_post_count" ]; then
                echo "ERROR: main mapped read count changed after renaming" >> {log}
                echo "Before rename: $main_pre_count" >> {log}
                echo "After rename: $main_post_count" >> {log}
                exit 1
            fi

            samtools quickcheck -v {output.main_bam} 2>> {log}
            samtools quickcheck -v {output.spikein_bam} 2>> {log}

            rm -f {output.main_bam}.tmp.bam {output.main_bam}.tmp.bam.bai

            echo "Done" >> {log}
            """

# ------------------------------------------------------------------------------
# 5. Genome-wide BAM summary matrix
# ------------------------------------------------------------------------------
# multiBamSummary computes read-count matrices across genomic bins from the final
# BAM files. These matrices are used later for sample correlation and PCA analyses.


rule multiBamSummary:
    input:
        expand(os.path.join(MAIN_BAM_DIR, "{sample}.bam"), sample=SAMPLES)
    output:
        npz = os.path.join(OUT_ROOT, "multiBamSummary", f"{config['GSE']}_multiBamSummary.npz"),
        raw_counts = os.path.join(OUT_ROOT, "multiBamSummary", f"{config['GSE']}_multiBamSummary.txt")
    params:
        job_name = f"{config['GSE']}_multiBamSummary"
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
        mkdir -p {OUT_ROOT}/multiBamSummary {LOG_ROOT}/multiBamSummary/benchmark

        multiBamSummary bins --bamfiles {input} -o {output.npz} \
          --smartLabels --binSize 100000 -p {resources.cpus} \
          --verbose --outRawCounts {output.raw_counts} --extendReads 75 \
          &> {log}
        """

# ------------------------------------------------------------------------------
# 6. Spike-in read counting
# ------------------------------------------------------------------------------
# Spike-in mapped reads are counted per sample from the spike-in BAM files.
# These counts are used to estimate sample-specific scaling factors.


if spikeIN == "yes":
    rule spikein_counts:
        input:
            bams = expand(os.path.join(OUT_ROOT, "split_bam/spikein", "{sample}.spikein.bam"), sample=SAMPLES),
            bais = expand(os.path.join(OUT_ROOT, "split_bam/spikein", "{sample}.spikein.bam.bai"), sample=SAMPLES)
        output:
            os.path.join(OUT_ROOT, "spikein_counts", "spike_reads_table.tsv")
        params:
            job_name = f"{config['GSE']}_SPIKECOUNT_TABLE",
            spikein_prefix = config[ORGANISM]["spikein_prefix"]
        resources:
            time_min = int(custom_param("spikein_counts_time", 10)),
            cpus = int(custom_param("spikein_counts_cpus", 1)),
            mem_gb = int(custom_param("spikein_counts_mem_gb", 16))
        log:
            os.path.join(LOG_ROOT, "spikein_counts", "spikein_counts_table.log")
        benchmark:
            os.path.join(LOG_ROOT, "spikein_counts", "benchmark", "spikein_counts_table.benchmark")
        conda:
            os.path.join(ENV_DIR, "samtools_picard.yaml")
        shell:
            r"""
            set -euo pipefail

            mkdir -p {OUT_ROOT}/spikein_counts {LOG_ROOT}/spikein_counts/benchmark

            : > {log}
            echo -e "sample\tspike_reads" > {output}

            for BAM in {input.bams}; do
                SAMPLE=$(basename "$BAM" .spikein.bam)

                echo "Counting spike-in reads for sample: $SAMPLE" >> {log}
                echo "BAM: $BAM" >> {log}
                echo "Only contigs starting with {params.spikein_prefix} will be counted" >> {log}

                SPIKE_READS=$(samtools idxstats "$BAM" \
                    | awk -v p="{params.spikein_prefix}" 'index($1,p)==1 {{sum += $3}} END {{print sum+0}}')

                echo -e "${{SAMPLE}}\t${{SPIKE_READS}}" >> {output}
                echo "Spike-in mapped reads: $SPIKE_READS" >> {log}
                echo "" >> {log}
            done
            """

# ------------------------------------------------------------------------------
# 7. Spike-in scaling factors
# ------------------------------------------------------------------------------
# Sample-specific scaling factors are computed from spike-in read counts.
# The resulting table includes both median-based and DESeq2-based normalization
# factors, which are used to generate spike-in-normalized coverage tracks.

if spikeIN == "yes":
    rule scale_factors:
        input:
            counts = os.path.join(OUT_ROOT, "spikein_counts", "spike_reads_table.tsv")
        output:
            os.path.join(OUT_ROOT, "spikein_counts", "scale_factors.tsv")
        params:
            script = os.path.join(SCRIPT_DIR, "spikein_scale_factors.R"),
            job_name = f"{config['GSE']}_SCALEFACTORS"
        resources:
            time_min = int(custom_param("scale_factors_time", 10)),
            cpus = int(custom_param("scale_factors_cpus", 1)),
            mem_gb = int(custom_param("scale_factors_mem_gb", 16))
        log:
            os.path.join(LOG_ROOT, "spikein_counts", f"scale_factors.{config['GSE']}.log")
        benchmark:
            os.path.join(LOG_ROOT, "spikein_counts", "benchmark", f"scale_factors.{config['GSE']}.benchmark")
        conda:
            os.path.join(ENV_DIR, "deseq2.yaml")
        shell:
            r"""
            mkdir -p {OUT_ROOT}/spikein_counts \
                     {LOG_ROOT}/spikein_counts/benchmark

            Rscript {params.script} \
                {input.counts} \
                {output} \
                &> {log}
            """

# ------------------------------------------------------------------------------
# 8. Coverage track generation
# ------------------------------------------------------------------------------
# Genome-wide bigWig coverage tracks are generated with deepTools bamCoverage.
# If spike-in normalization is enabled, external scaling factors are applied to
# the main-organism BAM files. Otherwise, CPM normalization is used.

if spikeIN == "yes":

    # Generate bigWig tracks using the median-based spike-in scaling factor.
    rule bamcoverage_bigwig_spikein_median:
        input:
            bam = os.path.join(OUT_ROOT, "split_bam/main", "{sample}.bam"),
            bai = os.path.join(OUT_ROOT, "split_bam/main", "{sample}.bam.bai"),
            scales = os.path.join(OUT_ROOT, "spikein_counts", "scale_factors.tsv")
        output:
            os.path.join(OUT_ROOT, "bigwig_spikein_median", "{sample}.coverage.spikein_median.bw")
        params:
            job_name = f"{config['GSE']}_{{sample}}_BAMCOV_SPIKEIN_MEDIAN",
            blacklist = config[ORGANISM]["blacklist_bed"],
            deeptools_genome_size = config[ORGANISM]["deeptools_genome_size"],
            ignoreForNormalization = config[ORGANISM]["ignoreForNormalization"]
        resources:
            time_min = int(custom_param("bigwig_time", 60)),
            cpus = int(custom_param("bigwig_cpus", 20)),
            mem_gb = int(custom_param("bigwig_mem_gb", 32))
        conda:
            os.path.join(ENV_DIR, "deeptools.yaml")
        log:
            os.path.join(LOG_ROOT, "bigwig_spikein_median", "bamCoverage.{sample}.spikein_median.log")
        benchmark:
            os.path.join(LOG_ROOT, "bigwig_spikein_median", "benchmark", "bamCoverage.{sample}.spikein_median.benchmark")
        shell:
            r"""
            set -euo pipefail

            mkdir -p {OUT_ROOT}/bigwig_spikein_median \
                     {LOG_ROOT}/bigwig_spikein_median/benchmark

            SCALE=$(awk -F '\t' -v s="{wildcards.sample}" 'NR>1 && $1==s {{print $3}}' {input.scales})

            if [ -z "$SCALE" ]; then
                echo "ERROR: scale_factor_median not found for sample {wildcards.sample}" > {log}
                exit 1
            fi

            echo "Sample: {wildcards.sample}" > {log}
            echo "Using scale_factor_median: $SCALE" >> {log}

            bamCoverage \
                -b {input.bam} \
                -o {output} \
                -p {resources.cpus} \
                --binSize 1 \
                --normalizeUsing RPGC \
                --effectiveGenomeSize {params.deeptools_genome_size} \
                --ignoreForNormalization {params.ignoreForNormalization} \
                --scaleFactor "$SCALE" \
                --blackListFileName {params.blacklist} \
                --skipNAs \
                &>> {log}
            """

    # Generate bigWig tracks using the DESeq2-based spike-in scaling factor.
    rule bamcoverage_bigwig_spikein_deseq2:
        input:
            bam = os.path.join(OUT_ROOT, "split_bam/main", "{sample}.bam"),
            bai = os.path.join(OUT_ROOT, "split_bam/main", "{sample}.bam.bai"),
            scales = os.path.join(OUT_ROOT, "spikein_counts", "scale_factors.tsv")
        output:
            os.path.join(OUT_ROOT, "bigwig_spikein_deseq2", "{sample}.coverage.spikein_deseq2.bw")
        params:
            job_name = f"{config['GSE']}_{{sample}}_BAMCOV_SPIKEIN_DESEQ2",
            blacklist = config[ORGANISM]["blacklist_bed"],
            deeptools_genome_size = config[ORGANISM]["deeptools_genome_size"],
            ignoreForNormalization = config[ORGANISM]["ignoreForNormalization"]
        resources:
            time_min = int(custom_param("bigwig_time", 60)),
            cpus = int(custom_param("bigwig_cpus", 20)),
            mem_gb = int(custom_param("bigwig_mem_gb", 32))
        conda:
            os.path.join(ENV_DIR, "deeptools.yaml")
        log:
            os.path.join(LOG_ROOT, "bigwig_spikein_deseq2", "bamCoverage.{sample}.spikein_deseq2.log")
        benchmark:
            os.path.join(LOG_ROOT, "bigwig_spikein_deseq2", "benchmark", "bamCoverage.{sample}.spikein_deseq2.benchmark")
        shell:
            r"""
            set -euo pipefail

            mkdir -p {OUT_ROOT}/bigwig_spikein_deseq2 \
                     {LOG_ROOT}/bigwig_spikein_deseq2/benchmark

            SCALE=$(awk -F '\t' -v s="{wildcards.sample}" 'NR>1 && $1==s {{print $5}}' {input.scales})

            if [ -z "$SCALE" ]; then
                echo "ERROR: scale_factor_deseq2 not found for sample {wildcards.sample}" > {log}
                exit 1
            fi

            echo "Sample: {wildcards.sample}" > {log}
            echo "Using scale_factor_deseq2: $SCALE" >> {log}

            bamCoverage \
                -b {input.bam} \
                -o {output} \
                -p {resources.cpus} \
                --binSize 1 \
                --normalizeUsing RPGC \
                --effectiveGenomeSize {params.deeptools_genome_size} \
                --ignoreForNormalization {params.ignoreForNormalization} \
                --scaleFactor "$SCALE" \
                --blackListFileName {params.blacklist} \
                --skipNAs \
                &>> {log}
            """
else:
    # For datasets without spike-in, coverage tracks are normalized by CPM.
    rule bamcoverage_bigwig:
        input:
            rules.index.output.bam
        output:
            os.path.join(OUT_ROOT, "bigwig", "{sample}.coverage.bw")
        params:
            job_name = f"{config['GSE']}_{{sample}}_BAMCOV",
            blacklist = config[ORGANISM]["blacklist_bed"],
            genomeSize = config[ORGANISM]["deeptools_genome_size"],
            ignoreForNormalization = config[ORGANISM]["ignoreForNormalization"]
        resources:
            time_min = int(custom_param("bigwig_time", 60)),
            cpus = int(custom_param("bigwig_cpus", 20)),
            mem_gb = int(custom_param("bigwig_mem_gb", 32))
        conda:
            os.path.join(ENV_DIR, "deeptools.yaml")
        log:
            os.path.join(LOG_ROOT, "bigwig", "bamCoverage.{sample}.coverage.log")
        benchmark:
            os.path.join(LOG_ROOT, "bigwig", "benchmark", "bamCoverage.{sample}.coverage.benchmark")
        shell:
            """
            mkdir -p {OUT_ROOT}/bigwig {LOG_ROOT}/bigwig/benchmark

            bamCoverage -b {input} -p {resources.cpus} -o {output} \
              --binSize 1 \
              --normalizeUsing CPM \
              --effectiveGenomeSize {params.genomeSize} \
              --ignoreForNormalization {params.ignoreForNormalization} \
              --blackListFileName {params.blacklist} \
              --skipNAs \
              &> {log}
            """

# ------------------------------------------------------------------------------
# 9. Genome-wide bigWig summary matrix
# ------------------------------------------------------------------------------
# multiBigwigSummary summarizes normalized coverage tracks across genomic bins.
# The resulting matrices are used for downstream correlation and PCA plots.

if spikeIN == "yes":
    rule multiBigwigSummary:
        input:
            median = expand(os.path.join(OUT_ROOT,"bigwig_spikein_median","{sample}.coverage.spikein_median.bw"),sample=SAMPLES),
            deseq2 = expand(os.path.join(OUT_ROOT,"bigwig_spikein_deseq2","{sample}.coverage.spikein_deseq2.bw"),sample=SAMPLES)
        output:
            median_npz = os.path.join(OUT_ROOT,"multiBigwigSummary",f"{config['GSE']}_multiBigwigSummary.spikein_median.npz"),
            median_txt = os.path.join(OUT_ROOT,"multiBigwigSummary",f"{config['GSE']}_multiBigwigSummary.spikein_median.txt"),
            deseq2_npz = os.path.join(OUT_ROOT,"multiBigwigSummary",f"{config['GSE']}_multiBigwigSummary.spikein_deseq2.npz"),
            deseq2_txt = os.path.join(OUT_ROOT,"multiBigwigSummary",f"{config['GSE']}_multiBigwigSummary.spikein_deseq2.txt")
        params:
            job_name = f"{config['GSE']}_multiBigwigSummary"
        resources:
            time_min = int(custom_param("multiBigwigSummary_time", 120)),
            cpus = int(custom_param("multiBigwigSummary_cpus", 20)),
            mem_gb = int(custom_param("multiBigwigSummary_mem_gb", 32))
        log:
            os.path.join(LOG_ROOT, "multiBigwigSummary", f"multiBigwigSummary.{config['GSE']}.log")
        benchmark:
            os.path.join(LOG_ROOT, "multiBigwigSummary", "benchmark", f"multiBigwigSummary.{config['GSE']}.benchmark")
        conda:
            os.path.join(ENV_DIR, "deeptools.yaml")
        shell:
            """
            set -euo pipefail

            mkdir -p {OUT_ROOT}/multiBigwigSummary \
                     {LOG_ROOT}/multiBigwigSummary/benchmark

            : > {log}

            echo "Running multiBigwigSummary with spike-in median normalized bigWigs" >> {log}

            multiBigwigSummary bins \
              --bwfiles {input.median} \
              -o {output.median_npz} \
              --smartLabels \
              --binSize 100000 \
              -p {resources.cpus} \
              --outRawCounts {output.median_txt} \
              &>> {log}

            echo "" >> {log}
            echo "Running multiBigwigSummary with spike-in DESeq2 normalized bigWigs" >> {log}

            multiBigwigSummary bins \
              --bwfiles {input.deseq2} \
              -o {output.deseq2_npz} \
              --smartLabels \
              --binSize 100000 \
              -p {resources.cpus} \
              --outRawCounts {output.deseq2_txt} \
              &>> {log}
            """

else:
    rule multiBigwigSummary:
        input:
            expand(os.path.join(OUT_ROOT, "bigwig", "{sample}.coverage.bw"), sample=SAMPLES)
        output:
            npz = os.path.join(
                OUT_ROOT,
                "multiBigwigSummary",
                f"{config['GSE']}_multiBigwigSummary.npz"
            ),
            txt = os.path.join(
                OUT_ROOT,
                "multiBigwigSummary",
                f"{config['GSE']}_multiBigwigSummary.txt"
            )
        params:
            job_name = f"{config['GSE']}_multiBigwigSummary"
        resources:
            time_min = int(custom_param("multiBigwigSummary_time", 60)),
            cpus = int(custom_param("multiBigwigSummary_cpus", 20)),
            mem_gb = int(custom_param("multiBigwigSummary_mem_gb", 32))
        log:
            os.path.join(LOG_ROOT, "multiBigwigSummary", f"multiBigwigSummary.{config['GSE']}.log")
        benchmark:
            os.path.join(LOG_ROOT, "multiBigwigSummary", "benchmark", f"multiBigwigSummary.{config['GSE']}.benchmark")
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
              --binSize 100000 \
              -p {resources.cpus} \
              --outRawCounts {output.txt} \
              &> {log}
            """