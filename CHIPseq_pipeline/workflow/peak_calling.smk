# ==============================================================================
# Peak calling and peak-level analysis
# ==============================================================================
# This module runs MACS3 peak calling and derives consensus and union peak sets for
# downstream quantification and quality control.
#
# Peak calling is controlled by the metadata table loaded in the main Snakefile,
# which defines IP samples, input controls, experimental conditions, ChIP targets
# and peak modes.
#
# Main steps:
#   1. Per-sample MACS3 peak calling
#   2. Consensus peak generation per condition and mark
#   3. Union peak set generation per mark
#   4. Read counting in union peaks
#   5. Peak-level QC and FRiP metrics
#   6. Signal heatmaps over union peaks
# ==============================================================================


if calling_peaks == "yes":

    # --------------------------------------------------------------------------
    # 1. Per-sample MACS3 peak calling
    # --------------------------------------------------------------------------
    # MACS3 is run independently for each IP sample using the matching input
    # control when available. Narrow and/or broad peak calling is selected from
    # the metadata table.

    rule calling_peaks:
        input:
            treatment = os.path.join(MAIN_BAM_DIR, "{sample}.bam"),
            control = get_control_bam
        output:
            peaks = os.path.join(
                OUT_ROOT, "peaks", "macs3", "{peak_mode}", "{sample}",
                "{sample}_peaks.{peak_mode}Peak"
            ),
            xls = os.path.join(
                OUT_ROOT, "peaks", "macs3", "{peak_mode}", "{sample}",
                "{sample}_peaks.xls"
            )
        params:
            job_name = f"{config['GSE']}_{{sample}}_MACS3_{{peak_mode}}",
            outdir = os.path.join(OUT_ROOT, "peaks", "macs3", "{peak_mode}", "{sample}"),
            genome_size = config[ORGANISM]["macs3_genome_size"],
            fmt = MACS_FORMAT,
            control_arg = lambda wildcards, input: f"-c {input.control}" if input.control else "",
            peak_mode_args = lambda wc: (
                f"--broad --broad-cutoff {config['macs3']['broad_cutoff']}"
                if wc.peak_mode == "broad"
                else f"-q {config['macs3']['narrow_cutoff']}"
            )
        wildcard_constraints:
            peak_mode = "narrow|broad"
        resources:
            time_min = int(custom_param("macs3_time", 120)),
            cpus = int(custom_param("macs3_cpus", 2)),
            mem_gb = int(custom_param("macs3_mem_gb", 16))
        log:
            os.path.join(LOG_ROOT, "peaks", "macs3", "{peak_mode}", "calling_peaks.{sample}.log")
        benchmark:
            os.path.join(
                LOG_ROOT, "peaks", "macs3", "benchmark",
                "calling_peaks.{peak_mode}.{sample}.benchmark"
            )
        conda:
            os.path.join(ENV_DIR, "macs3.yaml")
        shell:
            """
            mkdir -p {params.outdir} \
                     {LOG_ROOT}/peaks/macs3/{wildcards.peak_mode} \
                     {LOG_ROOT}/peaks/macs3/benchmark

            macs3 callpeak \
              -t {input.treatment} \
              {params.control_arg} \
              -f {params.fmt} \
              -g {params.genome_size} \
              -n {wildcards.sample} \
              --outdir {params.outdir} \
              {params.peak_mode_args} \
              &> {log}
            """

    # --------------------------------------------------------------------------
    # 2. Consensus peak generation
    # --------------------------------------------------------------------------
    # Replicate peak sets from the same condition, mark and peak mode are combined
    # to retain regions supported by a minimum number of biological replicates.
    # The replicate threshold is configurable and automatically capped by the
    # number of available peak files.

    rule consensus_peaks:
        input:
            peaks = get_peaks_for_consensus
        output:
            consensus  = os.path.join(OUT_ROOT, "peaks", "consensus", "{peak_mode}", "{mark}", "{condition}.{mark}.{peak_mode}.consensus.bed"),
            multiinter = os.path.join(OUT_ROOT, "peaks", "consensus", "{peak_mode}", "{mark}", "{condition}.{mark}.{peak_mode}.multiinter.bed"),
            summary    = os.path.join(QC_ROOT, "peaks", "consensus", "{condition}.{mark}.{peak_mode}.consensus.summary.tsv")
        params:
            job_name = f"{config['GSE']}_{{condition}}_{{mark}}_{{peak_mode}}_CONSENSUS",
            min_replicates = lambda wildcards, input: min(
                int(config.get("consensus_peaks", {}).get("min_replicates", 2)),
                len(input.peaks)),
            n_replicates   = lambda wildcards, input: len(input.peaks),
            merge_distance = 0
        wildcard_constraints:
            peak_mode = "narrow|broad"
        resources:
            time_min = int(custom_param("consensus_peaks_time", 30)),
            cpus = int(custom_param("consensus_peaks_cpus", 1)),
            mem_gb = int(custom_param("consensus_peaks_mem_gb", 16))
        log:
            os.path.join(LOG_ROOT, "peaks", "consensus", "{peak_mode}", "consensus_peaks.{condition}.{mark}.log")
        benchmark:
            os.path.join(LOG_ROOT, "peaks", "consensus", "benchmark", "consensus_peaks.{condition}.{mark}.{peak_mode}.benchmark")
        conda:
            os.path.join(ENV_DIR, "bedtools.yaml")
        shell:
            r"""
            set -euo pipefail

            mkdir -p {OUT_ROOT}/peaks/consensus/{wildcards.peak_mode}/{wildcards.mark} \
                     {QC_ROOT}/peaks/consensus \
                     {LOG_ROOT}/peaks/consensus/{wildcards.peak_mode} \
                     {LOG_ROOT}/peaks/consensus/benchmark

            : > {log}

            echo "Consensus peak generation" >> {log}
            echo "Condition: {wildcards.condition}" >> {log}
            echo "Mark: {wildcards.mark}" >> {log}
            echo "Peak mode: {wildcards.peak_mode}" >> {log}
            echo "Number of replicates: {params.n_replicates}" >> {log}
            echo "Minimum replicate support required: {params.min_replicates}" >> {log}
            echo "Merge distance: {params.merge_distance}" >> {log}
            echo "Input peak files:" >> {log}
            printf '%s\n' {input.peaks} >> {log}

            TMPDIR=$(mktemp -d)
            trap "rm -rf $TMPDIR" EXIT

            i=0
            for BED in {input.peaks}; do
                i=$((i + 1))

                awk 'BEGIN {{OFS="\t"}} !/^#/ {{print $1, $2, $3}}' "$BED" \
                    | sort -k1,1 -k2,2n \
                    > "$TMPDIR/rep${{i}}.bed"

                echo "Prepared sorted BED: $TMPDIR/rep${{i}}.bed" >> {log}
            done

            bedtools multiinter -i "$TMPDIR"/rep*.bed \
                > {output.multiinter} 2>> {log}

            awk -v min_rep="{params.min_replicates}" 'BEGIN {{OFS="\t"}} $4 >= min_rep {{print $1, $2, $3}}' {output.multiinter} \
                | sort -k1,1 -k2,2n \
                > "$TMPDIR/consensus.filtered.bed"

            if [ -s "$TMPDIR/consensus.filtered.bed" ]; then
                bedtools merge -i "$TMPDIR/consensus.filtered.bed" -d {params.merge_distance} \
                    > {output.consensus} 2>> {log}
            else
                : > {output.consensus}
                echo "WARNING: no consensus peaks passed the replicate support threshold." >> {log}
            fi

            N_CONSENSUS=$(wc -l < {output.consensus})

            {
                echo -e "condition\tmark\tpeak_mode\tn_replicates\tmin_replicates\tmerge_distance\tconsensus_peaks\tinput_peak_files"
                echo -e "{wildcards.condition}\t{wildcards.mark}\t{wildcards.peak_mode}\t{params.n_replicates}\t{params.min_replicates}\t{params.merge_distance}\t${{N_CONSENSUS}}\t{input.peaks}"
            } > {output.summary}

            echo "Consensus peaks generated: $N_CONSENSUS" >> {log}
            echo "Done" >> {log}
            """

    # --------------------------------------------------------------------------
    # 3. Union peak set generation
    # --------------------------------------------------------------------------
    # Consensus peaks from all conditions are merged for each mark and peak mode.
    # This creates a common peak reference used to quantify signal across samples.

    rule union_consensus_peaks:
        input:
            consensus = get_consensus_peaks_for_union
        output:
            union = os.path.join(OUT_ROOT, "peaks", "union", "{peak_mode}", "{mark}", "{mark}.{peak_mode}.union_peaks.bed"),
            summary = os.path.join(QC_ROOT, "peaks", "union", "{mark}.{peak_mode}.union_peaks.summary.tsv")
        params:
            job_name = f"{config['GSE']}_{{mark}}_{{peak_mode}}_UNION_PEAKS",
            merge_distance = 0,
            n_consensus_files = lambda wildcards, input: len(input.consensus)
        wildcard_constraints:
            peak_mode = "narrow|broad"
        resources:
            time_min = int(custom_param("union_peaks_time", 30)),
            cpus = int(custom_param("union_peaks_cpus", 1)),
            mem_gb = int(custom_param("union_peaks_mem_gb", 16))
        log:
            os.path.join(LOG_ROOT, "peaks", "union", "{peak_mode}", "union_consensus_peaks.{mark}.log")
        benchmark:
            os.path.join(LOG_ROOT, "peaks", "union", "benchmark", "union_consensus_peaks.{mark}.{peak_mode}.benchmark")
        conda:
            os.path.join(ENV_DIR, "bedtools.yaml")
        shell:
            r"""
            set -euo pipefail

            mkdir -p {OUT_ROOT}/peaks/union/{wildcards.peak_mode}/{wildcards.mark} \
                     {QC_ROOT}/peaks/union \
                     {LOG_ROOT}/peaks/union/{wildcards.peak_mode} \
                     {LOG_ROOT}/peaks/union/benchmark

            : > {log}

            echo "Union peak generation" >> {log}
            echo "Mark: {wildcards.mark}" >> {log}
            echo "Peak mode: {wildcards.peak_mode}" >> {log}
            echo "Number of consensus files: {params.n_consensus_files}" >> {log}
            echo "Merge distance: {params.merge_distance}" >> {log}
            echo "Input consensus peak files:" >> {log}
            printf '%s\n' {input.consensus} >> {log}

            TMPDIR=$(mktemp -d)
            trap "rm -rf $TMPDIR" EXIT

            cat {input.consensus} \
                | awk 'BEGIN {{OFS="\t"}} !/^#/ {{print $1, $2, $3}}' \
                | sort -k1,1 -k2,2n \
                > "$TMPDIR/all_consensus_peaks.sorted.bed"

            if [ -s "$TMPDIR/all_consensus_peaks.sorted.bed" ]; then
                bedtools merge \
                    -i "$TMPDIR/all_consensus_peaks.sorted.bed" \
                    -d {params.merge_distance} \
                    > {output.union} 2>> {log}
            else
                : > {output.union}
                echo "WARNING: no consensus peaks found for union." >> {log}
            fi

            N_UNION=$(wc -l < {output.union})
            TOTAL_BP=$(awk '{{sum += $3-$2}} END {{print sum+0}}' {output.union})

            {
                echo -e "mark\tpeak_mode\tn_consensus_files\tmerge_distance\tunion_peaks\ttotal_union_bp\tinput_consensus_files"
                echo -e "{wildcards.mark}\t{wildcards.peak_mode}\t{params.n_consensus_files}\t{params.merge_distance}\t${{N_UNION}}\t${{TOTAL_BP}}\t{input.consensus}"
            } > {output.summary}

            echo "Union peaks generated: $N_UNION" >> {log}
            echo "Total union peak bp: $TOTAL_BP" >> {log}
            echo "Done" >> {log}
            """

    # --------------------------------------------------------------------------
    # 4. Read counting in union peaks
    # --------------------------------------------------------------------------
    # Union peak BED files are converted to SAF format and used with featureCounts
    # to generate peak-by-sample count matrices for downstream quantitative analyses.

    rule count_reads_in_union_peaks:
        input:
            peaks = os.path.join(OUT_ROOT, "peaks", "union", "{peak_mode}", "{mark}", "{mark}.{peak_mode}.union_peaks.bed"),
            bams = lambda wc: [
                os.path.join(MAIN_BAM_DIR, f"{sample}.bam")
                for sample in get_samples_for_union_counts(wc)
            ],
            bais = lambda wc: [
                os.path.join(MAIN_BAM_DIR, f"{sample}.bam.bai")
                for sample in get_samples_for_union_counts(wc)
            ]
        output:
            saf = os.path.join(OUT_ROOT, "peaks", "counts", "{peak_mode}", "{mark}", "{mark}.{peak_mode}.union_peaks.saf"),
            raw = os.path.join(OUT_ROOT, "peaks", "counts", "{peak_mode}", "{mark}", "{mark}.{peak_mode}.featureCounts.raw.tsv"),
            featurecounts_summary = os.path.join(OUT_ROOT, "peaks", "counts", "{peak_mode}", "{mark}", "{mark}.{peak_mode}.featureCounts.raw.tsv.summary"),
            matrix = os.path.join(OUT_ROOT, "peaks", "counts", "{peak_mode}", "{mark}", "{mark}.{peak_mode}.counts_matrix.tsv"),
            summary = os.path.join(QC_ROOT, "peaks", "counts", "{mark}.{peak_mode}.counts.summary.tsv")
        params:
            job_name = f"{config['GSE']}_{{mark}}_{{peak_mode}}_COUNTS_UNION",
            samples = lambda wc: "\t".join(get_samples_for_union_counts(wc)),
            n_samples = lambda wc: len(get_samples_for_union_counts(wc)),
            paired_args = "-p --countReadPairs -B -C" if LAYOUT == "PAIRED" else ""
        wildcard_constraints:
            peak_mode = "narrow|broad"
        resources:
            time_min = int(custom_param("count_matrix_time", 60)),
            cpus = int(custom_param("count_matrix_cpus", 4)),
            mem_gb = int(custom_param("count_matrix_mem_gb", 16))
        log:
            os.path.join(LOG_ROOT, "peaks", "counts", "{peak_mode}", "count_reads_in_union_peaks.{mark}.log")
        benchmark:
            os.path.join(LOG_ROOT, "peaks", "counts", "benchmark", "count_reads_in_union_peaks.{mark}.{peak_mode}.benchmark")
        conda:
            os.path.join(ENV_DIR, "subread.yaml")
        shell:
            r"""
            set -euo pipefail

            mkdir -p {OUT_ROOT}/peaks/counts/{wildcards.peak_mode}/{wildcards.mark} \
                     {QC_ROOT}/peaks/counts \
                     {LOG_ROOT}/peaks/counts/{wildcards.peak_mode} \
                     {LOG_ROOT}/peaks/counts/benchmark

            : > {log}

            echo "Counting fragments/reads in union peaks" >> {log}
            echo "Mark: {wildcards.mark}" >> {log}
            echo "Peak mode: {wildcards.peak_mode}" >> {log}
            echo "Union peak file: {input.peaks}" >> {log}
            echo "Number of samples: {params.n_samples}" >> {log}
            echo "Samples:" >> {log}
            echo -e "{params.samples}" | tr '\t' '\n' >> {log}
            echo "BAM files:" >> {log}
            printf '%s\n' {input.bams} >> {log}

            awk 'BEGIN {{OFS="\t"; print "GeneID","Chr","Start","End","Strand"}}
                 {{
                    peak_id = "peak_" NR;
                    print peak_id, $1, $2 + 1, $3, ".";
                 }}' {input.peaks} > {output.saf}

            featureCounts \
              -T {resources.cpus} \
              -F SAF \
              -a {output.saf} \
              -o {output.raw} \
              {params.paired_args} \
              {input.bams} \
              &>> {log}

            {
                echo -e "peak_id\tchr\tstart\tend\t{params.samples}"
                awk 'BEGIN {{OFS="\t"}}
                     $1 !~ /^#/ && $1 != "Geneid" {{
                        printf "%s\t%s\t%s\t%s", $1, $2, $3 - 1, $4;
                        for (i=7; i<=NF; i++) printf "\t%s", $i;
                        printf "\n";
                     }}' {output.raw}
            } > {output.matrix}

            N_PEAKS=$(wc -l < {input.peaks})
            N_ROWS=$(($(wc -l < {output.matrix}) - 1))

            {
                echo -e "mark\tpeak_mode\tn_samples\tn_union_peaks\tn_count_rows\tcounting_tool\tpaired_mode\tcount_matrix"
                echo -e "{wildcards.mark}\t{wildcards.peak_mode}\t{params.n_samples}\t${{N_PEAKS}}\t${{N_ROWS}}\tfeatureCounts\t{LAYOUT}\t{output.matrix}"
            } > {output.summary}

            echo "Count matrix generated: {output.matrix}" >> {log}
            echo "Number of union peaks: $N_PEAKS" >> {log}
            echo "Number of count rows: $N_ROWS" >> {log}
            echo "Done" >> {log}
            """

    # --------------------------------------------------------------------------
    # 5. Peak-level QC and FRiP metrics
    # --------------------------------------------------------------------------
    # Peak-level QC summarizes the number, genomic span and width distribution of
    # individual and consensus peak sets. FRiP metrics quantify the fraction of
    # mapped reads falling within called peak regions.

    # Summarize per-sample MACS3 peak sets.
    rule individual_peak_qc:
        input:
            peaks = os.path.join(OUT_ROOT, "peaks", "macs3", "{peak_mode}", "{sample}", "{sample}_peaks.{peak_mode}Peak")
        output:
            tsv = os.path.join(QC_ROOT, "peaks", "individual", "{sample}.{peak_mode}.peak_qc.tsv")
        params:
            job_name = f"{config['GSE']}_{{sample}}_{{peak_mode}}_PEAK_QC"
        wildcard_constraints:
            peak_mode = "narrow|broad"
        resources:
            time_min = int(custom_param("peak_qc_time", 20)),
            cpus = int(custom_param("peak_qc_cpus", 1)),
            mem_gb = int(custom_param("peak_qc_mem_gb", 16))
        log:
            os.path.join(LOG_ROOT, "peaks", "qc", "individual_peak_qc.{sample}.{peak_mode}.log")
        benchmark:
            os.path.join(LOG_ROOT, "peaks", "qc", "benchmark", "individual_peak_qc.{sample}.{peak_mode}.benchmark")
        conda:
            os.path.join(ENV_DIR, "bedtools.yaml")
        shell:
            r"""
            set -euo pipefail

            mkdir -p {QC_ROOT}/peaks/individual {LOG_ROOT}/peaks/qc/benchmark

            awk -v sample="{wildcards.sample}" -v mode="{wildcards.peak_mode}" '
            BEGIN {{
                OFS="\t"; n=0; total_bp=0; min_width="NA"; max_width=0;
            }}
            !/^#/ {{
                width=$3-$2;
                if (width > 0) {{
                    n++;
                    total_bp += width;
                    if (min_width=="NA" || width < min_width) min_width=width;
                    if (width > max_width) max_width=width;
                }}
            }}
            END {{
                mean_width = (n > 0) ? total_bp/n : 0;
                print "sample","peak_mode","n_peaks","total_peak_bp","mean_peak_width","min_peak_width","max_peak_width";
                print sample,mode,n,total_bp,mean_width,min_width,max_width;
            }}' {input.peaks} > {output.tsv} 2> {log}
            """

    # Summarize condition-level consensus peak sets.
    rule consensus_peak_qc:
        input:
            consensus = os.path.join(OUT_ROOT, "peaks", "consensus", "{peak_mode}", "{mark}", "{condition}.{mark}.{peak_mode}.consensus.bed")
        output:
            tsv = os.path.join(QC_ROOT, "peaks", "consensus", "{condition}.{mark}.{peak_mode}.peak_qc.tsv")
        params:
            job_name = f"{config['GSE']}_{{condition}}_{{mark}}_{{peak_mode}}_CONS_QC"
        wildcard_constraints:
            peak_mode = "narrow|broad"
        resources:
            time_min = int(custom_param("peak_qc_time", 20)),
            cpus = int(custom_param("peak_qc_cpus", 1)),
            mem_gb = int(custom_param("peak_qc_mem_gb", 16))
        log:
            os.path.join(LOG_ROOT, "peaks", "qc", "consensus_peak_qc.{condition}.{mark}.{peak_mode}.log")
        benchmark:
            os.path.join(LOG_ROOT, "peaks", "qc", "benchmark", "consensus_peak_qc.{condition}.{mark}.{peak_mode}.benchmark")
        conda:
            os.path.join(ENV_DIR, "bedtools.yaml")
        shell:
            r"""
            set -euo pipefail

            mkdir -p {QC_ROOT}/peaks/consensus {LOG_ROOT}/peaks/qc/benchmark

            awk -v condition="{wildcards.condition}" -v mark="{wildcards.mark}" -v mode="{wildcards.peak_mode}" '
            BEGIN {{
                OFS="\t"; n=0; total_bp=0; min_width="NA"; max_width=0;
            }}
            !/^#/ {{
                width=$3-$2;
                if (width > 0) {{
                    n++;
                    total_bp += width;
                    if (min_width=="NA" || width < min_width) min_width=width;
                    if (width > max_width) max_width=width;
                }}
            }}
            END {{
                mean_width = (n > 0) ? total_bp/n : 0;
                print "condition","mark","peak_mode","n_consensus_peaks","total_consensus_bp","mean_peak_width","min_peak_width","max_peak_width";
                print condition,mark,mode,n,total_bp,mean_width,min_width,max_width;
            }}' {input.consensus} > {output.tsv} 2> {log}
            """

    # Compute FRiP using each sample's individual MACS3 peak set.
    rule frip_individual_peaks:
        input:
            bam = os.path.join(MAIN_BAM_DIR, "{sample}.bam"),
            bai = os.path.join(MAIN_BAM_DIR, "{sample}.bam.bai"),
            peaks = os.path.join(OUT_ROOT, "peaks", "macs3", "{peak_mode}", "{sample}", "{sample}_peaks.{peak_mode}Peak")
        output:
            tsv = os.path.join(QC_ROOT, "peaks", "frip", "individual", "{sample}.{peak_mode}.frip.tsv")
        params:
            job_name = f"{config['GSE']}_{{sample}}_{{peak_mode}}_FRIP"
        wildcard_constraints:
            peak_mode = "narrow|broad"
        resources:
            time_min = int(custom_param("frip_time", 30)),
            cpus = int(custom_param("frip_cpus", 2)),
            mem_gb = int(custom_param("frip_mem_gb", 16))
        log:
            os.path.join(LOG_ROOT, "peaks", "qc", "frip_individual.{sample}.{peak_mode}.log")
        benchmark:
            os.path.join(LOG_ROOT, "peaks", "qc", "benchmark", "frip_individual.{sample}.{peak_mode}.benchmark")
        conda:
            os.path.join(ENV_DIR, "samtools_picard.yaml")
        shell:
            r"""
            set -euo pipefail

            mkdir -p {QC_ROOT}/peaks/frip/individual {LOG_ROOT}/peaks/qc/benchmark

            total=$(samtools view -@ {resources.cpus} -c -F 4 {input.bam} 2> {log})
            in_peaks=$(samtools view -@ {resources.cpus} -c -F 4 -L {input.peaks} {input.bam} 2>> {log})

            awk -v sample="{wildcards.sample}" -v mode="{wildcards.peak_mode}" -v total="$total" -v in_peaks="$in_peaks" '
            BEGIN {{
                OFS="\t";
                frip = (total > 0) ? in_peaks / total : 0;
                print "sample","peak_mode","peak_set","total_mapped_reads","reads_in_peaks","FRiP";
                print sample,mode,"individual",total,in_peaks,frip;
            }}' > {output.tsv}
            """

    # Compute FRiP using the consensus peak set associated with each sample.
    rule frip_consensus_peaks:
        input:
            bam = os.path.join(MAIN_BAM_DIR, "{sample}.bam"),
            bai = os.path.join(MAIN_BAM_DIR, "{sample}.bam.bai"),
            consensus = get_consensus_for_sample
        output:
            tsv = os.path.join(QC_ROOT, "peaks", "frip", "consensus", "{sample}.{peak_mode}.consensus_frip.tsv")
        params:
            job_name = f"{config['GSE']}_{{condition}}_{{mark}}_{{peak_mode}}_FRIP_CONS"
        wildcard_constraints:
            peak_mode = "narrow|broad"
        resources:
            time_min = int(custom_param("frip_time", 30)),
            cpus = int(custom_param("frip_cpus", 2)),
            mem_gb = int(custom_param("frip_mem_gb", 16))
        log:
            os.path.join(LOG_ROOT, "peaks", "qc", "frip_consensus.{sample}.{peak_mode}.log")
        benchmark:
            os.path.join(LOG_ROOT, "peaks", "qc", "benchmark", "frip_consensus.{sample}.{peak_mode}.benchmark")
        conda:
            os.path.join(ENV_DIR, "samtools_picard.yaml")
        shell:
            r"""
            set -euo pipefail

            mkdir -p {QC_ROOT}/peaks/frip/consensus {LOG_ROOT}/peaks/qc/benchmark

            total=$(samtools view -@ {resources.cpus} -c -F 4 {input.bam} 2> {log})
            in_peaks=$(samtools view -@ {resources.cpus} -c -F 4 -L {input.consensus} {input.bam} 2>> {log})

            awk -v sample="{wildcards.sample}" -v mode="{wildcards.peak_mode}" -v total="$total" -v in_peaks="$in_peaks" '
            BEGIN {{
                OFS="\t";
                frip = (total > 0) ? in_peaks / total : 0;
                print "sample","peak_mode","peak_set","total_mapped_reads","reads_in_consensus_peaks","FRiP";
                print sample,mode,"consensus",total,in_peaks,frip;
            }}' > {output.tsv}
            """

# ------------------------------------------------------------------------------
# 6. Signal heatmaps over union peaks
# ------------------------------------------------------------------------------
# deepTools computeMatrix and plotHeatmap are used to visualize normalized signal
# around the center of union peaks. The input bigWig tracks depend on the selected
# normalization strategy.

if calling_peaks == "yes" and spikeIN == "yes":

    # For spike-in datasets, heatmaps are generated for both median-based and
    # DESeq2-based normalized bigWig tracks.

    rule heatmap_union_peaks_spikein:
        input:
            peaks = os.path.join(OUT_ROOT, "peaks", "union", "{peak_mode}", "{mark}", "{mark}.{peak_mode}.union_peaks.bed"),
            median_bws = lambda wc: [
                os.path.join(OUT_ROOT, "bigwig_spikein_median", f"{sample}.coverage.spikein_median.bw")
                for sample in get_samples_for_union_counts(wc)
            ],
            deseq2_bws = lambda wc: [
                os.path.join(OUT_ROOT, "bigwig_spikein_deseq2", f"{sample}.coverage.spikein_deseq2.bw")
                for sample in get_samples_for_union_counts(wc)
            ]
        output:
            median_matrix = os.path.join(QC_ROOT, "peaks", "signal_heatmap", "spikein_median", "{peak_mode}", "{mark}.{peak_mode}.computeMatrix.gz"),
            median_pdf = os.path.join(QC_ROOT, "peaks", "signal_heatmap", "spikein_median", "{peak_mode}", "{mark}.{peak_mode}.heatmap.pdf"),
            deseq2_matrix = os.path.join(QC_ROOT, "peaks", "signal_heatmap", "spikein_deseq2", "{peak_mode}", "{mark}.{peak_mode}.computeMatrix.gz"),
            deseq2_pdf = os.path.join(QC_ROOT, "peaks", "signal_heatmap", "spikein_deseq2", "{peak_mode}", "{mark}.{peak_mode}.heatmap.pdf")
        params:
            job_name = f"{config['GSE']}_{{mark}}_{{peak_mode}}_HEATMAP_SPIKEIN",
            samples = lambda wc: " ".join(get_samples_for_union_counts(wc))
        wildcard_constraints:
            peak_mode = "narrow|broad"
        resources:
            time_min = int(custom_param("heatmap_union_peaks_time", 120)),
            cpus = int(custom_param("heatmap_union_peaks_cpus", 8)),
            mem_gb = int(custom_param("heatmap_union_peaks_mem_gb", 32))
        log:
            os.path.join(LOG_ROOT, "peaks", "signal_heatmap", "{peak_mode}", "heatmap_union_peaks.{mark}.log")
        benchmark:
            os.path.join(LOG_ROOT, "peaks", "signal_heatmap", "benchmark", "heatmap_union_peaks.{mark}.{peak_mode}.benchmark")
        conda:
            os.path.join(ENV_DIR, "deeptools.yaml")
        shell:
            """
            mkdir -p {QC_ROOT}/peaks/signal_heatmap/spikein_median/{wildcards.peak_mode} \
                     {QC_ROOT}/peaks/signal_heatmap/spikein_deseq2/{wildcards.peak_mode} \
                     {LOG_ROOT}/peaks/signal_heatmap/{wildcards.peak_mode} \
                     {LOG_ROOT}/peaks/signal_heatmap/benchmark

            : > {log}

            echo "Computing heatmap matrix for spike-in median normalized bigWigs" >> {log}

            computeMatrix reference-point \
              --referencePoint center \
              -R {input.peaks} \
              -S {input.median_bws} \
              --samplesLabel {params.samples} \
              -b 3000 \
              -a 3000 \
              --skipZeros \
              --missingDataAsZero \
              -p {resources.cpus} \
              -o {output.median_matrix} \
              &>> {log}

            plotHeatmap \
              -m {output.median_matrix} \
              -out {output.median_pdf} \
              --sortUsing mean \
              --colorMap viridis \
              --whatToShow "heatmap and colorbar" \
              &>> {log}

            echo "Computing heatmap matrix for spike-in DESeq2 normalized bigWigs" >> {log}

            computeMatrix reference-point \
              --referencePoint center \
              -R {input.peaks} \
              -S {input.deseq2_bws} \
              --samplesLabel {params.samples} \
              -b 3000 \
              -a 3000 \
              --skipZeros \
              --missingDataAsZero \
              -p {resources.cpus} \
              -o {output.deseq2_matrix} \
              &>> {log}

            plotHeatmap \
              -m {output.deseq2_matrix} \
              -out {output.deseq2_pdf} \
              --sortUsing mean \
              --colorMap viridis \
              --whatToShow "heatmap and colorbar" \
              &>> {log}
            """

if calling_peaks == "yes" and spikeIN == "no":

    # For datasets without spike-in, heatmaps are generated from CPM-normalized
    # bigWig tracks.
    
    rule heatmap_union_peaks_cpm:
        input:
            peaks = os.path.join(OUT_ROOT, "peaks", "union", "{peak_mode}", "{mark}", "{mark}.{peak_mode}.union_peaks.bed"),
            bws = lambda wc: [
                os.path.join(OUT_ROOT, "bigwig", f"{sample}.coverage.bw")
                for sample in get_samples_for_union_counts(wc)
            ]
        output:
            matrix = os.path.join(QC_ROOT, "peaks", "signal_heatmap", "CPM", "{peak_mode}", "{mark}.{peak_mode}.computeMatrix.gz"),
            pdf = os.path.join(QC_ROOT, "peaks", "signal_heatmap", "CPM", "{peak_mode}", "{mark}.{peak_mode}.heatmap.pdf")
        params:
            job_name = f"{config['GSE']}_{{mark}}_{{peak_mode}}_HEATMAP_CPM",
            samples = lambda wc: " ".join(get_samples_for_union_counts(wc))
        wildcard_constraints:
            peak_mode = "narrow|broad"
        resources:
            time_min = int(custom_param("heatmap_union_peaks_time", 120)),
            cpus = int(custom_param("heatmap_union_peaks_cpus", 8)),
            mem_gb = int(custom_param("heatmap_union_peaks_mem_gb", 32))
        log:
            os.path.join(LOG_ROOT, "peaks", "signal_heatmap", "{peak_mode}", "heatmap_union_peaks.{mark}.log")
        benchmark:
            os.path.join(LOG_ROOT, "peaks", "signal_heatmap", "benchmark", "heatmap_union_peaks.{mark}.{peak_mode}.benchmark")
        conda:
            os.path.join(ENV_DIR, "deeptools.yaml")
        shell:
            """
            mkdir -p {QC_ROOT}/peaks/signal_heatmap/CPM/{wildcards.peak_mode} \
                     {LOG_ROOT}/peaks/signal_heatmap/{wildcards.peak_mode} \
                     {LOG_ROOT}/peaks/signal_heatmap/benchmark

            computeMatrix reference-point \
              --referencePoint center \
              -R {input.peaks} \
              -S {input.bws} \
              --samplesLabel {params.samples} \
              -b 3000 \
              -a 3000 \
              --skipZeros \
              --missingDataAsZero \
              -p {resources.cpus} \
              -o {output.matrix} \
              &> {log}

            plotHeatmap \
              -m {output.matrix} \
              -out {output.pdf} \
              --sortUsing mean \
              --colorMap viridis \
              --whatToShow "heatmap and colorbar" \
              &>> {log}
            """