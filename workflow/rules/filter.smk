# =============================================================================
# filter.smk — post-alignment filtering + Tn5 shift
# =============================================================================
# Standard ENCODE-style ATAC filtering:
#   1. keep properly-paired, primary alignments with MAPQ >= threshold
#   2. remove mitochondrial reads (chrM / MT)
#   3. mark & remove PCR/optical duplicates (Picard)
#   4. remove reads overlapping the ENCODE blacklist (Amemiya et al. 2019)
#   5. Tn5 shift: +4 (+ strand) / -5 (- strand) for base-pair-accurate cut sites
# Each step is logged for the MultiQC report.

rule filter_bam:
    input:
        bam=f"{RESULTS}/aligned/{{sample}}.sorted.bam",
    output:
        bam=temp(f"{RESULTS}/filtered/{{sample}}.namefilt.bam"),
    params:
        mapq=config["filtering"]["min_mapq"],
        proper="-f 2" if config["filtering"]["keep_proper_pairs"] else "",
    threads: config["resources"]["sort_threads"]
    log:
        f"{LOGS}/filter/{{sample}}.mapq.log",
    conda:
        "../envs/align.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.bam})
        # -F 1804: exclude unmapped, mate-unmapped, secondary, QC-fail, dup
        samtools view -@ {threads} -b {params.proper} -F 1804 \
            -q {params.mapq} {input.bam} > {output.bam} 2> {log}
        """


rule remove_mito:
    input:
        f"{RESULTS}/filtered/{{sample}}.namefilt.bam",
    output:
        bam=temp(f"{RESULTS}/filtered/{{sample}}.nomito.bam"),
        idx=temp(f"{RESULTS}/filtered/{{sample}}.nomito.bam.bai"),
    params:
        do_remove=config["filtering"]["remove_mito"],
    threads: config["resources"]["sort_threads"]
    log:
        f"{LOGS}/filter/{{sample}}.nomito.log",
    conda:
        "../envs/align.yaml"
    shell:
        r"""
        if [ "{params.do_remove}" = "True" ]; then
            samtools index {input}
            CHRS=$(samtools idxstats {input} | cut -f1 \
                   | grep -v -E '^(chrM|MT|chrMT|Mito)$' | tr '\n' ' ')
            samtools view -@ {threads} -b {input} $CHRS > {output.bam} 2> {log}
        else
            cp {input} {output.bam}
        fi
        samtools index {output.bam}
        """


rule mark_duplicates:
    input:
        f"{RESULTS}/filtered/{{sample}}.nomito.bam",
    output:
        bam=temp(f"{RESULTS}/filtered/{{sample}}.dedup.bam"),
        metrics=f"{RESULTS}/qc/picard/{{sample}}.dup_metrics.txt",
    params:
        do_remove="true" if config["filtering"]["remove_duplicates"] else "false",
    log:
        f"{LOGS}/filter/{{sample}}.markdup.log",
    conda:
        "../envs/align.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.metrics})
        picard MarkDuplicates I={input} O={output.bam} \
            M={output.metrics} REMOVE_DUPLICATES={params.do_remove} \
            VALIDATION_STRINGENCY=LENIENT 2> {log}
        """


rule remove_blacklist:
    input:
        bam=f"{RESULTS}/filtered/{{sample}}.dedup.bam",
        blacklist=blacklist_bed() if blacklist_bed() else [],
    output:
        bam=f"{RESULTS}/filtered/{{sample}}.filtered.bam",
        bai=f"{RESULTS}/filtered/{{sample}}.filtered.bam.bai",
    params:
        do_remove=config["filtering"]["remove_blacklist"] and bool(blacklist_bed()),
        bl=blacklist_bed(),
    threads: config["resources"]["sort_threads"]
    log:
        f"{LOGS}/filter/{{sample}}.blacklist.log",
    conda:
        "../envs/peaks.yaml"
    shell:
        r"""
        if [ "{params.do_remove}" = "True" ] && [ -s "{params.bl}" ]; then
            bedtools intersect -v -abam {input.bam} -b {params.bl} \
                > {output.bam} 2> {log}
        else
            cp {input.bam} {output.bam}
        fi
        samtools index {output.bam}
        """


# Tn5-shifted, name-sorted BAM and a BEDPE of cut sites for MACS3/Genrich.
rule tn5_shift:
    input:
        bam=f"{RESULTS}/filtered/{{sample}}.filtered.bam",
    output:
        bam=f"{RESULTS}/shifted/{{sample}}.shifted.bam",
        bai=f"{RESULTS}/shifted/{{sample}}.shifted.bam.bai",
    params:
        do_shift=config["filtering"]["tn5_shift"],
    threads: config["resources"]["sort_threads"]
    log:
        f"{LOGS}/filter/{{sample}}.tn5shift.log",
    conda:
        "../envs/deeptools.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.bam})
        if [ "{params.do_shift}" = "True" ]; then
            # deepTools alignmentSieve applies the canonical +4/-5 Tn5 shift.
            alignmentSieve -b {input.bam} -o {output.bam}.tmp \
                --ATACshift -p {threads} 2> {log}
            samtools sort -@ {threads} -o {output.bam} {output.bam}.tmp 2>> {log}
            rm -f {output.bam}.tmp
        else
            cp {input.bam} {output.bam}
        fi
        samtools index {output.bam}
        """
