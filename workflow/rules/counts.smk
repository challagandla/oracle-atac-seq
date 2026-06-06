# =============================================================================
# counts.smk — count fragments in consensus peaks (featureCounts)
# =============================================================================
# featureCounts (Subread) in paired-end fragment-counting mode produces the
# peak-by-sample matrix that feeds DESeq2 and chromVAR.

rule featurecounts:
    input:
        saf=f"{RESULTS}/consensus/consensus_peaks.saf",
        bams=expand(f"{RESULTS}/filtered/{{s}}.filtered.bam", s=SAMPLES),
    output:
        counts=f"{RESULTS}/counts/consensus_counts.tsv",
        summary=f"{RESULTS}/counts/consensus_counts.tsv.summary",
    threads: config["resources"]["general_threads"]
    log:
        f"{LOGS}/counts/featurecounts.log",
    conda:
        "../envs/peaks.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.counts})
        featureCounts -p --countReadPairs -F SAF \
            -a {input.saf} -o {output.counts} \
            -T {threads} {input.bams} 2> {log}
        """
