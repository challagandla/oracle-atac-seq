# =============================================================================
# coverage.smk — CPM-normalised bigWig signal tracks (deepTools bamCoverage)
# =============================================================================
# Browser bigWigs represent complete paired fragments at 10 bp resolution.
# Separate 1 bp Tn5 insertion tracks drive TSS enrichment and can be loaded in
# IGV when cut-site resolution is needed.

rule bamcoverage:
    input:
        bam=f"{PROCESSED}/filtered/{{sample}}.filtered.bam",
        bai=f"{PROCESSED}/filtered/{{sample}}.filtered.bam.bai",
    output:
        bw=f"{PROCESSED}/coverage/{{sample}}.cpm.bw",
    params:
        egs=effective_genome_size(),
    threads: config["resources"]["general_threads"]
    log:
        f"{LOGS}/coverage/{{sample}}.log",
    conda:
        "../envs/deeptools.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.bw})
        bamCoverage -b {input.bam} -o {output.bw} \
            --binSize 10 --normalizeUsing CPM \
            --effectiveGenomeSize {params.egs} \
            --extendReads --samFlagInclude 64 -p {threads} 2> {log}
        """


rule cutsite_coverage:
    input:
        bam=f"{PROCESSED}/shifted/{{sample}}.shifted.bam",
        bai=f"{PROCESSED}/shifted/{{sample}}.shifted.bam.bai",
    output:
        bw=f"{PROCESSED}/coverage/{{sample}}.cutsites.cpm.bw",
    params:
        egs=effective_genome_size(),
    threads: config["resources"]["general_threads"]
    log:
        f"{LOGS}/coverage/{{sample}}.cutsites.log",
    conda:
        "../envs/deeptools.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.bw})
        bamCoverage -b {input.bam} -o {output.bw} \
            --binSize 1 --Offset 1 --normalizeUsing CPM \
            --effectiveGenomeSize {params.egs} -p {threads} 2> {log}
        """
