# =============================================================================
# coverage.smk — CPM-normalised bigWig signal tracks (deepTools bamCoverage)
# =============================================================================
# bigWigs are built from the Tn5-shifted BAMs, normalised to counts-per-million
# with the genome's effective size, at 10 bp resolution. Load them in IGV/UCSC.

rule bamcoverage:
    input:
        bam=f"{RESULTS}/shifted/{{sample}}.shifted.bam",
        bai=f"{RESULTS}/shifted/{{sample}}.shifted.bam.bai",
    output:
        bw=f"{RESULTS}/coverage/{{sample}}.cpm.bw",
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
            --extendReads -p {threads} 2> {log}
        """
