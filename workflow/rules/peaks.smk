# =============================================================================
# peaks.smk — peak calling with MACS3 (per sample) and Genrich (per condition)
# =============================================================================
# MACS3 with --nomodel --shift -75 --extsize 150 is the de-facto standard for
# ATAC-seq (treats each read end as a Tn5 cut and builds 150 bp pileups around
# it). Genrich is run in ATAC mode (-j), pools replicates of a condition, and
# yields a single reproducible peak list — a useful cross-check.

rule macs3_callpeak:
    input:
        bam=f"{RESULTS}/shifted/{{sample}}.shifted.bam",
    output:
        narrowpeak=f"{RESULTS}/peaks/macs3/{{sample}}_peaks.narrowPeak",
        xls=f"{RESULTS}/peaks/macs3/{{sample}}_peaks.xls",
    params:
        gsize=macs_gsize(),
        q=config["peaks"]["macs3_qvalue"],
        extra=config["peaks"]["macs3_extra"],
        outdir=f"{RESULTS}/peaks/macs3",
        name=lambda wc: wc.sample,
    log:
        f"{LOGS}/peaks/macs3_{{sample}}.log",
    conda:
        "../envs/peaks.yaml"
    shell:
        r"""
        mkdir -p {params.outdir}
        macs3 callpeak -t {input.bam} -f BAMPE \
            -g {params.gsize} -q {params.q} {params.extra} \
            -n {params.name} --outdir {params.outdir} 2> {log}
        """


def genrich_input_bams(wildcards):
    # Genrich requires name-sorted BAMs; we sort the filtered BAMs by name.
    return expand(
        f"{RESULTS}/namesort/{{s}}.namesorted.bam",
        s=samples_in_condition(wildcards.cond),
    )


rule namesort_for_genrich:
    input:
        f"{RESULTS}/filtered/{{sample}}.filtered.bam",
    output:
        temp(f"{RESULTS}/namesort/{{sample}}.namesorted.bam"),
    threads: config["resources"]["sort_threads"]
    log:
        f"{LOGS}/peaks/namesort_{{sample}}.log",
    conda:
        "../envs/align.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output})
        samtools sort -n -@ {threads} -o {output} {input} 2> {log}
        """


rule genrich_callpeak:
    input:
        bams=genrich_input_bams,
    output:
        narrowpeak=f"{RESULTS}/peaks/genrich/{{cond}}.narrowPeak",
    params:
        extra=config["peaks"]["genrich_extra"],
        joined=lambda wc, input: ",".join(input.bams),
    log:
        f"{LOGS}/peaks/genrich_{{cond}}.log",
    conda:
        "../envs/peaks.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.narrowpeak})
        Genrich -t {params.joined} -o {output.narrowpeak} {params.extra} 2> {log}
        """
