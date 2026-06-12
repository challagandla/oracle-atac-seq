# =============================================================================
# align.smk — Bowtie2 alignment -> coordinate-sorted, indexed BAM
# =============================================================================
# Paired-end, -X 2000 to retain nucleosomal fragments. Output is the raw
# alignment BAM; all filtering happens in filter.smk so QC can compare before/
# after.

rule bowtie2_align:
    input:
        r1=f"{RESULTS}/trimmed/{{sample}}_R1.trimmed.fastq.gz",
        r2=f"{RESULTS}/trimmed/{{sample}}_R2.trimmed.fastq.gz",
        idx=expand(
            bowtie2_index_prefix() + ".{ext}",
            ext=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"],
        ),
    output:
        bam=f"{RESULTS}/aligned/{{sample}}.sorted.bam",
        bai=f"{RESULTS}/aligned/{{sample}}.sorted.bam.bai",
    params:
        prefix=bowtie2_index_prefix(),
        extra=config["alignment"]["bowtie2_extra"],
    threads: config["resources"]["align_threads"]
    log:
        bt2=f"{LOGS}/align/{{sample}}.bowtie2.log",
        sort=f"{LOGS}/align/{{sample}}.sort.log",
    conda:
        "../envs/align.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.bam})
        bowtie2 -x {params.prefix} -1 {input.r1} -2 {input.r2} \
            --rg-id {wildcards.sample} \
            --rg SM:{wildcards.sample} --rg LB:{wildcards.sample} --rg PL:ILLUMINA \
            -p {threads} {params.extra} 2> {log.bt2} \
          | samtools sort -@ {threads} -o {output.bam} - 2> {log.sort}
        samtools index {output.bam}
        """


rule flagstat_raw:
    input:
        f"{RESULTS}/aligned/{{sample}}.sorted.bam",
    output:
        f"{RESULTS}/qc/flagstat/{{sample}}.raw.flagstat",
    conda:
        "../envs/align.yaml"
    shell:
        "samtools flagstat {input} > {output}"
