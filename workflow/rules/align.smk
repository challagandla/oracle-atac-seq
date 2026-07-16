# =============================================================================
# align.smk — Bowtie2 alignment -> coordinate-sorted, indexed BAM
# =============================================================================
# Paired-end, -X 2000 to retain nucleosomal fragments. Output is the raw
# alignment BAM; all filtering happens in filter.smk so QC can compare before/
# after.

rule bowtie2_align:
    input:
        r1=f"{PROCESSED}/trimmed/{{sample}}_R1.trimmed.fastq.gz",
        r2=f"{PROCESSED}/trimmed/{{sample}}_R2.trimmed.fastq.gz",
        idx=bowtie2_index_files(),
    output:
        bam=f"{PROCESSED}/aligned/{{sample}}.sorted.bam",
        bai=f"{PROCESSED}/aligned/{{sample}}.sorted.bam.bai",
    params:
        prefix=bowtie2_index_prefix(),
        extra=config["alignment"]["bowtie2_extra"],
        tmp_root=f"{PROCESSED}/aligned/.tmp",
    threads: config["resources"]["align_threads"]
    log:
        bt2=f"{LOGS}/align/{{sample}}.bowtie2.log",
        sort=f"{LOGS}/align/{{sample}}.sort.log",
    conda:
        "../envs/align.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.bam}) {params.tmp_root}
        tmpdir=$(mktemp -d "{params.tmp_root}/{wildcards.sample}.XXXXXX")
        cleanup() {{ rm -rf -- "$tmpdir"; }}
        trap cleanup EXIT
        bowtie2 -x {params.prefix} -1 {input.r1} -2 {input.r2} \
            --rg-id {wildcards.sample} \
            --rg SM:{wildcards.sample} --rg LB:{wildcards.sample} --rg PL:ILLUMINA \
            -p {threads} {params.extra} 2> {log.bt2} \
          | samtools sort -@ {threads} -T "$tmpdir/chunk" \
                -o {output.bam} - 2> {log.sort}
        samtools quickcheck {output.bam} 2>> {log.sort}
        samtools index {output.bam}
        """


rule flagstat_raw:
    input:
        f"{PROCESSED}/aligned/{{sample}}.sorted.bam",
    output:
        f"{PROCESSED}/qc/flagstat/{{sample}}.raw.flagstat",
    conda:
        "../envs/align.yaml"
    shell:
        "samtools flagstat {input} > {output}"
