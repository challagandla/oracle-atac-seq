# =============================================================================
# download.smk — fetch paired FASTQs from SRA for sra-mode samples
# =============================================================================
# Uses the SRA Toolkit (prefetch + fasterq-dump) and gzips the result.
# Only samples whose row has an `sra` accession and no local fq1 trigger this.

SRA_SAMPLES = [s for s in SAMPLES if is_sra(s)]


rule sra_download:
    output:
        r1=f"{RESULTS}/fastq/{{sample}}_R1.fastq.gz",
        r2=f"{RESULTS}/fastq/{{sample}}_R2.fastq.gz",
    params:
        acc=lambda wc: sra_accession(wc.sample),
        outdir=f"{RESULTS}/fastq",
        tmp=lambda wc: f"{RESULTS}/fastq/_tmp_{wc.sample}",
    threads: config["resources"]["general_threads"]
    log:
        f"{LOGS}/download/{{sample}}.log",
    conda:
        "../envs/sra.yaml"
    shell:
        r"""
        mkdir -p {params.outdir} {params.tmp}
        prefetch --max-size 100G -O {params.tmp} {params.acc} 2> {log}
        fasterq-dump --split-files --threads {threads} \
            -O {params.tmp} {params.tmp}/{params.acc}/{params.acc}.sra 2>> {log}
        gzip -c {params.tmp}/{params.acc}_1.fastq > {output.r1} 2>> {log}
        gzip -c {params.tmp}/{params.acc}_2.fastq > {output.r2} 2>> {log}
        rm -rf {params.tmp}
        """
