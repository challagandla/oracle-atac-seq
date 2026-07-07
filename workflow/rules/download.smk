# =============================================================================
# download.smk — fetch paired FASTQs from SRA for sra-mode samples
# =============================================================================
# Uses the SRA Toolkit (prefetch + fasterq-dump) and gzips the result.
# Only samples whose row has an `sra` accession and no local fq1 trigger this.

SRA_SAMPLES = [s for s in SAMPLES if is_sra(s)]


rule sra_download:
    output:
        r1=f"{RAW}/{{sample}}_R1.fastq.gz",
        r2=f"{RAW}/{{sample}}_R2.fastq.gz",
    params:
        acc=lambda wc: sra_accession(wc.sample),
        outdir=RAW,
        tmp=lambda wc: f"{RAW}/_tmp_{wc.sample}",
        max_spots=config.get("sra", {}).get("max_spots", 0),
        max_size=config.get("sra", {}).get("max_size", "100G"),
    threads: config["resources"]["general_threads"]
    log:
        f"{LOGS}/download/{{sample}}.log",
    conda:
        "../envs/sra.yaml"
    shell:
        r"""
        mkdir -p {params.outdir} {params.tmp}
        if [ "{params.max_spots}" != "0" ]; then
            fastq-dump --split-files --skip-technical \
                --maxSpotId {params.max_spots} \
                -O {params.tmp} {params.acc} 2> {log}
        else
            prefetch --max-size {params.max_size} -O {params.tmp} {params.acc} 2> {log}
            fasterq-dump --split-files --threads {threads} \
                -O {params.tmp} {params.tmp}/{params.acc}/{params.acc}.sra 2>> {log}
        fi
        if [ ! -s {params.tmp}/{params.acc}_1.fastq ] || [ ! -s {params.tmp}/{params.acc}_2.fastq ]; then
            echo "Expected paired FASTQ outputs for {params.acc}, but one mate file is missing or empty." >> {log}
            exit 1
        fi
        gzip -c {params.tmp}/{params.acc}_1.fastq > {output.r1} 2>> {log}
        gzip -c {params.tmp}/{params.acc}_2.fastq > {output.r2} 2>> {log}
        rm -rf {params.tmp}
        """
