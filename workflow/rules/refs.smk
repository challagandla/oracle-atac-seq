# =============================================================================
# refs.smk — download & index the reference genome, blacklist, TSS annotation
# =============================================================================
# These rules only run if the corresponding files are not already supplied in
# config.yaml. Ensembl URLs are built from the genome preset.

def _ensembl_fasta_url():
    sp = GENOME["ensembl_species"]
    rel = GENOME["ensembl_release"]
    asm = GENOME["ensembl_assembly"]
    cap = sp.capitalize()
    return (
        f"https://ftp.ensembl.org/pub/release-{rel}/fasta/{sp}/dna/"
        f"{cap}.{asm}.dna_sm.primary_assembly.fa.gz"
    )


def _ensembl_gtf_url():
    sp = GENOME["ensembl_species"]
    rel = GENOME["ensembl_release"]
    asm = GENOME["ensembl_assembly"]
    cap = sp.capitalize()
    return (
        f"https://ftp.ensembl.org/pub/release-{rel}/gtf/{sp}/"
        f"{cap}.{asm}.{rel}.gtf.gz"
    )


rule download_genome_fasta:
    output:
        f"{REF}/genome.fa",
    params:
        url=_ensembl_fasta_url(),
    log:
        f"{LOGS}/refs/download_genome.log",
    conda:
        "../envs/sra.yaml"
    shell:
        r"""
        mkdir -p {REF}
        (curl -L --retry 3 -o {REF}/genome.fa.gz "{params.url}" \
          || wget -O {REF}/genome.fa.gz "{params.url}") 2> {log}
        gunzip -f {REF}/genome.fa.gz 2>> {log}
        """


rule download_gtf:
    output:
        f"{REF}/genes.gtf",
    params:
        url=_ensembl_gtf_url(),
    log:
        f"{LOGS}/refs/download_gtf.log",
    conda:
        "../envs/sra.yaml"
    shell:
        r"""
        mkdir -p {REF}
        (curl -L --retry 3 -o {REF}/genes.gtf.gz "{params.url}" \
          || wget -O {REF}/genes.gtf.gz "{params.url}") 2> {log}
        gunzip -f {REF}/genes.gtf.gz 2>> {log}
        """


rule download_blacklist:
    output:
        f"{REF}/blacklist.bed",
    params:
        url=GENOME.get("blacklist_url", ""),
    log:
        f"{LOGS}/refs/download_blacklist.log",
    conda:
        "../envs/sra.yaml"
    shell:
        r"""
        mkdir -p {REF}
        if [ -z "{params.url}" ]; then
            echo "No blacklist URL for this genome; writing empty file." > {log}
            : > {output}
        else
            (curl -L --retry 3 -o {REF}/blacklist.bed.gz "{params.url}" \
              || wget -O {REF}/blacklist.bed.gz "{params.url}") 2> {log}
            gunzip -f {REF}/blacklist.bed.gz 2>> {log}
        fi
        """


rule samtools_faidx:
    input:
        genome_fasta(),
    output:
        genome_fasta() + ".fai",
    log:
        f"{LOGS}/refs/faidx.log",
    conda:
        "../envs/align.yaml"
    shell:
        "samtools faidx {input} 2> {log}"


rule chrom_sizes:
    input:
        genome_fasta() + ".fai",
    output:
        f"{REF}/chrom.sizes",
    shell:
        "cut -f1,2 {input} > {output}"


rule bowtie2_build:
    input:
        genome_fasta(),
    output:
        expand(
            bowtie2_index_prefix() + ".{ext}",
            ext=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"],
        ),
    params:
        prefix=bowtie2_index_prefix(),
    threads: config["resources"]["align_threads"]
    log:
        f"{LOGS}/refs/bowtie2_build.log",
    conda:
        "../envs/align.yaml"
    shell:
        r"""
        mkdir -p $(dirname {params.prefix})
        bowtie2-build --threads {threads} {input} {params.prefix} 2> {log}
        """


# TSS BED (for TSS-enrichment QC) derived from the GTF.
rule tss_bed:
    input:
        gtf=genome_gtf(),
    output:
        f"{REF}/tss.bed",
    log:
        f"{LOGS}/refs/tss_bed.log",
    conda:
        "../envs/qc.yaml"
    shell:
        r"""
        python workflow/scripts/gtf_to_tss.py {input.gtf} {output} 2> {log}
        """
