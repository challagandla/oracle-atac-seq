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
        f"{REF}/blacklist.raw.bed",
    params:
        url=GENOME.get("blacklist_url", ""),
        md5=GENOME.get("blacklist_md5", ""),
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
            (curl -L --retry 3 -o {REF}/blacklist.raw.bed.gz "{params.url}" \
              || wget -O {REF}/blacklist.raw.bed.gz "{params.url}") 2> {log}
            if [ -n "{params.md5}" ]; then
                echo "{params.md5}  {REF}/blacklist.raw.bed.gz" \
                  | md5sum --check --status - || {{
                    echo "blacklist checksum verification failed" >> {log}; exit 1;
                  }}
            fi
            gunzip -f {REF}/blacklist.raw.bed.gz 2>> {log}
        fi
        """


# The ENCODE blacklists ship with UCSC names (chr1); Ensembl genomes call the
# same sequence 1. bedtools matches on the name, so the mismatch makes blacklist
# filtering a silent no-op. Rename to whatever the genome calls its chromosomes,
# and fail if the two cannot be reconciled.
rule harmonize_blacklist:
    input:
        blacklist=blacklist_source(),
        chrom=f"{REF}/chrom.sizes",
        scripts=script_inputs("harmonize_blacklist.py"),
    output:
        f"{REF}/blacklist.harmonized.bed",
    log:
        f"{LOGS}/refs/harmonize_blacklist.log",
    conda:
        "../envs/peaks.yaml"
    shell:
        r"""
        python workflow/scripts/harmonize_blacklist.py \
            --blacklist {input.blacklist} \
            --chrom-sizes {input.chrom} \
            --out {output} 2> {log}
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
        r"""
        mkdir -p {REF}
        cut -f1,2 {input} > {output}
        """


rule bowtie2_build:
    input:
        genome_fasta(),
    output:
        bowtie2_index_files(),
    params:
        prefix=bowtie2_index_prefix(),
        mode=bowtie2_build_mode(),
    threads: config["resources"]["align_threads"]
    log:
        f"{LOGS}/refs/bowtie2_build.log",
    conda:
        "../envs/align.yaml"
    shell:
        r"""
        mkdir -p $(dirname {params.prefix})
        bowtie2-build {params.mode} --threads {threads} {input} {params.prefix} 2> {log}
        """


# TSS BED (for TSS-enrichment QC) derived from the GTF.
rule tss_bed:
    input:
        gtf=genome_gtf(),
        chrom=f"{REF}/chrom.sizes",
        scripts=script_inputs("gtf_to_tss.py"),
    output:
        f"{REF}/tss.bed",
    log:
        f"{LOGS}/refs/tss_bed.log",
    conda:
        "../envs/qc.yaml"
    shell:
        r"""
        mkdir -p {REF}
        python workflow/scripts/gtf_to_tss.py {input.gtf} {output} \
            --chrom-sizes {input.chrom} 2> {log}
        test -s {output}
        """
