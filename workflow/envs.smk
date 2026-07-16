# Data-independent installer workflow for every unique ATAC-seq rule environment.
# setup.sh uses this file only to create and smoke-test packages; it never reads
# biological data or writes pipeline results.

rule all:
    input:
        ".snakemake/setup-env-checks/align.ok",
        ".snakemake/setup-env-checks/counts.ok",
        ".snakemake/setup-env-checks/deeptools.ok",
        ".snakemake/setup-env-checks/deseq2.ok",
        ".snakemake/setup-env-checks/footprint.ok",
        ".snakemake/setup-env-checks/motif.ok",
        ".snakemake/setup-env-checks/peaks.ok",
        ".snakemake/setup-env-checks/qc.ok",
        ".snakemake/setup-env-checks/r.ok",
        ".snakemake/setup-env-checks/sra.ok",

rule check_align_env:
    output:
        touch(".snakemake/setup-env-checks/align.ok")
    conda:
        "envs/align.yaml"
    shell:
        """
        mkdir -p .snakemake/setup-env-checks
        for tool in bowtie2 bowtie2-build samtools picard; do command -v "$tool" >/dev/null; done
        bowtie2 --version >/dev/null
        samtools --version >/dev/null
        touch {output}
        """

rule check_counts_env:
    output:
        touch(".snakemake/setup-env-checks/counts.ok")
    conda:
        "envs/counts.yaml"
    shell:
        """
        mkdir -p .snakemake/setup-env-checks
        export R_LIBS_USER="$CONDA_PREFIX/lib/R/library"
        Rscript --vanilla -e 'suppressPackageStartupMessages(library(Rsubread)); stopifnot(as.character(packageVersion("Rsubread")) == "2.20.0")'
        python -c 'import csv'
        touch {output}
        """

rule check_deeptools_env:
    output:
        touch(".snakemake/setup-env-checks/deeptools.ok")
    conda:
        "envs/deeptools.yaml"
    shell:
        """
        mkdir -p .snakemake/setup-env-checks
        for tool in bamCoverage alignmentSieve bamPEFragmentSize computeMatrix plotFingerprint plotHeatmap plotProfile; do
            command -v "$tool" >/dev/null
        done
        python -c 'import numpy; assert numpy.__version__.startswith("1.26.")'
        touch {output}
        """

rule check_deseq2_env:
    output:
        touch(".snakemake/setup-env-checks/deseq2.ok")
    conda:
        "envs/deseq2.yaml"
    shell:
        """
        mkdir -p .snakemake/setup-env-checks
        export R_LIBS_USER="$CONDA_PREFIX/lib/R/library"
        Rscript --vanilla -e 'suppressPackageStartupMessages({{library(DESeq2); library(ashr); library(ggplot2); library(pheatmap)}})'
        touch {output}
        """

rule check_footprint_env:
    output:
        touch(".snakemake/setup-env-checks/footprint.ok")
    conda:
        "envs/footprint.yaml"
    shell:
        """
        mkdir -p .snakemake/setup-env-checks
        command -v TOBIAS >/dev/null
        TOBIAS --version >/dev/null
        touch {output}
        """

rule check_motif_env:
    output:
        touch(".snakemake/setup-env-checks/motif.ok")
    conda:
        "envs/motif.yaml"
    shell:
        """
        mkdir -p .snakemake/setup-env-checks
        for tool in findMotifsGenome.pl samtools perl; do command -v "$tool" >/dev/null; done
        touch {output}
        """

rule check_peaks_env:
    output:
        touch(".snakemake/setup-env-checks/peaks.ok")
    conda:
        "envs/peaks.yaml"
    shell:
        """
        mkdir -p .snakemake/setup-env-checks
        for tool in macs3 Genrich bedtools samtools; do command -v "$tool" >/dev/null; done
        macs3 --version >/dev/null
        python -c 'import numpy, pysam'
        touch {output}
        """

rule check_qc_env:
    output:
        touch(".snakemake/setup-env-checks/qc.ok")
    conda:
        "envs/qc.yaml"
    shell:
        """
        mkdir -p .snakemake/setup-env-checks
        for tool in fastqc fastp multiqc gzip; do command -v "$tool" >/dev/null; done
        fastqc --version >/dev/null
        fastp --version >/dev/null 2>&1
        multiqc --version >/dev/null
        python -c 'import pysam'
        touch {output}
        """

rule check_r_env:
    output:
        touch(".snakemake/setup-env-checks/r.ok")
    conda:
        "envs/r.yaml"
    shell:
        """
        mkdir -p .snakemake/setup-env-checks
        export R_LIBS_USER="$CONDA_PREFIX/lib/R/library"
        Rscript --vanilla -e 'suppressPackageStartupMessages({{library(DESeq2); library(ChIPseeker); library(clusterProfiler); library(enrichplot); library(chromVAR); library(AnnotationDbi); library(BiocParallel)}})'
        touch {output}
        """

rule check_sra_env:
    output:
        touch(".snakemake/setup-env-checks/sra.ok")
    conda:
        "envs/sra.yaml"
    shell:
        """
        mkdir -p .snakemake/setup-env-checks
        for tool in prefetch fastq-dump fasterq-dump curl wget gzip md5sum; do
            command -v "$tool" >/dev/null
        done
        prefetch --version >/dev/null
        fastq-dump --version >/dev/null
        fasterq-dump --version >/dev/null
        touch {output}
        """
