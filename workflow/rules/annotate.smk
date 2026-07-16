# =============================================================================
# annotate.smk — annotate consensus peaks with ChIPseeker
# =============================================================================
# Assigns each peak to its genomic feature (promoter/UTR/exon/intron/intergenic)
# and nearest gene using the genome's TxDb + OrgDb, and emits feature-
# distribution plots. Falls back to GTF-based annotation if no TxDb is matched.

rule chipseeker_annotate:
    input:
        bed=f"{RESULTS}/consensus/consensus_peaks.bed",
        gtf=genome_gtf(),
        scripts=script_inputs("annotate_peaks.R", "atac_theme.R"),
    output:
        tsv=f"{RESULTS}/annotation/consensus_peaks.annotated.tsv",
        plot=f"{RESULTS}/annotation/feature_distribution.pdf",
    params:
        txdb=GENOME.get("txdb", ""),
        orgdb=GENOME.get("orgdb", ""),
        up=config["annotation"]["tss_upstream"],
        down=config["annotation"]["tss_downstream"],
    log:
        f"{LOGS}/annotation/chipseeker.log",
    conda:
        "../envs/r.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.tsv})
        # Use only the conda env's R library and ignore ~/.Rprofile/.Renviron (see diffacc.smk).
        export R_LIBS_USER="$CONDA_PREFIX/lib/R/library"
        Rscript --vanilla workflow/scripts/annotate_peaks.R \
            --bed {input.bed} \
            --gtf {input.gtf} \
            --txdb "{params.txdb}" \
            --orgdb "{params.orgdb}" \
            --tss-up {params.up} --tss-down {params.down} \
            --out-tsv {output.tsv} --out-plot {output.plot} 2> {log}
        """
