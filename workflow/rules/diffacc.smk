# =============================================================================
# diffacc.smk — differential accessibility with DESeq2
# =============================================================================
# DESeq2 on the consensus peak x sample count matrix. The design formula and
# contrast come from config.yaml and must reference columns in samples.tsv.
# Outputs a results table, MA plot, PCA, and a normalised-counts matrix.

rule deseq2:
    input:
        counts=f"{RESULTS}/counts/consensus_counts.tsv",
        samples=config["samples"],
    output:
        tsv=f"{RESULTS}/diffacc/differential_accessibility.tsv",
        ma=f"{RESULTS}/diffacc/MA_plot.pdf",
        pca=f"{RESULTS}/diffacc/PCA_plot.pdf",
        norm=f"{RESULTS}/diffacc/normalized_counts.tsv",
        up=f"{RESULTS}/diffacc/up_peaks.bed",
        down=f"{RESULTS}/diffacc/down_peaks.bed",
    params:
        design=config["diffacc"]["design"],
        factor=config["diffacc"]["contrast"][0],
        numerator=config["diffacc"]["contrast"][1],
        denominator=config["diffacc"]["contrast"][2],
        fdr=config["diffacc"]["fdr"],
        lfc=config["diffacc"]["lfc_threshold"],
    log:
        f"{LOGS}/diffacc/deseq2.log",
    conda:
        "../envs/r.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.tsv})
        Rscript workflow/scripts/differential_accessibility.R \
            --counts {input.counts} \
            --samples {input.samples} \
            --design "{params.design}" \
            --factor "{params.factor}" \
            --numerator "{params.numerator}" \
            --denominator "{params.denominator}" \
            --fdr {params.fdr} \
            --lfc {params.lfc} \
            --outdir $(dirname {output.tsv}) 2> {log}
        """
