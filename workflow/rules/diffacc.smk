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
        scripts=script_inputs("differential_accessibility.R", "atac_theme.R"),
    output:
        tsv=f"{RESULTS}/diffacc/differential_accessibility.tsv",
        ma=f"{RESULTS}/diffacc/MA_plot.pdf",
        pca=f"{RESULTS}/diffacc/PCA_plot.pdf",
        norm=f"{RESULTS}/diffacc/normalized_counts.tsv",
        tested=f"{RESULTS}/diffacc/tested_peaks.bed",
        up=f"{RESULTS}/diffacc/up_peaks.bed",
        down=f"{RESULTS}/diffacc/down_peaks.bed",
        volcano=f"{RESULTS}/diffacc/volcano_plot.pdf",
        corr=f"{RESULTS}/diffacc/sample_correlation_heatmap.pdf",
        dist=f"{RESULTS}/diffacc/sample_distance_heatmap.pdf",
        scree=f"{RESULTS}/diffacc/scree_plot.pdf",
        heat=f"{RESULTS}/diffacc/differential_peaks_heatmap.pdf",
        pval=f"{RESULTS}/diffacc/pvalue_histogram.pdf",
        summary=f"{RESULTS}/diffacc/diffacc_summary.tsv",
        png=expand(
            f"{RESULTS}/diffacc/{{f}}.png",
            f=["PCA_plot", "MA_plot", "volcano_plot", "scree_plot",
               "sample_correlation_heatmap", "sample_distance_heatmap",
               "differential_peaks_heatmap", "pvalue_histogram"],
        ),
        # Headline figures re-exported as *_mqc.png for the MultiQC report.
        mqc=expand(
            f"{RESULTS}/diffacc/{{f}}_mqc.png",
            f=["PCA_plot", "MA_plot", "volcano_plot",
               "sample_correlation_heatmap", "differential_peaks_heatmap"],
        ),
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
        "../envs/deseq2.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.tsv})
        # Use only the conda env's R library and ignore ~/.Rprofile/.Renviron, so a
        # host R_LIBS_USER or .libPaths() can't shadow env packages with an
        # ABI-incompatible build (e.g. SummarizedExperiment/Biobase load failures).
        export R_LIBS_USER="$CONDA_PREFIX/lib/R/library"
        Rscript --vanilla workflow/scripts/differential_accessibility.R \
            --counts {input.counts} \
            --samples {input.samples} \
            --design {params.design:q} \
            --factor {params.factor:q} \
            --numerator {params.numerator:q} \
            --denominator {params.denominator:q} \
            --fdr {params.fdr} \
            --lfc {params.lfc} \
            --outdir $(dirname {output.tsv}) 2> {log}
        """
