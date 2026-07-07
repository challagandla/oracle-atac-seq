# =============================================================================
# chromvar.smk — per-sample motif accessibility deviations (optional)
# =============================================================================
# chromVAR (Schep et al. 2017, Nat Methods) computes bias-corrected motif
# accessibility deviations from the consensus peak count matrix. Useful for
# ranking TF motifs by sample-to-sample variability without footprinting.

rule chromvar:
    input:
        counts=f"{RESULTS}/counts/consensus_counts.tsv",
        bed=f"{RESULTS}/consensus/consensus_peaks.bed",
        samples=config["samples"],
    output:
        dev=f"{RESULTS}/chromvar/chromvar_deviations.tsv",
        var=f"{RESULTS}/chromvar/chromvar_variability.tsv",
        varplot=f"{RESULTS}/chromvar/chromvar_variability.pdf",
        heatmap=f"{RESULTS}/chromvar/chromvar_deviation_heatmap.pdf",
    params:
        bsgenome=GENOME.get("bsgenome", ""),
        taxid=GENOME.get("taxid", 9606),
        top_n=config.get("chromvar", {}).get("top_n", 30),
    log:
        f"{LOGS}/chromvar/chromvar.log",
    conda:
        "../envs/r.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.dev})
        # Use only the conda env's R library and ignore ~/.Rprofile/.Renviron (see diffacc.smk).
        export R_LIBS_USER="$CONDA_PREFIX/lib/R/library"
        Rscript --vanilla workflow/scripts/chromvar_analysis.R \
            --counts {input.counts} \
            --bed {input.bed} \
            --samples {input.samples} \
            --bsgenome "{params.bsgenome}" \
            --taxid {params.taxid} --top-n {params.top_n} \
            --out-dev {output.dev} --out-var {output.var} \
            --out-var-plot {output.varplot} --out-heatmap {output.heatmap} 2> {log}
        """
