# =============================================================================
# enrich.smk — GO/KEGG functional enrichment of differentially-accessible genes
# =============================================================================
# Joins the ChIPseeker peak->gene annotation with the DESeq2 up/down peak sets
# and runs clusterProfiler (enrichGO offline via OrgDb; enrichKEGG online).
# Produces a multi-panel dot-plot PDF plus per-set result tables.

rule functional_enrichment:
    input:
        annotation=f"{RESULTS}/annotation/consensus_peaks.annotated.tsv",
        tested=f"{RESULTS}/diffacc/tested_peaks.bed",
        up=f"{RESULTS}/diffacc/up_peaks.bed",
        down=f"{RESULTS}/diffacc/down_peaks.bed",
        scripts=script_inputs("functional_enrichment.R", "atac_theme.R"),
    output:
        pdf=f"{RESULTS}/enrichment/enrichment_dotplots.pdf",
        status=f"{RESULTS}/enrichment/enrichment_status.tsv",
        tables=directory(f"{RESULTS}/enrichment/tables"),
    params:
        orgdb=GENOME.get("orgdb", ""),
        taxid=GENOME.get("taxid", 0),
        ontologies=",".join(config.get("functional_enrichment", {}).get("ontologies", ["BP", "MF"])),
        kegg=str(config.get("functional_enrichment", {}).get("kegg", False)).lower(),
        qvalue=config.get("functional_enrichment", {}).get("qvalue", 0.05),
        top_n=config.get("functional_enrichment", {}).get("top_n", 20),
    log:
        f"{LOGS}/enrichment/functional_enrichment.log",
    conda:
        "../envs/r.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.pdf})
        rm -rf {output.tables}
        mkdir -p {output.tables}
        # Use only the conda env's R library and ignore ~/.Rprofile/.Renviron (see diffacc.smk).
        export R_LIBS_USER="$CONDA_PREFIX/lib/R/library"
        Rscript --vanilla workflow/scripts/functional_enrichment.R \
            --annotation {input.annotation} \
            --tested {input.tested} \
            --up {input.up} --down {input.down} \
            --orgdb "{params.orgdb}" --taxid {params.taxid} \
            --ontologies "{params.ontologies}" --kegg {params.kegg} \
            --qvalue {params.qvalue} --top-n {params.top_n} \
            --out-pdf {output.pdf} --out-status {output.status} \
            --out-dir {output.tables} 2> {log}
        """
