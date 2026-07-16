# =============================================================================
# multiqc.smk — aggregate all QC into one HTML report
# =============================================================================
# MultiQC receives an exact manifest for the active sample sheet. It never scans
# a broad results directory, so stale analyses cannot leak into this report.

def qc_review_inputs(wildcards):
    """Outcome-blind library QC inputs; never include differential results."""
    inp = []
    inp += expand(f"{PROCESSED}/qc/fastqc/{{s}}_{{r}}_fastqc.zip", s=SAMPLES, r=["R1", "R2"])
    inp += expand(f"{PROCESSED}/qc/fastp/{{s}}.fastp.json", s=SAMPLES)
    inp += expand(f"{PROCESSED}/qc/flagstat/{{s}}.raw.flagstat", s=SAMPLES)
    inp += expand(f"{PROCESSED}/qc/flagstat/{{s}}.filtered.flagstat", s=SAMPLES)
    inp += expand(f"{LOGS}/align/{{s}}.bowtie2.log", s=SAMPLES)
    inp += expand(f"{PROCESSED}/qc/picard/{{s}}.dup_metrics.txt", s=SAMPLES)
    inp += [f"{PROCESSED}/qc/frip/frip_mqc.tsv"]
    inp += expand(f"{PROCESSED}/qc/fragsize/{{s}}.fragsize.txt", s=SAMPLES)
    if config.get("qc", {}).get("tss", {}).get("enabled", True):
        inp += expand(f"{PROCESSED}/qc/tss/{{s}}.tss_enrichment.png", s=SAMPLES)
    # --- ENCODE QC completeness + aggregate figures (report.smk) --------------
    qc = config.get("qc", {})
    if qc.get("library_complexity", True):
        inp += [f"{PROCESSED}/qc/complexity/complexity_mqc.tsv"]
    inp += [f"{PROCESSED}/qc/summary/qc_decisions_mqc.tsv"]
    if qc.get("fingerprint", True):
        inp += [f"{PROCESSED}/qc/fingerprint/fingerprint_metrics.txt"]
    if qc.get("sample_similarity", True):
        inp += [
            f"{PROCESSED}/qc/similarity/qc_sample_correlation_mqc.png",
            f"{PROCESSED}/qc/similarity/qc_sample_pca_mqc.png",
        ]
    if config.get("report", {}).get("enabled", True):
        inp += [f"{PROCESSED}/qc/fragsize/fragment_partitions_mqc.tsv"]
        inp += [f"{RESULTS}/figures/fragment_size_distribution_mqc.png"]
        if qc.get("tss", {}).get("enabled", True):
            inp += [f"{PROCESSED}/qc/tss/tss_enrichment_mqc.tsv"]
            inp += [f"{PROCESSED}/qc/tss/tss_enrichment_heatmap_mqc.png"]
    return inp


def multiqc_inputs(wildcards):
    """Final combined report inputs, after the included sample set is frozen."""
    inp = list(qc_review_inputs(wildcards))
    inp += [f"{RESULTS}/counts/consensus_counts.tsv.summary"]

    # --- differential-accessibility headline figures -------------------------
    if config["diffacc"].get("enabled", False):
        inp += expand(
            f"{RESULTS}/diffacc/{{f}}_mqc.png",
            f=["PCA_plot", "MA_plot", "volcano_plot",
               "sample_correlation_heatmap", "differential_peaks_heatmap"],
        )
    if config.get("chromvar", {}).get("enabled", False):
        inp += [
            f"{RESULTS}/chromvar/chromvar_variability_mqc.png",
            f"{RESULTS}/chromvar/chromvar_deviation_heatmap_mqc.png",
        ]
    return inp


rule qc_review_manifest:
    input:
        qc_review_inputs,
    output:
        f"{RESULTS}/qc/qc_review_file_list.txt",
    run:
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        paths = list(dict.fromkeys(map(str, input)))
        with open(output[0], "w") as handle:
            handle.write("\n".join(paths) + "\n")


rule qc_review:
    input:
        manifest=f"{RESULTS}/qc/qc_review_file_list.txt",
    output:
        f"{RESULTS}/qc/qc_review_report.html",
    params:
        outdir=f"{RESULTS}/qc",
    log:
        f"{LOGS}/qc/qc_review.log",
    conda:
        "../envs/qc.yaml"
    shell:
        r"""
        multiqc --file-list {input.manifest} -f -o {params.outdir} \
            -n qc_review_report.html 2> {log}
        """


rule multiqc_manifest:
    input:
        multiqc_inputs,
    output:
        f"{RESULTS}/qc/multiqc_file_list.txt",
    run:
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        paths = list(dict.fromkeys(map(str, input)))
        with open(output[0], "w") as handle:
            handle.write("\n".join(paths) + "\n")


rule multiqc:
    input:
        manifest=f"{RESULTS}/qc/multiqc_file_list.txt",
    output:
        f"{RESULTS}/qc/multiqc_report.html",
    params:
        outdir=f"{RESULTS}/qc",
    log:
        f"{LOGS}/qc/multiqc.log",
    conda:
        "../envs/qc.yaml"
    shell:
        r"""
        multiqc --file-list {input.manifest} -f -o {params.outdir} \
            -n multiqc_report.html 2> {log}
        """
