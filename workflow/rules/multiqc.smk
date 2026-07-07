# =============================================================================
# multiqc.smk — aggregate all QC into one HTML report
# =============================================================================
# MultiQC scans the results tree and assembles FastQC, fastp, Bowtie2, Picard,
# featureCounts, fragment-size and FRiP metrics into a single report.

def multiqc_inputs(wildcards):
    inp = []
    inp += expand(f"{PROCESSED}/qc/fastqc/{{s}}_{{r}}_fastqc.zip", s=SAMPLES, r=["R1", "R2"])
    inp += expand(f"{PROCESSED}/qc/fastp/{{s}}.fastp.json", s=SAMPLES)
    inp += expand(f"{PROCESSED}/qc/flagstat/{{s}}.raw.flagstat", s=SAMPLES)
    inp += expand(f"{PROCESSED}/qc/picard/{{s}}.dup_metrics.txt", s=SAMPLES)
    inp += [f"{PROCESSED}/qc/frip/frip_mqc.tsv"]
    inp += expand(f"{PROCESSED}/qc/fragsize/{{s}}.fragsize.txt", s=SAMPLES)
    inp += expand(f"{PROCESSED}/qc/tss/{{s}}.tss_enrichment.png", s=SAMPLES)
    inp += [f"{RESULTS}/counts/consensus_counts.tsv.summary"]

    # --- ENCODE QC completeness + aggregate figures (report.smk) --------------
    qc = config.get("qc", {})
    if qc.get("library_complexity", True):
        inp += [f"{PROCESSED}/qc/complexity/complexity_mqc.tsv"]
    if qc.get("fingerprint", True):
        inp += [f"{PROCESSED}/qc/fingerprint/fingerprint_metrics.txt"]
    if config.get("report", {}).get("enabled", True):
        inp += [f"{PROCESSED}/qc/fragsize/fragment_partitions_mqc.tsv"]
        if qc.get("tss", {}).get("enabled", True):
            inp += [f"{PROCESSED}/qc/tss/tss_enrichment_mqc.tsv"]
            inp += [f"{PROCESSED}/qc/tss/tss_enrichment_heatmap_mqc.png"]

    # --- differential-accessibility headline figures -------------------------
    if config["diffacc"].get("enabled", False):
        inp += expand(
            f"{RESULTS}/diffacc/{{f}}_mqc.png",
            f=["PCA_plot", "MA_plot", "volcano_plot",
               "sample_correlation_heatmap", "differential_peaks_heatmap"],
        )
    return inp


rule multiqc:
    input:
        multiqc_inputs,
    output:
        f"{RESULTS}/qc/multiqc_report.html",
    params:
        scan_dirs=f"{PROCESSED} {RESULTS}",
        outdir=f"{RESULTS}/qc",
    log:
        f"{LOGS}/qc/multiqc.log",
    conda:
        "../envs/qc.yaml"
    shell:
        r"""
        multiqc {params.scan_dirs} -f -o {params.outdir} \
            -n multiqc_report.html 2> {log}
        """
