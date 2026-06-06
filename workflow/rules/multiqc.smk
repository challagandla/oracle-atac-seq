# =============================================================================
# multiqc.smk — aggregate all QC into one HTML report
# =============================================================================
# MultiQC scans the results tree and assembles FastQC, fastp, Bowtie2, Picard,
# featureCounts, fragment-size and FRiP metrics into a single report.

def multiqc_inputs(wildcards):
    inp = []
    inp += expand(f"{RESULTS}/qc/fastqc/{{s}}_{{r}}_fastqc.zip", s=SAMPLES, r=["R1", "R2"])
    inp += expand(f"{RESULTS}/qc/fastp/{{s}}.fastp.json", s=SAMPLES)
    inp += expand(f"{RESULTS}/qc/flagstat/{{s}}.raw.flagstat", s=SAMPLES)
    inp += expand(f"{RESULTS}/qc/picard/{{s}}.dup_metrics.txt", s=SAMPLES)
    inp += expand(f"{RESULTS}/qc/frip/{{s}}.frip.txt", s=SAMPLES)
    inp += expand(f"{RESULTS}/qc/fragsize/{{s}}.fragsize.txt", s=SAMPLES)
    inp += expand(f"{RESULTS}/qc/tss/{{s}}.tss_enrichment.png", s=SAMPLES)
    inp += [f"{RESULTS}/counts/consensus_counts.tsv.summary"]
    return inp


rule multiqc:
    input:
        multiqc_inputs,
    output:
        f"{RESULTS}/qc/multiqc_report.html",
    params:
        scan_dir=RESULTS,
        outdir=f"{RESULTS}/qc",
    log:
        f"{LOGS}/qc/multiqc.log",
    conda:
        "../envs/qc.yaml"
    shell:
        r"""
        multiqc {params.scan_dir} -f -o {params.outdir} \
            -n multiqc_report.html 2> {log}
        """
