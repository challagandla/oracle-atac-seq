# =============================================================================
# report.smk — ENCODE QC completeness + aggregate publication figures
# =============================================================================
# These rules add the QC metrics and cross-sample figures that turn a bag of
# per-sample outputs into a publication-ready report:
#   - library complexity (NRF / PBC1 / PBC2)         [ENCODE core ATAC QC]
#   - deepTools fingerprint (signal-to-noise)
#   - aggregate fragment-size / nucleosome overlay
#   - aggregate TSS-enrichment metagene heatmap + profile + numeric score
# Headline figures are emitted both as vector PDFs (results/figures, for papers)
# and as *_mqc.png / *_mqc.tsv so MultiQC folds them into the single report.
# Everything here is gated by config switches so smoke tests stay cheap.

# Colour-blind-safe palette (fixed order; mirrors atac_theme.R) for deepTools.
_PALETTE = ["#2a78d6", "#1baf7a", "#eda100", "#008300", "#4a3aa7", "#e34948",
            "#e87ba4", "#eb6834", "#6d4b9f", "#00a2b3", "#8c8c00", "#a6611a"]


def _sample_colours(samples):
    return [_PALETTE[i % len(_PALETTE)] for i in range(len(samples))]


# -----------------------------------------------------------------------------
# Library complexity: NRF / PBC1 / PBC2 from the pre-dedup alignment.
# -----------------------------------------------------------------------------
rule library_complexity:
    input:
        bam=f"{PROCESSED}/aligned/{{sample}}.sorted.bam",
    output:
        txt=f"{PROCESSED}/qc/complexity/{{sample}}.complexity.txt",
    log:
        f"{LOGS}/qc/complexity_{{sample}}.log",
    conda:
        "../envs/qc.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.txt})
        python workflow/scripts/library_complexity.py \
            --bam {input.bam} --sample {wildcards.sample} \
            --out {output.txt} 2> {log}
        """


rule library_complexity_table:
    input:
        expand(f"{PROCESSED}/qc/complexity/{{s}}.complexity.txt", s=SAMPLES),
    output:
        f"{PROCESSED}/qc/complexity/complexity_mqc.tsv",
    log:
        f"{LOGS}/qc/complexity_table.log",
    shell:
        r"""
        {{
          echo "# id: 'library_complexity'"
          echo "# section_name: 'Library complexity (NRF / PBC)'"
          echo "# description: 'Non-redundant fraction and PCR-bottleneck coefficients from the pre-dedup alignment. ENCODE: NRF>0.9, PBC1>0.9, PBC2>3 indicate a complex library.'"
          echo "# plot_type: 'table'"
          echo "# pconfig:"
          echo "#     id: 'library_complexity_table'"
          echo "#     namespace: 'ATAC'"
          echo -e "Sample\ttotal_frags\tdistinct\tNRF\tPBC1\tPBC2"
          for f in {input}; do tail -n +2 "$f"; done
        }} > {output} 2> {log}
        """


# -----------------------------------------------------------------------------
# deepTools fingerprint: cumulative read-enrichment (signal-to-noise) across
# all libraries, plus JS-distance quality metrics MultiQC can ingest.
# -----------------------------------------------------------------------------
rule plot_fingerprint:
    input:
        bams=expand(f"{PROCESSED}/filtered/{{s}}.filtered.bam", s=SAMPLES),
        bais=expand(f"{PROCESSED}/filtered/{{s}}.filtered.bam.bai", s=SAMPLES),
    output:
        plot=f"{PROCESSED}/qc/fingerprint/fingerprint_mqc.png",
        metrics=f"{PROCESSED}/qc/fingerprint/fingerprint_metrics.txt",
    params:
        labels=" ".join(SAMPLES),
    threads: config["resources"]["general_threads"]
    log:
        f"{LOGS}/qc/fingerprint.log",
    conda:
        "../envs/deeptools.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.plot})
        plotFingerprint -b {input.bams} --labels {params.labels} \
            --plotFile {output.plot} --outQualityMetrics {output.metrics} \
            --minMappingQuality 30 --skipZeros -p {threads} \
            --plotTitle "ATAC fingerprint (signal-to-noise)" 2> {log}
        """


# -----------------------------------------------------------------------------
# Aggregate fragment-size / nucleosome overlay.
# -----------------------------------------------------------------------------
rule fragment_size_overlay:
    input:
        frag=expand(f"{PROCESSED}/qc/fragsize/{{s}}.fragsize.txt", s=SAMPLES),
        samples=config["samples"],
    output:
        plot=f"{RESULTS}/figures/fragment_size_distribution.png",
        table=f"{PROCESSED}/qc/fragsize/fragment_partitions_mqc.tsv",
    log:
        f"{LOGS}/qc/fragment_overlay.log",
    conda:
        "../envs/deeptools.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.plot}) $(dirname {output.table})
        python workflow/scripts/plot_fragment_sizes.py \
            --frag {input.frag} --samples {input.samples} \
            --out-plot {output.plot} --out-table {output.table} 2> {log}
        # Surface the figure in MultiQC too.
        cp -f $(dirname {output.plot})/fragment_size_distribution_mqc.png \
              $(dirname {output.table})/ 2>> {log} || true
        """


# -----------------------------------------------------------------------------
# Aggregate TSS enrichment: one matrix over all libraries -> metagene heatmap,
# overlay profile, and a numeric per-sample enrichment score.
# -----------------------------------------------------------------------------
rule tss_matrix_all:
    input:
        bws=expand(f"{PROCESSED}/coverage/{{s}}.cpm.bw", s=SAMPLES),
        tss=f"{REF}/tss.bed",
    output:
        matrix=f"{PROCESSED}/qc/tss/tss_all.matrix.gz",
    params:
        labels=" ".join(SAMPLES),
    threads: config["resources"]["general_threads"]
    log:
        f"{LOGS}/qc/tss_matrix_all.log",
    conda:
        "../envs/deeptools.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.matrix})
        computeMatrix reference-point --referencePoint center \
            -S {input.bws} --samplesLabel {params.labels} \
            -R {input.tss} -a 2000 -b 2000 --skipZeros -p {threads} \
            -o {output.matrix} 2> {log}
        """


rule tss_metagene:
    input:
        matrix=f"{PROCESSED}/qc/tss/tss_all.matrix.gz",
    output:
        heatmap=f"{RESULTS}/figures/tss_enrichment_heatmap.pdf",
        profile=f"{RESULTS}/figures/tss_enrichment_profile.pdf",
        mqc=f"{PROCESSED}/qc/tss/tss_enrichment_heatmap_mqc.png",
    params:
        colours=lambda wc: " ".join(_sample_colours(SAMPLES)),
    log:
        f"{LOGS}/qc/tss_metagene.log",
    conda:
        "../envs/deeptools.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.heatmap}) $(dirname {output.mqc})
        plotHeatmap -m {input.matrix} -o {output.heatmap} \
            --refPointLabel TSS --colorMap Blues --heatmapHeight 14 \
            --plotTitle "TSS enrichment" 2> {log}
        plotHeatmap -m {input.matrix} -o {output.mqc} \
            --refPointLabel TSS --colorMap Blues --heatmapHeight 14 \
            --plotTitle "TSS enrichment" 2>> {log}
        plotProfile -m {input.matrix} -o {output.profile} \
            --refPointLabel TSS --perGroup --colors {params.colours} \
            --plotTitle "TSS enrichment profile" 2>> {log}
        """


rule tss_enrichment_score:
    input:
        matrix=f"{PROCESSED}/qc/tss/tss_all.matrix.gz",
    output:
        tsv=f"{PROCESSED}/qc/tss/tss_enrichment_mqc.tsv",
    log:
        f"{LOGS}/qc/tss_score.log",
    conda:
        "../envs/deeptools.yaml"
    shell:
        r"""
        python workflow/scripts/tss_score.py \
            --matrix {input.matrix} --out {output.tsv} 2> {log}
        """


# -----------------------------------------------------------------------------
# Collector: the report/QC targets enabled in config.
# -----------------------------------------------------------------------------
def report_outputs():
    out = []
    qc = config.get("qc", {})
    if qc.get("library_complexity", True):
        out += [f"{PROCESSED}/qc/complexity/complexity_mqc.tsv"]
    if qc.get("fingerprint", True):
        out += [f"{PROCESSED}/qc/fingerprint/fingerprint_mqc.png"]
    if config.get("report", {}).get("enabled", True):
        out += [f"{RESULTS}/figures/fragment_size_distribution.png"]
        if qc.get("tss", {}).get("enabled", True):
            out += [
                f"{RESULTS}/figures/tss_enrichment_heatmap.pdf",
                f"{RESULTS}/figures/tss_enrichment_profile.pdf",
                f"{PROCESSED}/qc/tss/tss_enrichment_mqc.tsv",
            ]
    return out
