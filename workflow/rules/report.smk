# =============================================================================
# report.smk — ENCODE-informed QC + aggregate figures
# =============================================================================
# These rules add the QC metrics and cross-sample figures that turn a bag of
# per-sample outputs into a reviewable report:
#   - library complexity (NRF / PBC1 / PBC2)         [ENCODE core ATAC QC]
#   - deepTools fingerprint (signal-to-noise)
#   - aggregate fragment-size / nucleosome overlay
#   - aggregate TSS-enrichment metagene heatmap + profile + numeric score
# Headline figures are emitted both as vector PDFs (results/figures, for papers)
# and as *_mqc.png / *_mqc.tsv so MultiQC folds them into the single report.
# Everything here is gated by config switches so smoke tests stay cheap.

import os
import sys

# Colours come from workflow/scripts/palette.py, the single source of truth that
# atac_theme.R and plot_fragment_sizes.py also follow. Keyed by condition, so a
# hue means the same thing in every panel of the report.
sys.path.insert(0, os.path.join(workflow.basedir, "scripts"))
from palette import SEQUENTIAL_CMAP, shell_colours  # noqa: E402

_COND_OF = dict(zip(samples["sample"], samples["condition"]))


# -----------------------------------------------------------------------------
# Library complexity: NRF / PBC1 / PBC2 from the pre-dedup alignment.
# -----------------------------------------------------------------------------
rule library_complexity:
    input:
        # Proper-pair/MAPQ-filtered, nuclear, but not yet duplicate-removed.
        bam=f"{PROCESSED}/filtered/{{sample}}.nomito.bam",
        # library_complexity.py uses indexed pysam.fetch(); declaring the
        # temporary index keeps Snakemake from removing it after remove_mito.
        bai=f"{PROCESSED}/filtered/{{sample}}.nomito.bam.bai",
        scripts=script_inputs("library_complexity.py"),
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
            --min-mapq {config[filtering][min_mapq]} \
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
          echo "# description: 'Non-redundant fraction and PCR-bottleneck coefficients from proper-pair, MAPQ-filtered nuclear fragments before duplicate removal.'"
          echo "# plot_type: 'table'"
          echo "# pconfig:"
          echo "#     id: 'library_complexity_table'"
          echo "#     namespace: 'ATAC'"
          echo -e "Sample\ttotal_frags\tdistinct\tM1_singletons\tM2_doubletons\tNRF\tPBC1\tPBC2"
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
        mapq=config["filtering"]["min_mapq"],
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
            --minMappingQuality {params.mapq} --skipZeros -p {threads} \
            --plotTitle "ATAC fingerprint (signal-to-noise)" 2> {log}
        """


# -----------------------------------------------------------------------------
# Outcome-blind sample similarity from genome-bin fragment signal. This stays
# separate from DESeq2 so inclusion decisions never require MA/volcano results.
# -----------------------------------------------------------------------------
rule qc_sample_similarity:
    input:
        bws=expand(f"{PROCESSED}/coverage/{{s}}.cpm.bw", s=SAMPLES),
    output:
        matrix=f"{PROCESSED}/qc/similarity/coverage_bins.npz",
        raw=f"{PROCESSED}/qc/similarity/coverage_bins.tsv",
        corr=f"{PROCESSED}/qc/similarity/spearman_correlation.tsv",
        corr_pdf=f"{RESULTS}/figures/qc_sample_correlation.pdf",
        corr_mqc=f"{PROCESSED}/qc/similarity/qc_sample_correlation_mqc.png",
        pca=f"{PROCESSED}/qc/similarity/pca_loadings.tsv",
        pca_pdf=f"{RESULTS}/figures/qc_sample_pca.pdf",
        pca_mqc=f"{PROCESSED}/qc/similarity/qc_sample_pca_mqc.png",
    params:
        labels=" ".join(SAMPLES),
        colours=lambda wc: shell_colours(SAMPLES, _COND_OF),
        cmap=SEQUENTIAL_CMAP,
    threads: config["resources"]["general_threads"]
    log:
        f"{LOGS}/qc/sample_similarity.log",
    conda:
        "../envs/deeptools.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.matrix}) $(dirname {output.corr_pdf})
        multiBigwigSummary bins -b {input.bws} --labels {params.labels} \
            --binSize 10000 -p {threads} -o {output.matrix} \
            --outRawCounts {output.raw} > {log} 2>&1
        plotCorrelation -in {output.matrix} --corMethod spearman \
            --whatToPlot heatmap --skipZeros --plotNumbers \
            --colorMap {params.cmap} --plotTitle "Outcome-blind sample correlation" \
            --plotFile {output.corr_pdf} --outFileCorMatrix {output.corr} \
            >> {log} 2>&1
        plotCorrelation -in {output.matrix} --corMethod spearman \
            --whatToPlot heatmap --skipZeros --plotNumbers \
            --colorMap {params.cmap} --plotTitle "Outcome-blind sample correlation" \
            --plotFile {output.corr_mqc} >> {log} 2>&1
        plotPCA -in {output.matrix} --labels {params.labels} --log2 \
            --ntop 50000 --colors {params.colours} \
            --plotTitle "Outcome-blind sample PCA" \
            --plotFile {output.pca_pdf} --outFileNameData {output.pca} \
            >> {log} 2>&1
        plotPCA -in {output.matrix} --labels {params.labels} --log2 \
            --ntop 50000 --colors {params.colours} \
            --plotTitle "Outcome-blind sample PCA" \
            --plotFile {output.pca_mqc} >> {log} 2>&1
        for artifact in {output.matrix} {output.raw} {output.corr} \
                        {output.corr_pdf} {output.corr_mqc} {output.pca} \
                        {output.pca_pdf} {output.pca_mqc}; do
            if [ ! -s "$artifact" ]; then
                echo "error: sample-similarity artifact is empty: $artifact" >> {log}
                exit 1
            fi
        done
        """


# -----------------------------------------------------------------------------
# Aggregate fragment-size / nucleosome overlay.
# -----------------------------------------------------------------------------
rule fragment_size_overlay:
    input:
        frag=expand(f"{PROCESSED}/qc/fragsize/{{s}}.fragsize.txt", s=SAMPLES),
        samples=config["samples"],
        scripts=script_inputs("plot_fragment_sizes.py", "palette.py"),
    output:
        png=f"{RESULTS}/figures/fragment_size_distribution.png",
        pdf=f"{RESULTS}/figures/fragment_size_distribution.pdf",
        mqc=f"{RESULTS}/figures/fragment_size_distribution_mqc.png",
        table=f"{PROCESSED}/qc/fragsize/fragment_partitions_mqc.tsv",
    log:
        f"{LOGS}/qc/fragment_overlay.log",
    conda:
        "../envs/deeptools.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.png}) $(dirname {output.table})
        python workflow/scripts/plot_fragment_sizes.py \
            --frag {input.frag} --samples {input.samples} \
            --out-plot {output.png} --out-table {output.table} 2> {log}
        for artifact in {output.png} {output.pdf} {output.mqc} {output.table}; do
            if [ ! -s "$artifact" ]; then
                echo "error: $artifact was not produced by plot_fragment_sizes.py" >> {log}
                exit 1
            fi
        done
        """


# -----------------------------------------------------------------------------
# Aggregate TSS enrichment: one matrix over all libraries -> metagene heatmap,
# overlay profile, and a numeric per-sample enrichment score.
# -----------------------------------------------------------------------------
rule tss_qc_regions:
    input:
        tss=f"{REF}/tss.bed",
        scripts=script_inputs("select_tss_regions.py"),
    output:
        bed=f"{PROCESSED}/qc/tss/qc_regions.bed",
    params:
        max_regions=config.get("qc", {}).get("tss", {}).get("max_regions", 0),
    log:
        f"{LOGS}/qc/tss_regions.log",
    conda:
        "../envs/deeptools.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.bed})
        python workflow/scripts/select_tss_regions.py \
            --input {input.tss} --output {output.bed} \
            --max-regions {params.max_regions} > {log} 2>&1
        """


rule tss_matrix_all:
    input:
        bws=expand(f"{PROCESSED}/coverage/{{s}}.cutsites.cpm.bw", s=SAMPLES),
        tss=f"{PROCESSED}/qc/tss/qc_regions.bed",
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
            -R {input.tss} -a 2000 -b 2000 --binSize 10 \
            --skipZeros -p {threads} \
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
        # Already shell-quoted: a bare '#2a78d6' would open a shell comment and
        # swallow both the flag's argument and the redirect that follows it.
        colours=lambda wc: shell_colours(SAMPLES, _COND_OF),
        cmap=SEQUENTIAL_CMAP,
    log:
        f"{LOGS}/qc/tss_metagene.log",
    conda:
        "../envs/deeptools.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.heatmap}) $(dirname {output.mqc})
        plotHeatmap -m {input.matrix} -o {output.heatmap} \
            --refPointLabel TSS --colorMap {params.cmap} --heatmapHeight 14 \
            --plotTitle "TSS enrichment" 2> {log}
        plotHeatmap -m {input.matrix} -o {output.mqc} \
            --refPointLabel TSS --colorMap {params.cmap} --heatmapHeight 14 \
            --plotTitle "TSS enrichment" 2>> {log}
        plotProfile -m {input.matrix} -o {output.profile} \
            --refPointLabel TSS --perGroup --colors {params.colours} \
            --plotTitle "TSS enrichment profile" 2>> {log}
        for f in {output.heatmap} {output.profile} {output.mqc}; do
            if [ ! -s "$f" ]; then
                echo "error: deepTools exited 0 but $f is empty" >> {log}
                exit 1
            fi
        done
        """


rule tss_enrichment_score:
    input:
        matrix=f"{PROCESSED}/qc/tss/tss_all.matrix.gz",
        scripts=script_inputs("tss_score.py"),
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
# One auditable row per library: flow-through counts, assay QC, and review flag.
# This table never auto-excludes data. Review it, then set include=false in the
# sample sheet and rerun so exclusion is explicit and provenance-backed.
# -----------------------------------------------------------------------------
rule qc_summary_sample:
    input:
        raw=f"{PROCESSED}/aligned/{{sample}}.sorted.bam",
        filtered=f"{PROCESSED}/filtered/{{sample}}.filtered.bam",
        peaks=f"{PROCESSED}/peaks/macs3/{{sample}}_peaks.narrowPeak",
        frip=f"{PROCESSED}/qc/frip/{{sample}}.frip.txt",
        complexity=(f"{PROCESSED}/qc/complexity/{{sample}}.complexity.txt"
                    if config.get("qc", {}).get("library_complexity", True) else []),
        tss=(f"{PROCESSED}/qc/tss/tss_enrichment_mqc.tsv"
             if config.get("qc", {}).get("tss", {}).get("enabled", True) else []),
        scripts=script_inputs("qc_summary.py"),
    output:
        f"{PROCESSED}/qc/summary/{{sample}}.qc.tsv",
    params:
        complexity_arg=(f"--complexity {PROCESSED}/qc/complexity/{{sample}}.complexity.txt"
                        if config.get("qc", {}).get("library_complexity", True) else ""),
        tss_arg=(f"--tss {PROCESSED}/qc/tss/tss_enrichment_mqc.tsv"
                 if config.get("qc", {}).get("tss", {}).get("enabled", True) else ""),
        build=config.get("genome", {}).get("build", "custom"),
        mito_args=mitochondrial_contig_args("--mitochondrial-contig"),
    log:
        f"{LOGS}/qc/summary_{{sample}}.log",
    conda:
        "../envs/qc.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output})
        python workflow/scripts/qc_summary.py \
          --sample {wildcards.sample} --raw-bam {input.raw} \
          --filtered-bam {input.filtered} --peaks {input.peaks} \
          --frip {input.frip} {params.complexity_arg} {params.tss_arg} \
          --genome-build {params.build} {params.mito_args} \
          --out {output} 2> {log}
        """


rule qc_summary_table:
    input:
        expand(f"{PROCESSED}/qc/summary/{{s}}.qc.tsv", s=SAMPLES),
    output:
        f"{PROCESSED}/qc/summary/qc_decisions_mqc.tsv",
    log:
        f"{LOGS}/qc/qc_summary_table.log",
    shell:
        r"""
        mkdir -p $(dirname {output})
        {{
          echo "# id: 'atac_qc_decisions'"
          echo "# section_name: 'ATAC library QC review'"
          echo "# description: 'Transparent review flags; thresholds are guides, not automatic exclusions. Set include=false in the sample sheet after review.'"
          echo "# plot_type: 'table'"
          head -n 1 $(echo {input} | cut -d' ' -f1)
          for f in {input}; do tail -n 1 "$f"; done
        }} > {output} 2> {log}
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
    if qc.get("sample_similarity", True):
        out += [
            f"{RESULTS}/figures/qc_sample_correlation.pdf",
            f"{RESULTS}/figures/qc_sample_pca.pdf",
        ]
    if config.get("report", {}).get("enabled", True):
        out += [
            f"{RESULTS}/figures/fragment_size_distribution.png",
            f"{RESULTS}/figures/fragment_size_distribution.pdf",
        ]
        if qc.get("tss", {}).get("enabled", True):
            out += [
                f"{RESULTS}/figures/tss_enrichment_heatmap.pdf",
                f"{RESULTS}/figures/tss_enrichment_profile.pdf",
                f"{PROCESSED}/qc/tss/tss_enrichment_mqc.tsv",
            ]
    out += [f"{PROCESSED}/qc/summary/qc_decisions_mqc.tsv"]
    return out
