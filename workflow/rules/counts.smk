# =============================================================================
# counts.smk — count fragments in consensus peaks (featureCounts)
# =============================================================================
# Count one BAM per Rsubread process, then validate and merge the matrices.
# Isolating libraries makes this stage restartable per sample and prevents one
# failed library from leaving an unchecked partial cohort matrix.

COUNT_SAMPLE_DIR = f"{RESULTS}/counts/per_sample"


rule featurecounts_sample:
    input:
        saf=f"{RESULTS}/consensus/consensus_peaks.saf",
        bam=f"{PROCESSED}/filtered/{{sample}}.filtered.bam",
        script=script_inputs("run_featurecounts.R"),
    output:
        counts=f"{COUNT_SAMPLE_DIR}/{{sample}}.counts.tsv",
        summary=f"{COUNT_SAMPLE_DIR}/{{sample}}.counts.tsv.summary",
    threads: config["resources"]["general_threads"]
    log:
        f"{LOGS}/counts/{{sample}}.featurecounts.log",
    conda:
        "../envs/counts.yaml"
    shell:
        r"""
        mkdir -p "$(dirname {output.counts:q})" "$(dirname {log:q})"
        # An empty feature universe is a failed peak-construction stage, not a
        # valid zero-row matrix for differential analysis.
        if [ "$(tail -n +2 {input.saf:q} | wc -l)" -eq 0 ]; then
            echo "error: {input.saf} has no features; consensus peak calling failed." > {log:q}
            exit 1
        fi

        # Rsubread creates scratch files while pairing coordinate-sorted reads.
        # Give every concurrent library its own directory and clean it on exit.
        tmpdir=$(mktemp -d "${{TMPDIR:-/tmp}}/oracle-atac-featurecounts.{wildcards.sample}.XXXXXX")
        trap 'rm -rf "$tmpdir"' EXIT

        export R_LIBS_USER="$CONDA_PREFIX/lib/R/library"
        Rscript --vanilla {input.script:q} \
            {wildcards.sample:q} {input.saf:q} {input.bam:q} {threads} "$tmpdir" \
            {output.counts:q} {output.summary:q} > {log:q} 2>&1
        """


rule merge_featurecounts:
    input:
        counts=expand(f"{COUNT_SAMPLE_DIR}/{{s}}.counts.tsv", s=SAMPLES),
        summaries=expand(f"{COUNT_SAMPLE_DIR}/{{s}}.counts.tsv.summary", s=SAMPLES),
        script=script_inputs("merge_featurecounts.py"),
    output:
        counts=f"{RESULTS}/counts/consensus_counts.tsv",
        summary=f"{RESULTS}/counts/consensus_counts.tsv.summary",
    params:
        samples=" ".join(shlex.quote(sample) for sample in SAMPLES),
    log:
        f"{LOGS}/counts/merge_featurecounts.log",
    conda:
        "../envs/counts.yaml"
    shell:
        r"""
        mkdir -p "$(dirname {output.counts:q})" "$(dirname {log:q})"
        python {input.script:q} \
            --counts {input.counts:q} \
            --summaries {input.summaries:q} \
            --samples {params.samples} \
            --out-counts {output.counts:q} \
            --out-summary {output.summary:q} 2> {log:q}
        """
