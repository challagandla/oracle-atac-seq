# =============================================================================
# qc.smk — read QC, fragment-size distribution, TSS enrichment, FRiP
# =============================================================================
# ATAC-specific QC follows ENCODE standards and Yan et al. (2020):
#   - FastQC on raw reads
#   - fragment-size distribution (expect nucleosome-free + mono/di peaks)
#   - TSS enrichment score (computeMatrix on the TSS BED)
#   - FRiP (fraction of reads in peaks)
#   - flagstat before/after filtering

rule fastqc_raw:
    input:
        lambda wc: raw_fastqs(wc.sample)[0 if wc.read == "R1" else 1],
    output:
        html=f"{PROCESSED}/qc/fastqc/{{sample}}_{{read}}_fastqc.html",
        zip=f"{PROCESSED}/qc/fastqc/{{sample}}_{{read}}_fastqc.zip",
    params:
        outdir=f"{PROCESSED}/qc/fastqc",
    log:
        f"{LOGS}/qc/fastqc_{{sample}}_{{read}}.log",
    # FastQC needs ~10GB heap on very deep libraries (see --memory below). Declare
    # it so `snakemake --resources mem_mb=<node_total>` won't pack too many of
    # these JVMs onto one node (running many concurrently exhausts RAM and they
    # die silently at ~95%).
    resources:
        mem_mb=10000,
    conda:
        "../envs/qc.yaml"
    shell:
        r"""
        mkdir -p {params.outdir}
        # --memory 10000: FastQC's 512MB default makes its per-sequence
        # duplication module die near 95% on very deep (>150M-read) libraries;
        # ~10GB clears it. The JVM reserves but does not commit, so this is
        # harmless for small inputs.
        fastqc --memory 10000 -o {params.outdir} {input} 2> {log}
        # Normalise FastQC's auto-naming to <sample>_<read>_fastqc.*
        base=$(basename {input} | sed -E 's/\.(fastq|fq)(\.gz)?$//')
        for ext in html zip; do
            src="{params.outdir}/${{base}}_fastqc.$ext"
            dst="{params.outdir}/{wildcards.sample}_{wildcards.read}_fastqc.$ext"
            if [ -f "$src" ] && [ "$src" != "$dst" ]; then
                mv -f "{params.outdir}/${{base}}_fastqc.$ext" \
                      "$dst"
            fi
        done
        """


rule fragment_sizes:
    input:
        bam=f"{PROCESSED}/filtered/{{sample}}.filtered.bam",
        bai=f"{PROCESSED}/filtered/{{sample}}.filtered.bam.bai",
    output:
        hist=f"{PROCESSED}/qc/fragsize/{{sample}}.fragsize.txt",
        plot=f"{PROCESSED}/qc/fragsize/{{sample}}.fragsize.png",
    threads: config["resources"]["general_threads"]
    params:
        enabled=config.get("qc", {}).get("fragment_sizes", {}).get("enabled", True),
    log:
        f"{LOGS}/qc/fragsize_{{sample}}.log",
    conda:
        "../envs/deeptools.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.hist})
        if [ "{params.enabled}" = "True" ]; then
            bamPEFragmentSize -b {input.bam} -p {threads} \
                --outRawFragmentLengths {output.hist} \
                -o {output.plot} 2> {log}
        else
            printf "#bamPEFragmentSize\nSize\tOccurrences\tSample\n" > {output.hist}
            printf "Fragment-size QC disabled for this run.\n" > {log}
            printf "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mP8/x8AAwMCAO+/p9sAAAAASUVORK5CYII=\n" \
                | base64 -d > {output.plot}
        fi
        """


rule tss_enrichment:
    input:
        bw=f"{PROCESSED}/coverage/{{sample}}.cpm.bw",
        tss=f"{REF}/tss.bed",
    output:
        mat=f"{PROCESSED}/qc/tss/{{sample}}.tss_matrix.gz",
        plot=f"{PROCESSED}/qc/tss/{{sample}}.tss_enrichment.png",
    threads: config["resources"]["general_threads"]
    params:
        enabled=config.get("qc", {}).get("tss", {}).get("enabled", True),
        max_regions=config.get("qc", {}).get("tss", {}).get("max_regions", 0),
    log:
        f"{LOGS}/qc/tss_{{sample}}.log",
    conda:
        "../envs/deeptools.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.mat})
        if [ "{params.enabled}" = "True" ]; then
            regions="{input.tss}"
            if [ "{params.max_regions}" != "0" ]; then
                regions="{output.mat}.regions.bed"
                head -n {params.max_regions} {input.tss} > "$regions"
            fi
            computeMatrix reference-point --referencePoint center \
                -S {input.bw} -R "$regions" \
                -a 2000 -b 2000 -p {threads} \
                --skipZeros -o {output.mat} 2> {log}
            plotProfile -m {output.mat} -o {output.plot} \
                --refPointLabel TSS 2>> {log}
            if [ "$regions" != "{input.tss}" ]; then rm -f "$regions"; fi
        else
            printf "TSS enrichment QC disabled for this run.\n" | gzip -c > {output.mat}
            printf "TSS enrichment QC disabled for this run.\n" > {log}
            printf "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mP8/x8AAwMCAO+/p9sAAAAASUVORK5CYII=\n" \
                | base64 -d > {output.plot}
        fi
        """


rule frip:
    input:
        bam=f"{PROCESSED}/filtered/{{sample}}.filtered.bam",
        peaks=f"{PROCESSED}/peaks/macs3/{{sample}}_peaks.narrowPeak",
    output:
        f"{PROCESSED}/qc/frip/{{sample}}.frip.txt",
    log:
        f"{LOGS}/qc/frip_{{sample}}.log",
    conda:
        "../envs/peaks.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output})
        python workflow/scripts/frip.py \
            --bam {input.bam} --peaks {input.peaks} \
            --sample {wildcards.sample} --out {output} 2> {log}
        """


# Aggregate per-sample FRiP into a MultiQC custom-content table. The `_mqc.tsv`
# suffix and the embedded `# id:` header make MultiQC render it as its own
# section automatically (no extra config needed).
rule frip_table:
    input:
        expand(f"{PROCESSED}/qc/frip/{{s}}.frip.txt", s=SAMPLES),
    output:
        f"{PROCESSED}/qc/frip/frip_mqc.tsv",
    log:
        f"{LOGS}/qc/frip_table.log",
    shell:
        r"""
        {{
          echo "# id: 'frip'"
          echo "# section_name: 'FRiP (Fraction of Reads in Peaks)'"
          echo "# description: 'Fraction of properly-paired reads overlapping MACS3 peaks. ENCODE recommends > 0.2-0.3.'"
          echo "# plot_type: 'table'"
          echo "# pconfig:"
          echo "#     id: 'frip_table'"
          echo "#     namespace: 'ATAC'"
          echo -e "Sample\ttotal_reads\treads_in_peaks\tFRiP"
          # concatenate the data row (line 2) of every per-sample file
          for f in {input}; do tail -n +2 "$f"; done
        }} > {output} 2> {log}
        """
