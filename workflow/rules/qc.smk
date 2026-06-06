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
        html=f"{RESULTS}/qc/fastqc/{{sample}}_{{read}}_fastqc.html",
        zip=f"{RESULTS}/qc/fastqc/{{sample}}_{{read}}_fastqc.zip",
    params:
        outdir=f"{RESULTS}/qc/fastqc",
    log:
        f"{LOGS}/qc/fastqc_{{sample}}_{{read}}.log",
    conda:
        "../envs/qc.yaml"
    shell:
        r"""
        mkdir -p {params.outdir}
        fastqc -o {params.outdir} {input} 2> {log}
        # Normalise FastQC's auto-naming to <sample>_<read>_fastqc.*
        base=$(basename {input} | sed -E 's/\.(fastq|fq)(\.gz)?$//')
        for ext in html zip; do
            if [ -f "{params.outdir}/${{base}}_fastqc.$ext" ]; then
                mv -f "{params.outdir}/${{base}}_fastqc.$ext" \
                      "{params.outdir}/{wildcards.sample}_{wildcards.read}_fastqc.$ext"
            fi
        done
        """


rule fragment_sizes:
    input:
        bam=f"{RESULTS}/filtered/{{sample}}.filtered.bam",
        bai=f"{RESULTS}/filtered/{{sample}}.filtered.bam.bai",
    output:
        hist=f"{RESULTS}/qc/fragsize/{{sample}}.fragsize.txt",
        plot=f"{RESULTS}/qc/fragsize/{{sample}}.fragsize.png",
    threads: config["resources"]["general_threads"]
    log:
        f"{LOGS}/qc/fragsize_{{sample}}.log",
    conda:
        "../envs/deeptools.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.hist})
        bamPEFragmentSize -b {input.bam} -p {threads} \
            --outRawFragmentLengths {output.hist} \
            -o {output.plot} 2> {log}
        """


rule tss_enrichment:
    input:
        bw=f"{RESULTS}/coverage/{{sample}}.cpm.bw",
        tss=f"{REF}/tss.bed",
    output:
        mat=f"{RESULTS}/qc/tss/{{sample}}.tss_matrix.gz",
        plot=f"{RESULTS}/qc/tss/{{sample}}.tss_enrichment.png",
    threads: config["resources"]["general_threads"]
    log:
        f"{LOGS}/qc/tss_{{sample}}.log",
    conda:
        "../envs/deeptools.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.mat})
        computeMatrix reference-point --referencePoint center \
            -S {input.bw} -R {input.tss} \
            -a 2000 -b 2000 -p {threads} \
            --skipZeros -o {output.mat} 2> {log}
        plotProfile -m {output.mat} -o {output.plot} \
            --refPointLabel TSS 2>> {log}
        """


rule frip:
    input:
        bam=f"{RESULTS}/filtered/{{sample}}.filtered.bam",
        peaks=f"{RESULTS}/peaks/macs3/{{sample}}_peaks.narrowPeak",
    output:
        f"{RESULTS}/qc/frip/{{sample}}.frip.txt",
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
        expand(f"{RESULTS}/qc/frip/{{s}}.frip.txt", s=SAMPLES),
    output:
        f"{RESULTS}/qc/frip/frip_mqc.tsv",
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
