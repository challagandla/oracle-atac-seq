# =============================================================================
# consensus.smk — reproducible consensus peak set across samples
# =============================================================================
# Strategy (Yan et al. 2020; ENCODE): merge MACS3 peaks from all samples, then
# keep merged intervals supported by >= consensus_min_overlap samples. The
# result is a clean SAF/BED used for quantification and differential testing.

rule consensus_peaks:
    input:
        peaks=expand(
            f"{PROCESSED}/peaks/macs3/{{s}}_peaks.narrowPeak", s=SAMPLES
        ),
        chrom=f"{REF}/chrom.sizes",
    output:
        bed=f"{RESULTS}/consensus/consensus_peaks.bed",
        saf=f"{RESULTS}/consensus/consensus_peaks.saf",
    params:
        min_overlap=config["peaks"]["consensus_min_overlap"],
    log:
        f"{LOGS}/consensus/consensus.log",
    conda:
        "../envs/peaks.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.bed})
        python workflow/scripts/make_consensus.py \
            --peaks {input.peaks} \
            --chrom {input.chrom} \
            --min-overlap {params.min_overlap} \
            --bed {output.bed} --saf {output.saf} 2> {log}
        """
