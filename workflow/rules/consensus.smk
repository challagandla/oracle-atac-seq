# =============================================================================
# consensus.smk — reproducible consensus peak set across samples
# =============================================================================
# Fixed-width summit windows are retained when at least N distinct replicates
# from one condition support the locus. This preserves condition-specific peaks
# while preventing broad intervals created by transitive overlap chains.

rule consensus_peaks:
    input:
        peaks=expand(
            f"{PROCESSED}/peaks/macs3/{{s}}_peaks.narrowPeak", s=SAMPLES
        ),
        chrom=f"{REF}/chrom.sizes",
        blacklist=(blacklist_bed()
                   if config["filtering"]["remove_blacklist"] and blacklist_bed()
                   else []),
        scripts=script_inputs("make_consensus.py"),
    output:
        bed=f"{RESULTS}/consensus/consensus_peaks.bed",
        saf=f"{RESULTS}/consensus/consensus_peaks.saf",
    params:
        min_replicates=config["peaks"]["consensus_min_replicates"],
        width=config["peaks"]["consensus_peak_width"],
        conditions=" ".join(samples.loc[SAMPLES, "condition"]),
        blacklist_arg=(f"--blacklist {blacklist_bed()}"
                       if config["filtering"]["remove_blacklist"] and blacklist_bed()
                       else ""),
    log:
        f"{LOGS}/consensus/consensus.log",
    conda:
        "../envs/peaks.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.bed})
        python workflow/scripts/make_consensus.py \
            --peaks {input.peaks} \
            --conditions {params.conditions} \
            --chrom {input.chrom} \
            --min-replicates {params.min_replicates} \
            --width {params.width} \
            {params.blacklist_arg} \
            --bed {output.bed} --saf {output.saf} 2> {log}

        # An empty consensus set silently turns every downstream stage into a
        # no-op: featureCounts counts nothing, DESeq2 tests nothing, and the
        # figures are drawn from an empty frame. Stop here instead.
        n=$(wc -l < {output.bed})
        echo "consensus: $n fixed-width peaks (>= {params.min_replicates} replicate(s) within a condition)" >> {log}
        if [ "$n" -eq 0 ]; then
            echo "error: consensus peak set is empty." >> {log}
            echo "       No locus was supported by >= {params.min_replicates} replicates within a condition." >> {log}
            echo "       Review per-sample peak calls and biological replication." >> {log}
            exit 1
        fi
        if [ "{params.blacklist_arg}" != "" ]; then
            overlap=$(bedtools intersect -u -a {output.bed} -b {input.blacklist} | wc -l)
            if [ "$overlap" -ne 0 ]; then
                echo "error: consensus blacklist postcondition failed ($overlap overlaps)." >> {log}
                exit 1
            fi
        fi
        """
