# =============================================================================
# footprint.smk — TF footprinting with TOBIAS (optional)
# =============================================================================
# TOBIAS (Bentsen et al. 2020, Nat Commun) corrects Tn5 sequence bias
# (ATACorrect), scores footprints (FootprintScores), and detects bound vs
# unbound motifs (BINDetect) across all conditions in one comparative call.
# Runs only when footprinting.enabled.
# Requires a motif database (JASPAR/MEME) set in config: footprinting.motif_db.

rule tobias_merge_condition:
    input:
        bams=lambda wc: expand(
            f"{PROCESSED}/filtered/{{s}}.filtered.bam",
            s=samples_in_condition(wc.cond),
        ),
    output:
        bam=temp(f"{RESULTS}/footprint/{{cond}}.merged.bam"),
        bai=temp(f"{RESULTS}/footprint/{{cond}}.merged.bam.bai"),
    threads: config["resources"]["general_threads"]
    log:
        f"{LOGS}/footprint/merge_{{cond}}.log",
    conda:
        "../envs/align.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.bam})
        samtools merge -f -@ {threads} {output.bam} {input.bams} 2> {log}
        samtools index {output.bam}
        """


rule tobias_atacorrect:
    input:
        bam=f"{RESULTS}/footprint/{{cond}}.merged.bam",
        bai=f"{RESULTS}/footprint/{{cond}}.merged.bam.bai",
        fasta=genome_fasta(),
        peaks=f"{RESULTS}/consensus/consensus_peaks.bed",
        blacklist=blacklist_bed() if blacklist_bed() else [],
    output:
        corrected=f"{RESULTS}/footprint/{{cond}}_corrected.bw",
        uncorrected=f"{RESULTS}/footprint/{{cond}}_uncorrected.bw",
        bias=f"{RESULTS}/footprint/{{cond}}_bias.bw",
        expected=f"{RESULTS}/footprint/{{cond}}_expected.bw",
        diagnostic=f"{RESULTS}/footprint/{{cond}}_atacorrect.pdf",
    params:
        outdir=f"{RESULTS}/footprint",
        prefix=lambda wc: wc.cond,
        blacklist_arg=(f"--blacklist {blacklist_bed()}" if blacklist_bed() else ""),
    threads: config["resources"]["general_threads"]
    log:
        f"{LOGS}/footprint/atacorrect_{{cond}}.log",
    conda:
        "../envs/footprint.yaml"
    shell:
        r"""
        TOBIAS ATACorrect --bam {input.bam} --genome {input.fasta} \
            --peaks {input.peaks} --outdir {params.outdir} \
            --prefix {params.prefix} {params.blacklist_arg} \
            --cores {threads} 2> {log}
        for artifact in {output.corrected} {output.uncorrected} {output.bias} \
                        {output.expected} {output.diagnostic}; do
            if [ ! -s "$artifact" ]; then
                echo "error: TOBIAS ATACorrect did not produce $artifact" >> {log}
                exit 1
            fi
        done
        """


rule tobias_footprint_scores:
    input:
        corrected=f"{RESULTS}/footprint/{{cond}}_corrected.bw",
        peaks=f"{RESULTS}/consensus/consensus_peaks.bed",
    output:
        scores=f"{RESULTS}/footprint/{{cond}}_footprints.bw",
    threads: config["resources"]["general_threads"]
    log:
        f"{LOGS}/footprint/scores_{{cond}}.log",
    conda:
        "../envs/footprint.yaml"
    shell:
        r"""
        TOBIAS FootprintScores --signal {input.corrected} \
            --regions {input.peaks} --output {output.scores} \
            --cores {threads} 2> {log}
        """


rule tobias_bindetect:
    input:
        scores=expand(f"{RESULTS}/footprint/{{c}}_footprints.bw", c=conditions()),
        fasta=genome_fasta(),
        peaks=f"{RESULTS}/consensus/consensus_peaks.bed",
        motifs=(config["footprinting"]["motif_db"]
                if config["footprinting"].get("enabled", False) else []),
    output:
        results=f"{RESULTS}/footprint/bindetect/bindetect_results.txt",
        figures=f"{RESULTS}/footprint/bindetect/bindetect_figures.pdf",
    params:
        outdir=f"{RESULTS}/footprint/bindetect",
        cond_names=" ".join(conditions()),
    threads: config["resources"]["general_threads"]
    log:
        f"{LOGS}/footprint/bindetect.log",
    conda:
        "../envs/footprint.yaml"
    shell:
        r"""
        # BINDetect writes motif-specific files beyond the two summary outputs.
        # Recreate the directory so removed motifs cannot leave stale results.
        rm -rf {params.outdir}
        mkdir -p {params.outdir}
        TOBIAS BINDetect --motifs {input.motifs} \
            --signals {input.scores} --genome {input.fasta} \
            --peaks {input.peaks} --outdir {params.outdir} \
            --cond_names {params.cond_names} \
            --cores {threads} 2> {log}
        test -s {output.results}
        test -s {output.figures}
        """
