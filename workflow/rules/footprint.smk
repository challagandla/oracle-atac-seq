# =============================================================================
# footprint.smk — TF footprinting with TOBIAS (optional)
# =============================================================================
# TOBIAS (Bentsen et al. 2020, Nat Commun) corrects Tn5 sequence bias
# (ATACorrect), scores footprints (FootprintScores), and detects bound vs
# unbound motifs (BINDetect) per condition. Runs only when footprinting.enabled.
# Requires a motif database (JASPAR/MEME) set in config: footprinting.motif_db.

rule tobias_merge_condition:
    input:
        bams=lambda wc: expand(
            f"{PROCESSED}/shifted/{{s}}.shifted.bam",
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
        fasta=genome_fasta(),
        peaks=f"{RESULTS}/consensus/consensus_peaks.bed",
    output:
        corrected=f"{RESULTS}/footprint/{{cond}}_corrected.bw",
    params:
        outdir=f"{RESULTS}/footprint",
        prefix=lambda wc: wc.cond,
    threads: config["resources"]["general_threads"]
    log:
        f"{LOGS}/footprint/atacorrect_{{cond}}.log",
    conda:
        "../envs/footprint.yaml"
    shell:
        r"""
        TOBIAS ATACorrect --bam {input.bam} --genome {input.fasta} \
            --peaks {input.peaks} --outdir {params.outdir} \
            --prefix {params.prefix} --cores {threads} 2> {log}
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
        scores=f"{RESULTS}/footprint/{{cond}}_footprints.bw",
        fasta=genome_fasta(),
        peaks=f"{RESULTS}/consensus/consensus_peaks.bed",
    output:
        bed=f"{RESULTS}/footprint/{{cond}}_footprints.bed",
    params:
        motifs=config["footprinting"]["motif_db"],
        outdir=f"{RESULTS}/footprint/{{cond}}_bindetect",
    threads: config["resources"]["general_threads"]
    log:
        f"{LOGS}/footprint/bindetect_{{cond}}.log",
    conda:
        "../envs/footprint.yaml"
    shell:
        r"""
        TOBIAS BINDetect --motifs {params.motifs} \
            --signals {input.scores} --genome {input.fasta} \
            --peaks {input.peaks} --outdir {params.outdir} \
            --cores {threads} 2> {log}
        # Flatten the per-TF bound sites into a single BED for the target list.
        cat {params.outdir}/*/beds/*_bound.bed > {output.bed} 2>> {log} || \
            : > {output.bed}
        """
