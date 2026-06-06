# =============================================================================
# motif.smk — HOMER motif enrichment on differentially-accessible peaks
# =============================================================================
# findMotifsGenome.pl runs de-novo + known motif enrichment on the up- and
# down-regulated peak sets from DESeq2, against the matched HOMER genome.
# NOTE: the HOMER genome package must be installed once:
#   perl $(dirname $(which homer))/../share/homer*/configureHomer.pl -install <genome>

rule homer_motifs:
    input:
        peaks=f"{RESULTS}/diffacc/{{direction}}_peaks.bed",
    output:
        html=f"{RESULTS}/motif/{{direction}}/homerResults.html",
    params:
        genome=GENOME.get("homer_genome", ""),
        size=config["motif"]["homer_size"],
        outdir=f"{RESULTS}/motif/{{direction}}",
    threads: config["resources"]["general_threads"]
    log:
        f"{LOGS}/motif/homer_{{direction}}.log",
    conda:
        "../envs/motif.yaml"
    shell:
        r"""
        mkdir -p {params.outdir}
        # HOMER needs a 6-column BED (name + score + strand); add if missing.
        awk 'BEGIN{{OFS="\t"}}{{print $1,$2,$3,"peak"NR,0,"+"}}' {input.peaks} \
            > {params.outdir}/input.bed
        findMotifsGenome.pl {params.outdir}/input.bed {params.genome} \
            {params.outdir} -size {params.size} -p {threads} 2> {log}
        """
