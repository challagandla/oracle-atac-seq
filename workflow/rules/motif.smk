# =============================================================================
# motif.smk — HOMER motif enrichment on differentially-accessible peaks
# =============================================================================
# findMotifsGenome.pl runs de-novo + known motif enrichment on the up- and
# down-regulated peak sets from DESeq2.
#
# We hand HOMER the genome FASTA rather than a HOMER genome *package*. Two
# reasons, both practical:
#   * the package is a separate ~3 GB `configureHomer.pl -install hg38` step that
#     conda cannot express, so a fresh checkout could never run this rule;
#   * the packaged genomes are UCSC-named. Our peaks carry the alignment genome's
#     names (Ensembl "1"). Passing the same FASTA the reads were aligned to means
#     the names agree by construction.
# HOMER writes a preparsed index beside the peaks; give it its own directory so
# the two directions do not race.

rule homer_motifs:
    input:
        peaks=f"{RESULTS}/diffacc/{{direction}}_peaks.bed",
        tested=f"{RESULTS}/diffacc/tested_peaks.bed",
        fasta=genome_fasta(),
    output:
        html=f"{RESULTS}/motif/{{direction}}/homerResults.html",
        known=f"{RESULTS}/motif/{{direction}}/knownResults.txt",
    params:
        size=config["motif"]["homer_size"],
        outdir=f"{RESULTS}/motif/{{direction}}",
        preparsed=f"{RESULTS}/motif/_preparsed_{{direction}}",
    threads: config["resources"]["general_threads"]
    log:
        f"{LOGS}/motif/homer_{{direction}}.log",
    conda:
        "../envs/motif.yaml"
    shell:
        r"""
        # HOMER creates many files that Snakemake cannot enumerate. Recreate the
        # direction-specific directory on every run so an empty or changed
        # contrast cannot leave apparently current files from an older result.
        rm -rf {params.outdir} {params.preparsed}
        mkdir -p {params.outdir} {params.preparsed}
        # An empty peak set here means DESeq2 found nothing in this direction.
        # That is a legitimate result (and the calibration figure reports it), so
        # emit an honest placeholder rather than failing the DAG -- but say the
        # count, so "no motifs" can never be mistaken for "motif analysis ran".
        n=$(wc -l < {input.peaks})
        if [ "$n" -eq 0 ]; then
            echo "No {wildcards.direction}-regulated peaks; HOMER not run." > {log}
            printf '<html><body><h3>No %s-regulated peaks at the configured FDR.</h3>\n<p>HOMER was not run. This is a result, not a failure.</p></body></html>\n' \
                "{wildcards.direction}" > {output.html}
            printf 'Motif Name\tConsensus\tP-value\tq-value (Benjamini)\n' > {output.known}
        else
            echo "HOMER on $n {wildcards.direction} peaks" > {log}
            # HOMER needs a 6-column BED (name + score + strand).
            awk 'BEGIN{{OFS="\t"}}{{print $1,$2,$3,"peak"NR,0,"+"}}' {input.peaks} \
                > {params.outdir}/input.bed
            # Compare against the other tested accessible loci, not arbitrary
            # genomic sequence. This reduces open-chromatin composition bias.
            awk 'BEGIN{{OFS="\t"}} NR==FNR{{hit[$4]=1; next}} \
                 !($4 in hit){{print $1,$2,$3,$4,0,"+"}}' \
                 {input.peaks} {input.tested} > {params.outdir}/background.bed
            if [ ! -s {params.outdir}/background.bed ]; then
                echo "No non-target tested peaks are available for motif background." >> {log}
                exit 1
            fi
            findMotifsGenome.pl {params.outdir}/input.bed {input.fasta} \
                {params.outdir} -size {params.size} -p {threads} \
                -bg {params.outdir}/background.bed \
                -preparsedDir {params.preparsed} >> {log} 2>&1
        fi
        """
