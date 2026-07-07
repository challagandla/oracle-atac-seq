# =============================================================================
# counts.smk — count fragments in consensus peaks (featureCounts)
# =============================================================================
# featureCounts (Subread) in paired-end fragment-counting mode produces the
# peak-by-sample matrix that feeds DESeq2 and chromVAR.

rule featurecounts:
    input:
        saf=f"{RESULTS}/consensus/consensus_peaks.saf",
        bams=expand(f"{PROCESSED}/filtered/{{s}}.filtered.bam", s=SAMPLES),
    output:
        counts=f"{RESULTS}/counts/consensus_counts.tsv",
        summary=f"{RESULTS}/counts/consensus_counts.tsv.summary",
    threads: config["resources"]["general_threads"]
    log:
        f"{LOGS}/counts/featurecounts.log",
    conda:
        "../envs/peaks.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.counts})
        if [ "$(tail -n +2 {input.saf} | wc -l)" -eq 0 ]; then
            {{
              echo '# Program:featureCounts; Command:"featureCounts skipped because consensus SAF has no features"'
              printf "Geneid\tChr\tStart\tEnd\tStrand\tLength"
              for bam in {input.bams}; do printf "\t%s" "$bam"; done
              printf "\n"
            }} > {output.counts}
            {{
              printf "Status"
              for bam in {input.bams}; do printf "\t%s" "$bam"; done
              printf "\nAssigned"
              for bam in {input.bams}; do printf "\t0"; done
              printf "\nUnassigned_NoFeatures"
              for bam in {input.bams}; do
                  n=$(samtools view -c -f 2 -F 1804 "$bam")
                  printf "\t%s" "$n"
              done
              printf "\n"
            }} > {output.summary}
            echo "featureCounts skipped: consensus SAF has no features." > {log}
        else
            featureCounts -p --countReadPairs -F SAF \
                -a {input.saf} -o {output.counts} \
                -T {threads} {input.bams} 2> {log}
        fi
        """
