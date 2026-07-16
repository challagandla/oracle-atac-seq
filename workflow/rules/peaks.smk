# =============================================================================
# peaks.smk — peak calling with MACS3 (per sample) and Genrich (per condition)
# =============================================================================
# MACS3 receives the unshifted, filtered paired-end BAM in BAMPE mode and models
# observed fragments directly. Genrich is an optional, condition-pooled ATAC
# cross-check; it is not a substitute for IDR.

rule macs3_callpeak:
    input:
        bam=f"{PROCESSED}/filtered/{{sample}}.filtered.bam",
        scripts=script_inputs("normalize_narrowpeak.py"),
    output:
        narrowpeak=f"{PROCESSED}/peaks/macs3/{{sample}}_peaks.narrowPeak",
        xls=f"{PROCESSED}/peaks/macs3/{{sample}}_peaks.xls",
        summits=f"{PROCESSED}/peaks/macs3/{{sample}}_summits.bed",
    params:
        gsize=macs_gsize(),
        q=config["peaks"]["macs3_qvalue"],
        extra=config["peaks"]["macs3_extra"],
        outdir=f"{PROCESSED}/peaks/macs3",
        name=lambda wc: wc.sample,
    log:
        f"{LOGS}/peaks/macs3_{{sample}}.log",
    conda:
        "../envs/peaks.yaml"
    shell:
        r"""
        mkdir -p {params.outdir}
        macs3 callpeak -t {input.bam} -f BAMPE \
            -g {params.gsize} -q {params.q} {params.extra} \
            -n {params.name} --outdir {params.outdir} 2> {log}

        # MACS3 writes an empty narrowPeak and exits 0 when it finds nothing, so
        # the DAG goes green on a run that produced no peaks. An ATAC library with
        # zero peaks is a failed library or an under-powered subsample, never a
        # result. Say which.
        n=$(wc -l < {output.narrowpeak})
        reads=$(samtools view -c {input.bam})
        echo "macs3: $n peaks from $reads reads" >> {log}
        if [ "$n" -eq 0 ]; then
            echo "error: MACS3 called 0 peaks for {wildcards.sample} from $reads reads." >> {log}
            echo "       Too few reads to model a background, or q-value too strict" >> {log}
            echo "       (peaks.macs3_qvalue = {params.q}). A real ATAC library needs" >> {log}
            echo "       enough usable fragments for stable peak detection." >> {log}
            exit 1
        fi
        # UCSC narrowPeak column 5 is an integer display score in [0, 1000].
        # MACS3 can emit larger values for very significant peaks, so validate
        # every field and clamp only that display column; signal/p/q values stay
        # unchanged for ranking and downstream analysis.
        python workflow/scripts/normalize_narrowpeak.py \
            --input {output.narrowpeak} --output {output.narrowpeak} 2>> {log}
        if [ ! -s {output.summits} ]; then
            echo "error: MACS3 called peaks but did not publish a nonempty summit BED." >> {log}
            exit 1
        fi
        """


def genrich_input_bams(wildcards):
    # Genrich requires name-sorted BAMs; we sort the filtered BAMs by name.
    return expand(
        f"{PROCESSED}/namesort/{{s}}.namesorted.bam",
        s=samples_in_condition(wildcards.cond),
    )


rule namesort_for_genrich:
    input:
        f"{PROCESSED}/filtered/{{sample}}.filtered.bam",
    output:
        temp(f"{PROCESSED}/namesort/{{sample}}.namesorted.bam"),
    threads: config["resources"]["sort_threads"]
    log:
        f"{LOGS}/peaks/namesort_{{sample}}.log",
    conda:
        "../envs/align.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output})
        samtools sort -n -@ {threads} -o {output} {input} 2> {log}
        """


rule genrich_callpeak:
    input:
        bams=genrich_input_bams,
        scripts=script_inputs("normalize_narrowpeak.py"),
    output:
        narrowpeak=f"{PROCESSED}/peaks/genrich/{{cond}}.narrowPeak",
    params:
        extra=config["peaks"]["genrich_extra"],
        joined=lambda wc, input: ",".join(input.bams),
    log:
        f"{LOGS}/peaks/genrich_{{cond}}.log",
    conda:
        "../envs/peaks.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.narrowpeak})
        Genrich -t {params.joined} -o {output.narrowpeak} {params.extra} 2> {log}
        n=$(awk 'NF >= 3 && $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ && $3 > $2 {{ n++ }} END {{ print n+0 }}' \
            {output.narrowpeak})
        total=$(wc -l < {output.narrowpeak})
        echo "Genrich: $n valid peaks ($total nonempty records) for {wildcards.cond}" >> {log}
        if [ "$n" -eq 0 ] || [ "$n" -ne "$total" ]; then
            echo "error: Genrich produced an empty or invalid narrowPeak for {wildcards.cond}." >> {log}
            exit 1
        fi
        python workflow/scripts/normalize_narrowpeak.py \
            --input {output.narrowpeak} --output {output.narrowpeak} 2>> {log}
        """
