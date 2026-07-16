# =============================================================================
# trim.smk — adapter & quality trimming with fastp
# =============================================================================
# fastp auto-detects Nextera/Tn5 adapters for paired-end data and emits a JSON
# report that MultiQC ingests. If trimming.enabled is false, reads are passed
# through unchanged (symlinked) so downstream rules have a stable path.

rule fastp:
    input:
        r1=lambda wc: raw_fastqs(wc.sample)[0],
        r2=lambda wc: raw_fastqs(wc.sample)[1],
    output:
        r1=f"{PROCESSED}/trimmed/{{sample}}_R1.trimmed.fastq.gz",
        r2=f"{PROCESSED}/trimmed/{{sample}}_R2.trimmed.fastq.gz",
        json=f"{PROCESSED}/qc/fastp/{{sample}}.fastp.json",
        html=f"{PROCESSED}/qc/fastp/{{sample}}.fastp.html",
    params:
        enabled=config["trimming"]["enabled"],
        extra=config["trimming"]["extra"],
    threads: config["resources"]["general_threads"]
    log:
        f"{LOGS}/trim/{{sample}}.log",
    conda:
        "../envs/qc.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.r1}) $(dirname {output.json})
        if [ "{params.enabled}" = "True" ]; then
            fastp -i {input.r1} -I {input.r2} \
                  -o {output.r1} -O {output.r2} \
                  --thread {threads} {params.extra} \
                  -j {output.json} -h {output.html} 2> {log}
        else
            # Downstream paths always end in .gz. Preserve real gzip inputs as
            # symlinks, but compress plain FASTQ rather than giving it a false
            # suffix that Bowtie2 would try to decompress.
            r1_part="{output.r1}.part.$$"
            r2_part="{output.r2}.part.$$"
            cleanup() {{ rm -f -- "$r1_part" "$r2_part"; }}
            trap cleanup EXIT
            pass_fastq() {{
                src="$1"; dest="$2"; part="$3"
                case "$src" in
                    *.gz)
                        gzip -t "$src"
                        resolved=$(python -c \
                            'import os,sys; sys.stdout.write(os.path.realpath(sys.argv[1]))' \
                            "$src")
                        ln -sfn "$resolved" "$dest"
                        ;;
                    *)
                        gzip -c "$src" > "$part"
                        gzip -t "$part"
                        mv -f "$part" "$dest"
                        ;;
                esac
            }}
            pass_fastq {input.r1:q} {output.r1:q} "$r1_part"
            pass_fastq {input.r2:q} {output.r2:q} "$r2_part"
            echo '{{"summary":{{"note":"trimming disabled"}}}}' > {output.json}
            echo "<html><body>trimming disabled</body></html>" > {output.html}
        fi
        """
