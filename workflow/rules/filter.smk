# =============================================================================
# filter.smk — post-alignment filtering + Tn5 shift
# =============================================================================
# ENCODE-informed bulk ATAC filtering:
#   1. keep properly-paired, primary alignments with MAPQ >= threshold
#   2. remove complete fragments on configured mitochondrial contigs
#   3. mark & remove PCR/optical duplicates (Picard)
#   4. remove complete fragments when either mate overlaps the blacklist
#   5. Tn5 shift: +4 (+ strand) / -5 (- strand) for base-pair-accurate cut sites
# Each step is logged for the MultiQC report.

rule filter_bam:
    input:
        bam=f"{PROCESSED}/aligned/{{sample}}.sorted.bam",
    output:
        bam=temp(f"{PROCESSED}/filtered/{{sample}}.namefilt.bam"),
    params:
        mapq=config["filtering"]["min_mapq"],
        proper="-f 2",
        exclude_flags=(3852 if config["filtering"]["remove_duplicates"] else 2828),
    threads: config["resources"]["sort_threads"]
    log:
        f"{LOGS}/filter/{{sample}}.mapq.log",
    conda:
        "../envs/align.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.bam})
        # Always exclude unmapped (4), mate-unmapped (8), secondary (256),
        # QC-fail (512), and supplementary (2048) alignments: -F 2828.
        # Add duplicate (1024) only when duplicate removal is requested: -F 3852.
        samtools view -@ {threads} -b {params.proper} -F {params.exclude_flags} \
            -q {params.mapq} {input.bam} > {output.bam} 2> {log}
        """


rule remove_mito:
    input:
        f"{PROCESSED}/filtered/{{sample}}.namefilt.bam",
    output:
        bam=temp(f"{PROCESSED}/filtered/{{sample}}.nomito.bam"),
        idx=temp(f"{PROCESSED}/filtered/{{sample}}.nomito.bam.bai"),
    params:
        do_remove=config["filtering"]["remove_mito"],
        mito_contigs=" ".join(shlex.quote(c) for c in MITOCHONDRIAL_CONTIGS),
    threads: config["resources"]["sort_threads"]
    log:
        f"{LOGS}/filter/{{sample}}.nomito.log",
    conda:
        "../envs/peaks.yaml"
    shell:
        r"""
        if [ "{params.do_remove}" = "True" ]; then
            tmp_prefix="{output.bam}.mito"
            mito_bed="$tmp_prefix.bed"
            rm -f "$mito_bed" "$tmp_prefix.namesort.bam" "$tmp_prefix.clean.bam"
            : > "$mito_bed"

            # Derive full-contig intervals from this BAM header. Exact string
            # matching supports custom accessions and avoids regex aliases.
            while IFS=$'\t' read -r contig length; do
                case " {params.mito_contigs} " in
                    *" $contig "*)
                        if [ "$length" -gt 0 ]; then
                            printf '%s\t0\t%s\n' "$contig" "$length" >> "$mito_bed"
                        fi
                        ;;
                esac
            done < <(
                samtools view -H {input} \
                  | awk -F '\t' '$1 == "@SQ" {{
                        sn=""; ln="";
                        for (i=2; i<=NF; i++) {{
                            if ($i ~ /^SN:/) sn=substr($i, 4);
                            if ($i ~ /^LN:/) ln=substr($i, 4);
                        }}
                        if (sn != "" && ln != "") print sn "\t" ln;
                    }}'
            )

            if [ -s "$mito_bed" ]; then
                before=$(samtools view -c -f 64 {input})
                samtools sort -n -@ {threads} -o "$tmp_prefix.namesort.bam" \
                    {input} 2> {log}
                bedtools pairtobed -abam "$tmp_prefix.namesort.bam" \
                    -b "$mito_bed" -type neither -ubam 2>> {log} \
                  | samtools sort -@ {threads} -o "$tmp_prefix.clean.bam" - 2>> {log}
                mv "$tmp_prefix.clean.bam" {output.bam}
                overlap=$(samtools view -c -L "$mito_bed" {output.bam})
                pair_counts=$(samtools flagstat {output.bam})
                read1=$(awk '$NF == "read1" {{print $1}}' <<< "$pair_counts")
                read2=$(awk '$NF == "read2" {{print $1}}' <<< "$pair_counts")
                if [ -z "$read1" ] || [ -z "$read2" ]; then
                    echo "mitochondrial filtering could not parse mate counts from samtools flagstat" >> {log}
                    exit 1
                fi
                after=$read1
                printf 'mitochondrial fragments removed: %s; retained: %s\n' \
                    "$((before - after))" "$after" >> {log}
                if [ "$overlap" -ne 0 ] || [ "$read1" -ne "$read2" ]; then
                    echo "mitochondrial filtering postcondition failed: overlap=$overlap read1=$read1 read2=$read2" >> {log}
                    exit 1
                fi
                rm -f "$tmp_prefix.namesort.bam"
            else
                printf 'Mitochondrial filtering is enabled, but none of the configured contigs were found in the BAM header: %s\n' \
                    "{params.mito_contigs}" > {log}
                printf 'Add the exact assembly mitochondrial accession to filtering.mitochondrial_contigs, or set remove_mito=false only for a reference that truly lacks one.\n' >> {log}
                rm -f "$mito_bed"
                exit 1
            fi
            rm -f "$mito_bed"
        else
            cp {input} {output.bam}
            printf 'Mitochondrial filtering disabled.\n' > {log}
        fi
        samtools index {output.bam}
        """


rule normalize_read_groups:
    input:
        bam=f"{PROCESSED}/filtered/{{sample}}.nomito.bam",
    output:
        bam=temp(f"{PROCESSED}/filtered/{{sample}}.nomito.rg.bam"),
    threads: config["resources"]["sort_threads"]
    log:
        f"{LOGS}/filter/{{sample}}.readgroups.log",
    conda:
        "../envs/align.yaml"
    shell:
        r"""
        # The sample sheet represents one biological library per BAM. Reapply
        # that contract here so retained BAMs made by an older workflow cannot
        # reach Picard without an ID/SM/LB/PL read group. overwrite_all also
        # repairs records whose RG tag points to a missing header entry.
        samtools addreplacerg -@ {threads} -m overwrite_all -w \
            -r $'@RG\tID:{wildcards.sample}\tSM:{wildcards.sample}\tLB:{wildcards.sample}\tPL:ILLUMINA' \
            -o {output.bam} {input.bam} 2> {log}
        samtools quickcheck {output.bam} 2>> {log}
        """


rule mark_duplicates:
    input:
        f"{PROCESSED}/filtered/{{sample}}.nomito.rg.bam",
    output:
        bam=temp(f"{PROCESSED}/filtered/{{sample}}.dedup.bam"),
        metrics=f"{PROCESSED}/qc/picard/{{sample}}.dup_metrics.txt",
    params:
        do_remove="true" if config["filtering"]["remove_duplicates"] else "false",
    log:
        f"{LOGS}/filter/{{sample}}.markdup.log",
    conda:
        "../envs/align.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.metrics})
        picard MarkDuplicates I={input} O={output.bam} \
            M={output.metrics} REMOVE_DUPLICATES={params.do_remove} \
            VALIDATION_STRINGENCY=LENIENT 2> {log}
        """


rule remove_blacklist:
    input:
        bam=f"{PROCESSED}/filtered/{{sample}}.dedup.bam",
        blacklist=(blacklist_bed()
                   if config["filtering"]["remove_blacklist"] and blacklist_bed()
                   else []),
    output:
        bam=f"{PROCESSED}/filtered/{{sample}}.filtered.bam",
        bai=f"{PROCESSED}/filtered/{{sample}}.filtered.bam.bai",
    params:
        do_remove=config["filtering"]["remove_blacklist"] and bool(blacklist_bed()),
        bl=blacklist_bed(),
    threads: config["resources"]["sort_threads"]
    log:
        f"{LOGS}/filter/{{sample}}.blacklist.log",
    conda:
        "../envs/peaks.yaml"
    shell:
        r"""
        if [ "{params.do_remove}" = "True" ]; then
            if [ ! -s "{params.bl}" ]; then
                echo "error: blacklist removal is enabled but the harmonized BED is empty." > {log}
                echo "       Provide a valid assembly-matched blacklist or disable removal explicitly." >> {log}
                exit 1
            fi
            before=$(samtools view -c -f 64 {input.bam})
            tmp_prefix="{output.bam}.blacklist"
            rm -f "$tmp_prefix.namesort.bam" "$tmp_prefix.clean.bam"

            # Filtering alignments independently leaves apparently proper pairs
            # with one mate missing. pairToBed examines name-grouped fragments
            # and retains a pair only when neither end overlaps the blacklist.
            samtools sort -n -@ {threads} -o "$tmp_prefix.namesort.bam" \
                {input.bam} 2> {log}
            bedtools pairtobed -abam "$tmp_prefix.namesort.bam" \
                -b {params.bl} -type neither -ubam 2>> {log} \
              | samtools sort -@ {threads} -o "$tmp_prefix.clean.bam" - 2>> {log}
            mv "$tmp_prefix.clean.bam" {output.bam}
            rm -f "$tmp_prefix.namesort.bam"
            samtools index {output.bam}

            # Verify the filter did what it says. bedtools matches intervals by
            # chromosome NAME: a UCSC blacklist against an Ensembl BAM removes
            # nothing and exits 0. Assert the postcondition instead of trusting it.
            left=$(samtools view -c -L {params.bl} {output.bam})
            pair_counts=$(samtools flagstat {output.bam})
            read1=$(awk '$NF == "read1" {{print $1}}' <<< "$pair_counts")
            read2=$(awk '$NF == "read2" {{print $1}}' <<< "$pair_counts")
            if [ -z "$read1" ] || [ -z "$read2" ]; then
                echo "error: could not parse mate counts from samtools flagstat." >> {log}
                exit 1
            fi
            after=$read1
            echo "blacklist: $before -> $after fragments ($((before - after)) removed)" >> {log}
            if [ "$left" -ne 0 ]; then
                echo "error: $left reads still overlap the blacklist after filtering." >> {log}
                echo "       The blacklist and the BAM disagree about chromosome names." >> {log}
                exit 1
            fi
            if [ "$read1" -ne "$read2" ]; then
                echo "error: blacklist filtering produced unequal mate counts ($read1 vs $read2)." >> {log}
                exit 1
            fi
        else
            cp {input.bam} {output.bam}
            samtools index {output.bam}
        fi
        """


rule flagstat_filtered:
    input:
        f"{PROCESSED}/filtered/{{sample}}.filtered.bam",
    output:
        f"{PROCESSED}/qc/flagstat/{{sample}}.filtered.flagstat",
    log:
        f"{LOGS}/qc/flagstat_filtered_{{sample}}.log",
    conda:
        "../envs/align.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output})
        samtools flagstat {input} > {output} 2> {log}
        """


# Tn5-shifted, coordinate-sorted BAM for cut-site QC and tracks.
# Peak callers deliberately use complete fragments from the unshifted BAM.
rule tn5_shift:
    input:
        bam=f"{PROCESSED}/filtered/{{sample}}.filtered.bam",
    output:
        bam=f"{PROCESSED}/shifted/{{sample}}.shifted.bam",
        bai=f"{PROCESSED}/shifted/{{sample}}.shifted.bam.bai",
    params:
        do_shift=config["filtering"]["tn5_shift"],
        chunk=config["filtering"].get("tn5_shift_genome_chunk_length", 50000000),
        tmpdir=lambda wc: f"{PROCESSED}/shifted/_tmp_{wc.sample}",
    threads: config["resources"]["sort_threads"]
    log:
        f"{LOGS}/filter/{{sample}}.tn5shift.log",
    conda:
        "../envs/deeptools.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.bam}) {params.tmpdir}
        if [ "{params.do_shift}" = "True" ]; then
            rm -rf {params.tmpdir}/deeptools_tmp
            mkdir -p {params.tmpdir}/deeptools_tmp
            export TMPDIR="{params.tmpdir}/deeptools_tmp"
            export TMP="$TMPDIR"
            export TEMP="$TMPDIR"
            ulimit -n 65536 || true
            # deepTools alignmentSieve applies the canonical +4/-5 Tn5 shift.
            alignmentSieve -b {input.bam} -o {output.bam}.tmp \
                --ATACshift --genomeChunkLength {params.chunk} \
                -p {threads} 2> {log}
            # alignmentSieve preserves the @RG header but drops RG tags from
            # records. Restore the one-library-per-BAM contract before these
            # files reach cut-site coverage and insertion-level QC.
            samtools addreplacerg -@ {threads} -m overwrite_all -w \
                -r $'@RG\tID:{wildcards.sample}\tSM:{wildcards.sample}\tLB:{wildcards.sample}\tPL:ILLUMINA' \
                -O BAM -o {output.bam}.rg.tmp.bam {output.bam}.tmp 2>> {log}
            samtools sort -@ {threads} -o {output.bam} {output.bam}.rg.tmp.bam 2>> {log}
            rm -f {output.bam}.tmp {output.bam}.rg.tmp.bam
            rm -rf {params.tmpdir}
        else
            cp {input.bam} {output.bam}
        fi
        samtools index {output.bam}
        samtools quickcheck {output.bam} 2>> {log}
        input_total=$(samtools view -c {input.bam})
        total=$(samtools view -c {output.bam})
        tagged=$(samtools view -c -d RG:{wildcards.sample} {output.bam})
        if [ "$total" -eq 0 ] || [ "$total" -ne "$input_total" ] || [ "$tagged" -ne "$total" ]; then
            echo "error: shifted BAM postcondition failed: input=$input_total output=$total tagged=$tagged" >> {log}
            exit 1
        fi
        """
