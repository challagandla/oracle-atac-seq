import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def rule_text(relative_path, rule_name):
    text = (ROOT / relative_path).read_text(encoding="utf-8")
    match = re.search(
        rf"^rule {re.escape(rule_name)}:\n(?P<body>.*?)(?=^rule \w+:|\Z)",
        text,
        flags=re.MULTILINE | re.DOTALL,
    )
    assert match, f"rule {rule_name!r} not found in {relative_path}"
    return match.group("body")


def test_new_alignments_receive_library_read_group_metadata():
    align = rule_text("workflow/rules/align.smk", "bowtie2_align")

    assert "--rg-id {wildcards.sample}" in align
    assert "--rg SM:{wildcards.sample}" in align
    assert "--rg LB:{wildcards.sample}" in align
    assert "--rg PL:ILLUMINA" in align


def test_alignment_sort_uses_an_isolated_restart_safe_temp_directory():
    align = rule_text("workflow/rules/align.smk", "bowtie2_align")

    assert 'tmp_root=f"{PROCESSED}/aligned/.tmp"' in align
    assert 'mktemp -d "{params.tmp_root}/{wildcards.sample}.XXXXXX"' in align
    assert 'trap cleanup EXIT' in align
    assert 'samtools sort -@ {threads} -T "$tmpdir/chunk"' in align
    assert "samtools quickcheck {output.bam}" in align


def test_legacy_bams_are_normalized_before_picard():
    normalize = rule_text("workflow/rules/filter.smk", "normalize_read_groups")
    mark_duplicates = rule_text("workflow/rules/filter.smk", "mark_duplicates")

    assert "samtools addreplacerg" in normalize
    assert "-m overwrite_all" in normalize
    # Fresh Bowtie2 BAMs already contain this ID; -w replaces the existing
    # header line while overwrite_all normalizes every record's RG tag.
    assert "-m overwrite_all -w" in normalize
    assert "SM:{wildcards.sample}" in normalize
    assert "LB:{wildcards.sample}" in normalize
    assert "PL:ILLUMINA" in normalize
    assert "samtools quickcheck" in normalize
    assert ".nomito.rg.bam" in mark_duplicates


def test_alignment_filter_flags_follow_duplicate_removal_setting():
    filter_bam = rule_text("workflow/rules/filter.smk", "filter_bam")

    assert 'proper="-f 2"' in filter_bam
    assert 'keep_proper_pairs"] else' not in filter_bam
    # 2828 = unmapped + mate-unmapped + secondary + QC-fail + supplementary.
    # Duplicate (1024) is added only when the user requested duplicate removal.
    assert '3852 if config["filtering"]["remove_duplicates"] else 2828' in filter_bam
    assert "-F {params.exclude_flags}" in filter_bam


def test_library_complexity_retains_the_index_used_by_pysam_fetch():
    complexity = rule_text("workflow/rules/report.smk", "library_complexity")

    assert 'bam=f"{PROCESSED}/filtered/{{sample}}.nomito.bam"' in complexity
    assert 'bai=f"{PROCESSED}/filtered/{{sample}}.nomito.bam.bai"' in complexity


def test_mitochondrial_filtering_and_qc_share_the_resolved_contig_list():
    remove_mito = rule_text("workflow/rules/filter.smk", "remove_mito")
    report = rule_text("workflow/rules/report.smk", "qc_summary_sample")
    peaks = (ROOT / "config/config.yaml").read_text(encoding="utf-8")

    assert "MITOCHONDRIAL_CONTIGS" in remove_mito
    assert "bedtools pairtobed" in remove_mito
    assert "-type neither" in remove_mito
    assert "mitochondrial_contig_args" in report
    assert "--mitochondrial-contig" in report
    assert "-e chrM" not in peaks


def test_mitochondrial_filtering_fails_when_no_configured_contig_is_present():
    remove_mito = rule_text("workflow/rules/filter.smk", "remove_mito")

    message = remove_mito.index(
        "Mitochondrial filtering is enabled, but none of the configured contigs"
    )
    failure = remove_mito.index("exit 1", message)
    disabled = remove_mito.index("Mitochondrial filtering disabled", message)
    assert message < failure < disabled


def test_fragment_filters_check_both_mates_in_one_flagstat_scan():
    for rule_name in ("remove_mito", "remove_blacklist"):
        rule = rule_text("workflow/rules/filter.smk", rule_name)
        assert "pair_counts=$(samtools flagstat {output.bam})" in rule
        assert "'$NF == \"read1\" {{print $1}}'" in rule
        assert "'$NF == \"read2\" {{print $1}}'" in rule
        assert "after=$read1" in rule
        assert "samtools view -c -f 64 {output.bam}" not in rule


def test_enabled_blacklist_filtering_fails_closed_on_an_empty_bed():
    remove_blacklist = rule_text("workflow/rules/filter.smk", "remove_blacklist")

    enabled = remove_blacklist.index('if [ "{params.do_remove}" = "True" ]')
    empty = remove_blacklist.index('if [ ! -s "{params.bl}" ]', enabled)
    failure = remove_blacklist.index("exit 1", empty)
    disabled = remove_blacklist.index("cp {input.bam} {output.bam}", failure)
    assert enabled < empty < failure < disabled
