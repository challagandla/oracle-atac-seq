import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def source(relative_path):
    return (ROOT / relative_path).read_text(encoding="utf-8")


def rule_text(relative_path, rule_name):
    text = source(relative_path)
    match = re.search(
        rf"^rule {re.escape(rule_name)}:\n(?P<body>.*?)(?=^rule \w+:|\Z)",
        text,
        flags=re.MULTILINE | re.DOTALL,
    )
    assert match, f"rule {rule_name!r} not found in {relative_path}"
    return match.group("body")


def test_fastqc_jobs_are_isolated_and_publish_declared_outputs():
    fastqc = rule_text("workflow/rules/qc.smk", "fastqc_raw")

    assert "mktemp -d" in fastqc
    assert "{wildcards.sample}_{wildcards.read}.XXXXXX" in fastqc
    assert 'staged="$tmpdir/{wildcards.sample}_{wildcards.read}${{suffix}}"' in fastqc
    assert "os.path.realpath(sys.argv[1])" in fastqc
    assert 'ln -s "$resolved" "$staged"' in fastqc
    assert 'fastqc --memory 10000 -o "$tmpdir" "$staged"' in fastqc
    assert 'src_html="$tmpdir/{wildcards.sample}_{wildcards.read}_fastqc.html"' in fastqc
    assert 'src_zip="$tmpdir/{wildcards.sample}_{wildcards.read}_fastqc.zip"' in fastqc
    assert 'mv -f "$src_html" {output.html:q}' in fastqc
    assert 'mv -f "$src_zip" {output.zip:q}' in fastqc
    assert "fastqc --memory 10000 -o {params.outdir}" not in fastqc


def test_every_supported_raw_fastq_suffix_is_ignored_except_test_fixtures():
    lines = source(".gitignore").splitlines()
    patterns = ["*.fastq", "*.fq", "*.fastq.gz", "*.fq.gz"]

    for pattern in patterns:
        assert pattern in lines
    last_pattern = max(lines.index(pattern) for pattern in patterns)
    assert lines.index("!/.test/data/") > last_pattern
    assert lines.index("!/.test/data/**") > last_pattern


def test_shifted_bams_restore_and_verify_per_record_read_groups():
    shift = rule_text("workflow/rules/filter.smk", "tn5_shift")

    alignment_sieve = shift.index("alignmentSieve")
    restore = shift.index("samtools addreplacerg", alignment_sieve)
    sort = shift.index("samtools sort", restore)
    verify = shift.index("samtools view -c -d RG:{wildcards.sample}", sort)
    assert alignment_sieve < restore < sort < verify
    assert "-m overwrite_all -w" in shift
    assert "-O BAM -o {output.bam}.rg.tmp.bam" in shift
    assert "input_total=$(samtools view -c {input.bam})" in shift
    assert '"$total" -ne "$input_total"' in shift
    assert "shifted BAM postcondition failed" in shift


def test_peak_callers_validate_and_normalize_narrowpeak_outputs():
    for rule_name in ("macs3_callpeak", "genrich_callpeak"):
        rule = rule_text("workflow/rules/peaks.smk", rule_name)
        assert 'scripts=script_inputs("normalize_narrowpeak.py")' in rule
        assert "python workflow/scripts/normalize_narrowpeak.py" in rule
        assert "--input {output.narrowpeak} --output {output.narrowpeak}" in rule


def test_fragment_size_rule_declares_every_generated_artifact():
    overlay = rule_text("workflow/rules/report.smk", "fragment_size_overlay")

    assert 'png=f"{RESULTS}/figures/fragment_size_distribution.png"' in overlay
    assert 'pdf=f"{RESULTS}/figures/fragment_size_distribution.pdf"' in overlay
    assert 'mqc=f"{RESULTS}/figures/fragment_size_distribution_mqc.png"' in overlay
    assert 'table=f"{PROCESSED}/qc/fragsize/fragment_partitions_mqc.tsv"' in overlay
    assert "{output.png} {output.pdf} {output.mqc} {output.table}" in overlay


def test_fragment_size_pdf_is_a_final_report_target():
    report = source("workflow/rules/report.smk")

    collector = report.split("def report_outputs():", 1)[1]
    assert 'f"{RESULTS}/figures/fragment_size_distribution.pdf"' in collector


def test_aggregate_tss_matrix_honours_the_configured_region_bound():
    regions = rule_text("workflow/rules/report.smk", "tss_qc_regions")
    matrix = rule_text("workflow/rules/report.smk", "tss_matrix_all")
    per_sample = rule_text("workflow/rules/qc.smk", "tss_enrichment")

    assert 'max_regions=config.get("qc", {}).get("tss", {}).get("max_regions", 0)' in regions
    assert 'scripts=script_inputs("select_tss_regions.py")' in regions
    assert "--max-regions {params.max_regions}" in regions
    shared = 'f"{PROCESSED}/qc/tss/qc_regions.bed"'
    assert shared in matrix
    assert shared in per_sample
    assert "head -n" not in matrix
    assert "head -n" not in per_sample


def test_fragment_size_multiqc_image_is_in_exact_manifest():
    multiqc = source("workflow/rules/multiqc.smk")

    review_inputs = multiqc.split("def qc_review_inputs(wildcards):", 1)[1].split(
        "def multiqc_inputs(wildcards):", 1
    )[0]
    assert 'f"{RESULTS}/figures/fragment_size_distribution_mqc.png"' in review_inputs
    assert "multiqc --file-list {input.manifest}" in multiqc


def test_qc_review_report_is_outcome_blind_and_a_final_target():
    multiqc = source("workflow/rules/multiqc.smk")
    common = source("workflow/rules/common.smk")

    review_inputs = multiqc.split("def qc_review_inputs(wildcards):", 1)[1].split(
        "def multiqc_inputs(wildcards):", 1
    )[0]
    final_inputs = multiqc.split("def multiqc_inputs(wildcards):", 1)[1].split(
        "rule qc_review_manifest:", 1
    )[0]
    assert "diffacc" not in review_inputs
    assert "chromvar" not in review_inputs
    assert "consensus_counts" not in review_inputs
    assert "qc_review_inputs(wildcards)" in final_inputs
    assert "consensus_counts.tsv.summary" in final_inputs
    assert 'f"{RESULTS}/qc/qc_review_report.html"' in common


def test_enrichment_tables_are_declared_and_recreated_per_run():
    enrichment = rule_text("workflow/rules/enrich.smk", "functional_enrichment")

    assert 'tables=directory(f"{RESULTS}/enrichment/tables")' in enrichment
    assert "rm -rf {output.tables}" in enrichment
    assert "--out-dir {output.tables}" in enrichment


def test_enrichment_mapping_coverage_uses_the_exact_tested_peak_universe():
    enrichment = source("workflow/scripts/functional_enrichment.R")

    assert "tested_peaks <- unique(as.character(tested_bed[[4]]))" in enrichment
    assert "source_gene %in% valid_entrez_ids" in enrichment
    assert (
        "mapping_coverage <- length(mapped_tested_peaks) / length(tested_peaks)"
        in enrichment
    )
    assert "sum(entrez_like) / nrow(ann)" not in enrichment
    assert "if (mapping_coverage < 0.5)" in enrichment


def test_enrichment_requires_orgdb_and_validates_actual_entrez_keys():
    enrichment = source("workflow/scripts/functional_enrichment.R")

    assert "if (!have_org)" in enrichment
    assert "functional enrichment requires a matching OrgDb" in enrichment
    assert 'AnnotationDbi::keytypes(orgdb)' in enrichment
    assert 'AnnotationDbi::keys(orgdb, keytype = "ENTREZID")' in enrichment
    assert "source_gene %in% valid_entrez_ids" in enrichment


def test_enrichment_checks_mapping_coverage_within_each_target_direction():
    enrichment = source("workflow/scripts/functional_enrichment.R")

    assert "gene_of <- function(bed, direction)" in enrichment
    assert "p <- unique(as.character" in enrichment
    assert (
        "peak_mapping_coverage <- length(intersect(p, mapped_peaks)) / length(p)"
        in enrichment
    )
    assert "if (peak_mapping_coverage < 0.5)" in enrichment
    assert "if (peak_mapping_coverage < 0.8)" in enrichment
    assert 'gene_of(opt$up, "up")' in enrichment
    assert 'gene_of(opt$down, "down")' in enrichment


def test_homer_direction_directory_is_recreated_before_writing_results():
    motif = rule_text("workflow/rules/motif.smk", "homer_motifs")

    cleanup = motif.index("rm -rf {params.outdir} {params.preparsed}")
    create = motif.index("mkdir -p {params.outdir} {params.preparsed}")
    inspect_peaks = motif.index("n=$(wc -l < {input.peaks})")
    assert cleanup < create < inspect_peaks


def test_tested_peak_opportunity_universe_is_declared_and_shared_downstream():
    deseq2_rule = rule_text("workflow/rules/diffacc.smk", "deseq2")
    deseq2_script = source("workflow/scripts/differential_accessibility.R")
    motif = rule_text("workflow/rules/motif.smk", "homer_motifs")
    enrichment = rule_text("workflow/rules/enrich.smk", "functional_enrichment")

    assert 'tested=f"{RESULTS}/diffacc/tested_peaks.bed"' in deseq2_rule
    assert "tested <- res_df[is.finite(res_df$pvalue)," in deseq2_script
    assert 'P("tested_peaks.bed")' in deseq2_script
    assert 'tested=f"{RESULTS}/diffacc/tested_peaks.bed"' in motif
    assert "{input.peaks} {input.tested}" in motif
    assert "-bg {params.outdir}/background.bed" in motif
    assert 'tested=f"{RESULTS}/diffacc/tested_peaks.bed"' in enrichment
    assert "--tested {input.tested}" in enrichment


def test_lfc_shrinkage_is_required_and_never_silently_substituted():
    deseq2 = source("workflow/scripts/differential_accessibility.R")

    assert "res_shrunk <- lfcShrink(" in deseq2
    assert 'type = "ashr", res = res' in deseq2
    assert "using unshrunken LFC" not in deseq2
    assert "res_shrunk <- tryCatch(" not in deseq2


def test_deseq2_figures_use_dynamic_factor_labels_and_finite_pvalue_floors():
    deseq2 = source("workflow/scripts/differential_accessibility.R")

    assert "colnames(ann) <- opt$factor" in deseq2
    assert "ann_colours <- setNames(list(cond_col), opt$factor)" in deseq2
    assert "positive_padj <- vol$padj[is.finite(vol$padj) & vol$padj > 0]" in deseq2
    assert ".Machine$double.xmin" in deseq2
    assert "min(vol$padj[vol$padj > 0]" not in deseq2


def test_volcano_does_not_draw_raw_lfc_cutoffs_on_the_shrunken_lfc_axis():
    deseq2 = source("workflow/scripts/differential_accessibility.R")
    volcano = deseq2.split("# ---- volcano plot", 1)[1]

    assert "geom_vline" not in volcano
    assert 'x = "log2 fold-change (shrunken)"' in volcano
    assert 'Calls use padj < %.3g and |raw LFC| >= %.3g' in volcano


def test_raw_pvalue_histogram_is_independent_of_adjusted_pvalue_availability():
    deseq2 = source("workflow/scripts/differential_accessibility.R")

    volcano_branch = deseq2.split("vol <- res_df[!is.na(res_df$padj), ]", 1)[1]
    after_branch = volcano_branch.split(
        "# Use every finite raw p-value.", 1
    )[1]
    assert "pvals <- res_df[is.finite(res_df$pvalue)," in after_branch
    assert 'P("pvalue_histogram.pdf")' in after_branch
    assert "No finite raw p-values" in after_branch


def test_zero_raw_lfc_is_never_misclassified_as_down():
    deseq2 = source("workflow/scripts/differential_accessibility.R")

    assert 'res_df$direction <- "ns"' in deseq2
    assert 'res_df$direction[eligible & res_df$log2FoldChange > 0] <- "up"' in deseq2
    assert 'res_df$direction[eligible & res_df$log2FoldChange < 0] <- "down"' in deseq2
    assert 'ifelse(res_df$log2FoldChange > 0, "up", "down")' not in deseq2


def test_bindetect_directory_is_recreated_before_writing_results():
    bindetect = rule_text("workflow/rules/footprint.smk", "tobias_bindetect")

    cleanup = bindetect.index("rm -rf {params.outdir}")
    create = bindetect.index("mkdir -p {params.outdir}")
    execute = bindetect.index("TOBIAS BINDetect")
    assert cleanup < create < execute


def test_sra_download_uses_fresh_attempt_directory_and_atomic_outputs():
    download = rule_text("workflow/rules/download.smk", "sra_download")

    assert "mktemp -d" in download
    assert "trap cleanup EXIT" in download
    assert ".part.$$" in download
    assert "gzip -t" in download
    assert "_tmp_{wc.sample}" not in download


def test_trim_disabled_preserves_gzip_and_atomically_compresses_plain_fastq():
    fastp = rule_text("workflow/rules/trim.smk", "fastp")

    assert 'case "$src" in' in fastp
    assert "*.gz)" in fastp
    assert 'gzip -t "$src"' in fastp
    assert "os.path.realpath(sys.argv[1])" in fastp
    assert "readlink -f" not in fastp
    assert 'ln -sfn "$resolved" "$dest"' in fastp
    assert 'gzip -c "$src" > "$part"' in fastp
    assert 'gzip -t "$part"' in fastp
    assert 'mv -f "$part" "$dest"' in fastp
    assert "trap cleanup EXIT" in fastp


def test_disabled_qc_placeholders_use_the_declared_python_runtime():
    qc = source("workflow/rules/qc.smk")

    assert "base64 -d" not in qc
    assert qc.count("base64.b64decode(sys.stdin.read())") == 2


def test_footprint_motif_database_is_a_tracked_input():
    bindetect = rule_text("workflow/rules/footprint.smk", "tobias_bindetect")
    provenance = source("workflow/rules/provenance.smk")

    assert 'motifs=(config["footprinting"]["motif_db"]' in bindetect
    assert 'if config["footprinting"].get("enabled", False) else [])' in bindetect
    assert "--motifs {input.motifs}" in bindetect
    assert 'paths.append(motif_db)' in provenance
    assert 'reference_inputs.append(motif_db)' in provenance


def test_pinned_tobias_atacorrect_publishes_and_checks_all_documented_outputs():
    atacorrect = rule_text("workflow/rules/footprint.smk", "tobias_atacorrect")
    environment = source("workflow/envs/footprint.yaml")

    assert "- tobias=0.16.1" in environment
    expected_outputs = {
        "corrected": "{{cond}}_corrected.bw",
        "uncorrected": "{{cond}}_uncorrected.bw",
        "bias": "{{cond}}_bias.bw",
        "expected": "{{cond}}_expected.bw",
        "diagnostic": "{{cond}}_atacorrect.pdf",
    }
    for label, filename in expected_outputs.items():
        assert f'{label}=f"{{RESULTS}}/footprint/{filename}"' in atacorrect

    assert "--outdir {params.outdir}" in atacorrect
    assert "--prefix {params.prefix}" in atacorrect
    output_check = atacorrect.split("for artifact in", 1)[1]
    for label in expected_outputs:
        assert f"{{output.{label}}}" in output_check
    assert 'if [ ! -s "$artifact" ]' in output_check


def test_tobias_atacorrect_tracks_the_temporary_merged_bam_index():
    merge = rule_text("workflow/rules/footprint.smk", "tobias_merge_condition")
    atacorrect = rule_text("workflow/rules/footprint.smk", "tobias_atacorrect")

    index_path = 'f"{RESULTS}/footprint/{{cond}}.merged.bam.bai"'
    assert f"bai=temp({index_path})" in merge
    assert f"bai={index_path}" in atacorrect


def test_fingerprint_uses_the_configured_mapping_quality():
    fingerprint = rule_text("workflow/rules/report.smk", "plot_fingerprint")

    assert 'mapq=config["filtering"]["min_mapq"]' in fingerprint
    assert "--minMappingQuality {params.mapq}" in fingerprint


def test_outcome_blind_sample_similarity_is_in_qc_review():
    similarity = rule_text("workflow/rules/report.smk", "qc_sample_similarity")
    multiqc = source("workflow/rules/multiqc.smk")

    assert "multiBigwigSummary bins" in similarity
    assert "--corMethod spearman" in similarity
    assert "plotPCA" in similarity
    assert "diffacc" not in similarity
    assert "qc_sample_correlation_mqc.png" in multiqc
    assert "qc_sample_pca_mqc.png" in multiqc


def test_genrich_must_publish_a_nonempty_valid_peak_file():
    genrich = rule_text("workflow/rules/peaks.smk", "genrich_callpeak")

    assert "Genrich: $n valid peaks" in genrich
    assert 'if [ "$n" -eq 0 ] || [ "$n" -ne "$total" ]' in genrich


def test_macs3_summits_are_declared_and_must_be_nonempty():
    macs3 = rule_text("workflow/rules/peaks.smk", "macs3_callpeak")

    assert 'summits=f"{PROCESSED}/peaks/macs3/{{sample}}_summits.bed"' in macs3
    assert "macs3 callpeak" in macs3
    assert "-f BAMPE" in macs3
    assert "if [ ! -s {output.summits} ]" in macs3


def test_diffacc_declares_every_regular_png_created_by_the_figure_helpers():
    deseq2 = rule_text("workflow/rules/diffacc.smk", "deseq2")

    expected = [
        "PCA_plot", "MA_plot", "volcano_plot", "scree_plot",
        "sample_correlation_heatmap", "sample_distance_heatmap",
        "differential_peaks_heatmap", "pvalue_histogram",
    ]
    assert 'png=expand(' in deseq2
    for stem in expected:
        assert f'"{stem}"' in deseq2


def test_chromvar_declares_pngs_and_multiqc_tracks_its_report_images():
    chromvar = rule_text("workflow/rules/chromvar.smk", "chromvar")
    multiqc = source("workflow/rules/multiqc.smk")

    for name in (
        "chromvar_variability.png",
        "chromvar_variability_mqc.png",
        "chromvar_deviation_heatmap.png",
        "chromvar_deviation_heatmap_mqc.png",
    ):
        assert name in chromvar
    for name in (
        "chromvar_variability_mqc.png",
        "chromvar_deviation_heatmap_mqc.png",
    ):
        assert name in multiqc


def test_chromvar_exports_motif_ids_and_uses_a_seeded_serial_backend():
    chromvar = source("workflow/scripts/chromvar_analysis.R")
    environment = source("workflow/envs/r.yaml")

    assert "variability$motif <- rownames(variability)" in chromvar
    assert "write.table(variability, opt$out_var" in chromvar
    assert "BiocParallel::register(" in chromvar
    assert "BiocParallel::SerialParam(RNGseed = opt$seed" in chromvar
    assert "bioconductor-biocparallel" in environment


def test_provenance_hashes_each_selected_fastq_and_publishes_the_table():
    provenance = source("workflow/rules/provenance.smk")
    common = source("workflow/rules/common.smk")

    assert "return [path for sample in SAMPLES for path in raw_fastqs(sample)]" in provenance
    assert "raw=_provenance_raw_inputs" in provenance
    assert 'raw=f"{RESULTS}/provenance/raw_inputs.sha256.tsv"' in provenance
    assert '"sha256": _sha256_file(path)' in provenance
    assert '"raw_inputs": raw_records' in provenance
    assert 'f"{RESULTS}/provenance/raw_inputs.sha256.tsv"' in common


def test_provenance_reads_snakemake_version_from_the_running_package():
    provenance = source("workflow/rules/provenance.smk")

    assert "import importlib.metadata" in provenance
    assert '_installed_package_version("snakemake")' in provenance
    assert 'command_version("snakemake", "--version")' not in provenance


def test_reference_derived_outputs_create_the_reference_directory_explicitly():
    chrom_sizes = rule_text("workflow/rules/refs.smk", "chrom_sizes")
    tss = rule_text("workflow/rules/refs.smk", "tss_bed")

    assert "mkdir -p {REF}" in chrom_sizes
    assert "mkdir -p {REF}" in tss
    assert 'chrom=f"{REF}/chrom.sizes"' in tss
    assert "--chrom-sizes {input.chrom}" in tss
    assert "test -s {output}" in tss


def test_featurecounts_isolated_per_sample_then_validated_before_merge():
    per_sample = rule_text("workflow/rules/counts.smk", "featurecounts_sample")
    merge = rule_text("workflow/rules/counts.smk", "merge_featurecounts")

    assert "{{sample}}.counts.tsv" in per_sample
    assert "oracle-atac-featurecounts.{wildcards.sample}.XXXXXX" in per_sample
    assert 'script=script_inputs("run_featurecounts.R")' in per_sample
    assert '"../envs/counts.yaml"' in per_sample
    assert 'export R_LIBS_USER="$CONDA_PREFIX/lib/R/library"' in per_sample
    assert "Rscript --vanilla {input.script:q}" in per_sample
    assert "{input.bam:q}" in per_sample
    assert "{input.bams}" not in per_sample
    assert 'script=script_inputs("merge_featurecounts.py")' in merge
    assert "--counts {input.counts:q}" in merge
    assert "--summaries {input.summaries:q}" in merge
    assert "--samples {params.samples}" in merge
    count_environment = source("workflow/envs/counts.yaml")
    assert "r-base=4.4" in count_environment
    assert "bioconductor-rsubread=2.20.0" in count_environment
    assert "subread=" not in source("workflow/envs/peaks.yaml")


def test_rsubread_counter_preserves_paired_fragment_semantics_and_validates_outputs():
    counter = source("workflow/scripts/run_featurecounts.R")

    for contract in (
        "isPairedEnd = TRUE",
        "countReadPairs = TRUE",
        "requireBothEndsMapped = TRUE",
        "countChimericFragments = FALSE",
        "allowMultiOverlap = FALSE",
        "countMultiMappingReads = FALSE",
        'match("Assigned", fc$stat$Status)',
        "sum(fc$counts[, 1]) != stat_values[[assigned_index]]",
        "!identical(fc$annotation, expected_annotation)",
    ):
        assert contract in counter


def test_count_matrix_consumers_use_authoritative_merged_sample_ids():
    for script in (
        "workflow/scripts/differential_accessibility.R",
        "workflow/scripts/chromvar_analysis.R",
    ):
        text = source(script)
        assert 'count_samples <- colnames(' in text
        assert "setequal(count_samples, expected_samples)" in text
        assert "match(count_samples," in text
        assert "non-negative integer counts" in text
        assert 'sub("\\\\.filtered\\\\.bam$"' not in text


def test_bowtie2_large_index_mode_tracks_the_correct_suffix():
    common = source("workflow/rules/common.smk")
    refs = rule_text("workflow/rules/refs.smk", "bowtie2_build")
    align = rule_text("workflow/rules/align.smk", "bowtie2_align")

    assert 'suffix = "bt2l"' in common
    assert 'return "--large-index"' in common
    assert "bowtie2_index_files()" in refs
    assert "{params.mode}" in refs
    assert "idx=bowtie2_index_files()" in align


def test_tss_qc_description_matches_ten_base_bin_size():
    qc = source("workflow/rules/qc.smk")

    assert "10-bp-binned Tn5 insertion signal" in qc
    assert "base-resolution Tn5 insertion" not in qc
    assert "--binSize 10" in rule_text("workflow/rules/qc.smk", "tss_enrichment")


def test_annotation_does_not_drop_valid_custom_contigs_as_nonstandard():
    annotation = source("workflow/scripts/annotate_peaks.R")

    assert "keepStandardChromosomes" not in annotation
    assert "shared_peak_count" in annotation
    assert "for (style in unique(target_styles))" in annotation
    assert "keepSeqlevels(best, shared" in annotation


def test_annotation_fails_on_severe_txdb_contig_loss():
    annotation = source("workflow/scripts/annotate_peaks.R")

    assert "retained_fraction <- length(peaks) / n_before" in annotation
    assert "if (retained_fraction < 0.5)" in annotation
    assert "if (retained_fraction < 0.8)" in annotation
    assert "refusing a severely biased subset" in annotation


def test_chromvar_does_not_drop_valid_custom_contigs_as_nonstandard():
    chromvar = source("workflow/scripts/chromvar_analysis.R")

    assert "keepStandardChromosomes" not in chromvar
    assert "shared_peak_count" in chromvar
    assert "for (style in unique(target_styles))" in chromvar
    assert "keepSeqlevels(best, shared" in chromvar


def test_chromvar_aligns_peak_ranges_to_count_ids_before_constructing_assay():
    chromvar = source("workflow/scripts/chromvar_analysis.R")

    uniqueness = chromvar.index("if (anyDuplicated(bed_ids) || anyDuplicated(count_ids))")
    same_set = chromvar.index("if (!setequal(bed_ids, count_ids))")
    reorder = chromvar.index("bed <- bed[match(count_ids, bed_ids), , drop = FALSE]")
    construct = chromvar.index("se <- SummarizedExperiment")
    assert uniqueness < same_set < reorder < construct


def test_chromvar_fails_on_severe_bsgenome_contig_loss():
    chromvar = source("workflow/scripts/chromvar_analysis.R")

    assert "retained_fraction <- nrow(se) / n_before" in chromvar
    assert "if (retained_fraction < 0.5)" in chromvar
    assert "if (retained_fraction < 0.8)" in chromvar
    assert "refusing a severely biased subset" in chromvar
