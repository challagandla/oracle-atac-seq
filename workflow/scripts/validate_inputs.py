#!/usr/bin/env python3
"""Validate the workflow configuration and biological sample sheet.

This module is imported while Snakemake builds the DAG.  Errors therefore stop
the run before references are downloaded or compute-heavy jobs are submitted.
"""
from __future__ import annotations

import math
import re
from pathlib import Path


SAFE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
SAFE_PATH = re.compile(r"^[A-Za-z0-9_./:+,@%=-]+$")
SAFE_CONTIG = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._:+-]*$")
SRA_RUN = re.compile(r"^[SED]RR[0-9]+$", re.IGNORECASE)
FASTQ_SUFFIX = re.compile(r".+\.(?:fastq|fq)(?:\.gz)?$")
TRUE_VALUES = {"1", "true", "yes", "y"}
FALSE_VALUES = {"0", "false", "no", "n"}
DEFAULT_MITOCHONDRIAL_CONTIGS = ["chrM", "MT", "chrMT", "Mito", "M"]


def _normalized_path(path):
    """Collapse relative components and symlink aliases for duplicate checks."""
    return str(Path(path).expanduser().resolve(strict=False))


def _is_gzip_reference(path, inspect_file):
    """Detect gzip references by common suffix or magic bytes when available."""
    text = str(path).strip()
    if text.lower().endswith((".gz", ".bgz", ".bgzf")):
        return True
    if inspect_file and text:
        try:
            with Path(text).open("rb") as handle:
                return handle.read(2) == b"\x1f\x8b"
        except OSError:
            return False
    return False


def _check_numeric(
    errors, label, value, *, integer=False, lower=None, upper=None,
    lower_open=False, upper_open=False,
):
    """Append one clear error when a configurable numeric value is unsafe."""
    try:
        if integer:
            # Keep the preflight and every CLI consumer on one contract. YAML
            # 100.0 is mathematically integral but argparse/head/tool thread
            # flags expect the lexical integer 100, so reject it early.
            if isinstance(value, bool) or not isinstance(value, int):
                raise ValueError
            parsed = value
        else:
            if isinstance(value, bool):
                raise ValueError
            parsed = float(value)
            if not math.isfinite(parsed):
                raise ValueError
    except (TypeError, ValueError, OverflowError):
        kind = "an integer" if integer else "a finite number"
        errors.append(f"{label} must be {kind} within its documented range")
        return None

    below = lower is not None and (parsed <= lower if lower_open else parsed < lower)
    above = upper is not None and (parsed >= upper if upper_open else parsed > upper)
    if below or above:
        left = "<" if lower_open else "<="
        right = "<" if upper_open else "<="
        if lower is not None and upper is not None:
            expected = f"{lower} {left} {label} {right} {upper}"
        elif lower is not None:
            expected = f"{label} {'>' if lower_open else '>='} {lower}"
        else:
            expected = f"{label} {'<' if upper_open else '<='} {upper}"
        errors.append(f"{expected}; got {value!r}")
        return None
    return parsed


def _enabled(section, default=False):
    if isinstance(section, dict):
        return bool(section.get("enabled", default))
    return bool(section)


def parse_include(value):
    """Return a strict boolean for the optional sample-sheet include column."""
    text = str(value).strip().lower()
    if text == "":
        return True
    if text in TRUE_VALUES:
        return True
    if text in FALSE_VALUES:
        return False
    raise ValueError(
        f"include must be one of {sorted(TRUE_VALUES | FALSE_VALUES)}; got {value!r}"
    )


def validate_config_and_samples(config, table, check_files=True):
    """Validate inputs and return the rows explicitly included for analysis."""
    errors = []
    required = {"sample", "condition", "replicate", "fq1", "fq2", "sra"}
    missing = sorted(required - set(table.columns))
    if missing:
        errors.append("sample sheet is missing columns: " + ", ".join(missing))
        raise ValueError("Invalid ATAC-seq inputs:\n- " + "\n- ".join(errors))

    if table.empty:
        errors.append("sample sheet has no data rows")

    for col in required:
        table[col] = table[col].fillna("").astype(str).str.strip()

    if table["sample"].duplicated().any():
        dup = sorted(table.loc[table["sample"].duplicated(False), "sample"].unique())
        errors.append("sample identifiers are not unique: " + ", ".join(dup))

    include = []
    for row_no, (_, row) in enumerate(table.iterrows(), start=2):
        sample = row["sample"]
        condition = row["condition"]
        try:
            selected_for_run = parse_include(row.get("include", ""))
        except ValueError as exc:
            errors.append(f"row {row_no} ({sample}): {exc}")
            selected_for_run = False
        include.append(selected_for_run)
        if not SAFE_ID.fullmatch(sample):
            errors.append(
                f"row {row_no}: sample {sample!r} must use only letters, numbers, ., _, or -"
            )
        if not SAFE_ID.fullmatch(condition):
            errors.append(
                f"row {row_no}: condition {condition!r} must use only letters, numbers, ., _, or -"
            )
        if not row["replicate"]:
            errors.append(f"row {row_no} ({sample}): replicate is empty")

        fq1, fq2, sra = row["fq1"], row["fq2"], row["sra"]
        local = bool(fq1 or fq2)
        remote = bool(sra)
        if local and remote:
            errors.append(f"row {row_no} ({sample}): provide FASTQs or an SRA run, not both")
        elif local and not (fq1 and fq2):
            errors.append(f"row {row_no} ({sample}): paired-end input needs both fq1 and fq2")
        elif not local and not remote:
            errors.append(f"row {row_no} ({sample}): provide fq1/fq2 or an SRA run accession")
        elif remote and not SRA_RUN.fullmatch(sra):
            errors.append(f"row {row_no} ({sample}): {sra!r} is not an SRR/ERR/DRR run accession")
        if fq1 and fq1 == fq2:
            errors.append(f"row {row_no} ({sample}): fq1 and fq2 point to the same file")
        for label, path in (("fq1", fq1), ("fq2", fq2)):
            if path and not SAFE_PATH.fullmatch(path):
                errors.append(
                    f"row {row_no} ({sample}): {label} contains whitespace or shell-special "
                    "characters; rename or link the file to a simple path"
                )
            elif path and not FASTQ_SUFFIX.fullmatch(path):
                errors.append(
                    f"row {row_no} ({sample}): {label} must end in .fastq, .fq, "
                    ".fastq.gz, or .fq.gz"
                )
        if check_files and selected_for_run and fq1:
            for label, path in (("fq1", fq1), ("fq2", fq2)):
                if not Path(path).is_file():
                    errors.append(f"row {row_no} ({sample}): {label} does not exist: {path}")

    table = table.copy()
    table["include"] = include
    selected = table.loc[table["include"]].copy()
    if selected.empty:
        errors.append("all libraries are excluded by the include column")
    # SRA Toolkit names its output from the canonical uppercase accession.
    # Normalize here so a valid lowercase user entry cannot create a path that
    # differs from fasterq-dump's output directory or FASTQ basename.
    selected.loc[:, "sra"] = selected["sra"].str.upper()
    if selected.duplicated(["condition", "replicate"]).any():
        pairs = selected.loc[
            selected.duplicated(["condition", "replicate"], False),
            ["condition", "replicate"],
        ].drop_duplicates()
        errors.append(
            "replicate identifiers repeat within a condition: "
            + ", ".join(f"{r.condition}/{r.replicate}" for r in pairs.itertuples())
        )

    # A row represents one independent biological library. Reusing the same
    # accession or paired FASTQ input under another sample name would turn a
    # duplicated file into an apparent replicate and invalidate consensus and
    # differential analyses. Only selected rows participate in this check so a
    # deliberately excluded audit row can remain documented in the sheet.
    sra_owners = {}
    fastq_pair_owners = {}
    fastq_file_owners = {}
    for row in selected.itertuples(index=False):
        if row.sra:
            sra_owners.setdefault(row.sra.upper(), []).append(row.sample)
        elif row.fq1 and row.fq2:
            pair = tuple(sorted((_normalized_path(row.fq1), _normalized_path(row.fq2))))
            fastq_pair_owners.setdefault(pair, []).append(row.sample)
            for path in set(pair):
                fastq_file_owners.setdefault(path, []).append(row.sample)
    for accession, owners in sorted(sra_owners.items()):
        if len(owners) > 1:
            errors.append(
                f"SRA run {accession} is reused by selected samples: "
                + ", ".join(sorted(owners))
            )
    for pair, owners in sorted(fastq_pair_owners.items()):
        if len(owners) > 1:
            errors.append(
                "the same paired FASTQ input is reused by selected samples "
                f"{', '.join(sorted(owners))}: {pair[0]}, {pair[1]}"
            )
    for path, owners in sorted(fastq_file_owners.items()):
        if len(owners) > 1:
            errors.append(
                f"FASTQ input file is reused by selected samples "
                f"{', '.join(sorted(owners))}: {path}"
            )

    diff = config.get("diffacc", {})
    if _enabled(diff):
        if str(diff.get("method", "DESeq2")) != "DESeq2":
            errors.append("diffacc.method must be 'DESeq2'")
        contrast_factor = None
        contrast = diff.get("contrast", [])
        if not isinstance(contrast, list) or len(contrast) != 3:
            errors.append("diffacc.contrast must be [factor, numerator, denominator]")
        else:
            factor, numerator, denominator = map(str, contrast)
            contrast_factor = factor
            if factor not in selected.columns:
                errors.append(f"diffacc contrast factor {factor!r} is not a sample-sheet column")
            else:
                levels = set(selected[factor].astype(str))
                absent = sorted({numerator, denominator} - levels)
                if absent:
                    errors.append("diffacc contrast levels absent from selected rows: " + ", ".join(absent))
                else:
                    level_counts = selected[factor].astype(str).value_counts()
                    insufficient = [
                        (level, int(level_counts.get(level, 0)))
                        for level in (numerator, denominator)
                        if int(level_counts.get(level, 0)) < 2
                    ]
                    if insufficient:
                        errors.append(
                            "differential contrast requires at least two independent "
                            f"libraries in each {factor!r} level; "
                            + ", ".join(f"{level}={count}" for level, count in insufficient)
                        )
                if numerator == denominator:
                    errors.append("diffacc numerator and denominator must differ")
        design = str(diff.get("design", "")).strip()
        if not design.startswith("~"):
            errors.append("diffacc.design must be an R formula beginning with '~'")
        else:
            # Accept ordinary R formula operators and functions, but make every
            # bare variable resolve to sample metadata before DESeq2 starts.
            missing_design = []
            for match in re.finditer(r"[A-Za-z][A-Za-z0-9_.]*", design[1:]):
                name = match.group(0)
                tail = design[1:][match.end():].lstrip()
                if name not in selected.columns and not tail.startswith("(") \
                        and name not in {"TRUE", "FALSE"}:
                    missing_design.append(name)
            if missing_design:
                errors.append(
                    "diffacc.design references missing sample-sheet columns: "
                    + ", ".join(sorted(set(missing_design)))
                )
            design_names = set(
                re.findall(r"[A-Za-z][A-Za-z0-9_.]*", design[1:])
            )
            if contrast_factor and contrast_factor not in design_names:
                errors.append(
                    f"diffacc contrast factor {contrast_factor!r} is absent from diffacc.design"
                )
            for name in design_names & set(selected.columns):
                if selected[name].astype(str).str.strip().eq("").any():
                    errors.append(f"diffacc design column {name!r} contains missing values")
        if "condition" in selected.columns:
            counts = selected.groupby("condition")["sample"].nunique()
            low = counts[counts < 2]
            if not low.empty:
                errors.append(
                    "differential analysis requires at least two biological replicates per condition; "
                    + ", ".join(f"{k}={v}" for k, v in low.items())
                )

    motif = _enabled(config.get("motif", {}))
    annotation = _enabled(config.get("annotation", {}))
    enrichment = _enabled(config.get("functional_enrichment", {}))
    footprinting = _enabled(config.get("footprinting", {}))
    chromvar = _enabled(config.get("chromvar", {}))
    filtering = config.setdefault("filtering", {})
    for key in (
        "remove_mito", "remove_blacklist", "keep_proper_pairs",
        "remove_duplicates", "tn5_shift",
    ):
        if key in filtering and not isinstance(filtering[key], bool):
            errors.append(f"filtering.{key} must be true or false (not a quoted string)")
    boolean_settings = (
        ("trimming", "enabled"),
        ("peaks", "run_genrich"),
        ("diffacc", "enabled"),
        ("annotation", "enabled"),
        ("functional_enrichment", "enabled"),
        ("functional_enrichment", "kegg"),
        ("motif", "enabled"),
        ("footprinting", "enabled"),
        ("chromvar", "enabled"),
        ("qc", "library_complexity"),
        ("qc", "fingerprint"),
        ("qc", "sample_similarity"),
        ("report", "enabled"),
    )
    for section_name, key in boolean_settings:
        section = config.get(section_name, {})
        if isinstance(section, dict) and key in section \
                and not isinstance(section[key], bool):
            errors.append(
                f"{section_name}.{key} must be true or false (not a quoted string)"
            )
    sample_similarity = config.get("qc", {}).get("sample_similarity", True)
    if sample_similarity is True and len(selected) < 2:
        errors.append(
            "qc.sample_similarity requires at least two included libraries; "
            "set it to false for a single-library QC run"
        )
    report_enabled = config.get("report", {}).get("enabled", True)
    condition_count = selected["condition"].nunique()
    if condition_count > 12 and (sample_similarity is True or report_enabled is True):
        errors.append(
            f"categorical QC/report overlays support at most 12 conditions, "
            f"but {condition_count} are included; consolidate conditions or set "
            "qc.sample_similarity and report.enabled to false"
        )
    tss_section = config.get("qc", {}).get("tss", {})
    if not isinstance(tss_section, dict):
        errors.append("qc.tss must be a mapping with enabled/max_regions fields")
        tss_section = {}
    if isinstance(tss_section, dict) and "enabled" in tss_section \
            and not isinstance(tss_section["enabled"], bool):
        errors.append("qc.tss.enabled must be true or false (not a quoted string)")
    fragment_section = config.get("qc", {}).get("fragment_sizes", {})
    if not isinstance(fragment_section, dict):
        errors.append("qc.fragment_sizes must be a mapping with an enabled field")
        fragment_section = {}
    if isinstance(fragment_section, dict) and "enabled" in fragment_section \
            and not isinstance(fragment_section["enabled"], bool):
        errors.append(
            "qc.fragment_sizes.enabled must be true or false (not a quoted string)"
        )

    _check_numeric(
        errors, "filtering.min_mapq", filtering.get("min_mapq", 30),
        integer=True, lower=0, upper=255,
    )
    _check_numeric(
        errors, "filtering.tn5_shift_genome_chunk_length",
        filtering.get("tn5_shift_genome_chunk_length", 50_000_000),
        integer=True, lower=1,
    )
    _check_numeric(
        errors, "sra.max_spots", config.get("sra", {}).get("max_spots", 0),
        integer=True, lower=0,
    )
    _check_numeric(
        errors, "qc.tss.max_regions", tss_section.get("max_regions", 0),
        integer=True, lower=0,
    )
    _check_numeric(
        errors, "peaks.macs3_qvalue",
        config.get("peaks", {}).get("macs3_qvalue", 0.01),
        lower=0, upper=1, lower_open=True, upper_open=True,
    )
    _check_numeric(
        errors, "diffacc.fdr", config.get("diffacc", {}).get("fdr", 0.05),
        lower=0, upper=1, lower_open=True, upper_open=True,
    )
    _check_numeric(
        errors, "diffacc.lfc_threshold",
        config.get("diffacc", {}).get("lfc_threshold", 0), lower=0,
    )
    _check_numeric(
        errors, "annotation.tss_upstream",
        config.get("annotation", {}).get("tss_upstream", 3000),
        integer=True, lower=1,
    )
    _check_numeric(
        errors, "annotation.tss_downstream",
        config.get("annotation", {}).get("tss_downstream", 3000),
        integer=True, lower=1,
    )
    _check_numeric(
        errors, "functional_enrichment.qvalue",
        config.get("functional_enrichment", {}).get("qvalue", 0.05),
        lower=0, upper=1, lower_open=True, upper_open=True,
    )
    _check_numeric(
        errors, "functional_enrichment.top_n",
        config.get("functional_enrichment", {}).get("top_n", 20),
        integer=True, lower=1,
    )
    _check_numeric(
        errors, "motif.homer_size", config.get("motif", {}).get("homer_size", 200),
        integer=True, lower=1,
    )
    _check_numeric(
        errors, "chromvar.top_n", config.get("chromvar", {}).get("top_n", 30),
        integer=True, lower=1,
    )
    _check_numeric(
        errors, "chromvar.seed", config.get("chromvar", {}).get("seed", 1),
        integer=True, lower=0, upper=2_147_483_647,
    )
    mito_contigs = filtering.get(
        "mitochondrial_contigs", DEFAULT_MITOCHONDRIAL_CONTIGS
    )
    if not isinstance(mito_contigs, list):
        errors.append("filtering.mitochondrial_contigs must be a YAML list")
        mito_contigs = []
    else:
        clean_contigs = []
        for contig in mito_contigs:
            if not isinstance(contig, str) or not SAFE_CONTIG.fullmatch(contig):
                errors.append(
                    "filtering.mitochondrial_contigs entries must be non-empty, "
                    "shell-safe contig names"
                )
                continue
            clean_contigs.append(contig)
        if len(clean_contigs) != len(set(clean_contigs)):
            errors.append("filtering.mitochondrial_contigs contains duplicate names")
        mito_contigs = clean_contigs
    if not mito_contigs:
        errors.append(
            "filtering.mitochondrial_contigs must be a non-empty YAML list"
        )
    # Resolve the fallback once so every downstream consumer sees the same list.
    filtering["mitochondrial_contigs"] = mito_contigs
    if filtering.get("keep_proper_pairs", True) is not True:
        errors.append(
            "filtering.keep_proper_pairs must be true: peak calling, fragment "
            "counting, FRiP, library complexity, and usable-fragment QC must use "
            "the same properly paired fragment universe"
        )
    if filtering.get("tn5_shift", True) is not True:
        errors.append(
            "filtering.tn5_shift must be true: the workflow publishes Tn5 cut-site "
            "tracks and insertion-based QC, which are invalid from unshifted reads"
        )
    if motif and not _enabled(diff):
        errors.append("motif.enabled requires diffacc.enabled because it uses up/down peak sets")
    if enrichment and not (_enabled(diff) and annotation):
        errors.append("functional_enrichment.enabled requires diffacc.enabled and annotation.enabled")
    if chromvar and len(selected) < 2:
        errors.append(
            "chromvar.enabled requires at least two included libraries for "
            "cross-sample variability and heatmap outputs"
        )
    if enrichment:
        enrichment_config = config.get("functional_enrichment", {})
        ontologies = enrichment_config.get("ontologies", ["BP", "MF"])
        if not isinstance(ontologies, list) or any(
            not isinstance(value, str) or value not in {"BP", "MF", "CC"}
            for value in ontologies
        ):
            errors.append(
                "functional_enrichment.ontologies must be a YAML list containing only BP, MF, or CC"
            )
        elif len(ontologies) != len(set(ontologies)):
            errors.append("functional_enrichment.ontologies contains duplicates")
        elif not ontologies and not enrichment_config.get("kegg", False):
            errors.append(
                "functional_enrichment.enabled requires at least one GO ontology or kegg=true"
            )
    if footprinting:
        motif_db = str(config.get("footprinting", {}).get("motif_db", "")).strip()
        if not motif_db:
            errors.append("footprinting.enabled requires footprinting.motif_db")
        elif not SAFE_PATH.fullmatch(motif_db):
            errors.append(
                "footprinting.motif_db contains whitespace or shell-special characters; "
                "use a simple path"
            )
        elif check_files and not Path(motif_db).is_file():
            errors.append(f"footprinting motif database does not exist: {motif_db}")
        if selected["condition"].nunique() < 2:
            errors.append("comparative footprinting requires at least two conditions")

    resources = config.get("resources", {})
    for key in ("align_threads", "sort_threads", "general_threads"):
        _check_numeric(
            errors, f"resources.{key}", resources.get(key, 0),
            integer=True, lower=1,
        )

    large_index = config.get("alignment", {}).get("large_index", False)
    if not isinstance(large_index, bool):
        errors.append("alignment.large_index must be true or false")

    for key in ("samples", "results_dir", "raw_dir", "processed_dir", "logs_dir", "reference_dir"):
        path = str(config.get(key, "")).strip()
        if path and not SAFE_PATH.fullmatch(path):
            errors.append(
                f"{key} contains whitespace or shell-special characters; use a simple path"
            )

    peaks = config.get("peaks", {})
    extra = str(peaks.get("macs3_extra", ""))
    if "--shift" in extra or "--extsize" in extra:
        errors.append("MACS3 BAMPE mode uses complete fragments; remove --shift/--extsize from peaks.macs3_extra")
    min_replicates = _check_numeric(
        errors, "peaks.consensus_min_replicates",
        peaks.get("consensus_min_replicates", 2), integer=True, lower=1,
    )
    _check_numeric(
        errors, "peaks.consensus_peak_width",
        peaks.get("consensus_peak_width", 500), integer=True, lower=50,
    )
    if min_replicates is not None and not selected.empty:
        available = selected.groupby("condition")["sample"].nunique()
        insufficient = available[available < min_replicates]
        if not insufficient.empty:
            detail = ", ".join(
                f"{condition}={count}" for condition, count in insufficient.items()
            )
            errors.append(
                "every condition must satisfy peaks.consensus_min_replicates="
                f"{min_replicates}; insufficient selected libraries: {detail}"
            )

    genome = config.get("genome", {})
    build = str(genome.get("build", ""))
    if build not in {"human", "mouse", "mouse_mm10", "rat", "custom"}:
        errors.append(
            "genome.build must be human, mouse, mouse_mm10, rat, or custom"
        )
    if filtering.get("remove_blacklist", True) and not genome.get("blacklist") \
            and build in {"mouse", "rat", "custom"}:
        errors.append(
            f"filtering.remove_blacklist=true but {build!r} has no bundled "
            "assembly-matched blacklist; provide genome.blacklist or set removal false"
        )
    for key in ("fasta", "gtf", "blacklist"):
        path = str(genome.get(key, "")).strip()
        if path and not SAFE_PATH.fullmatch(path):
            errors.append(
                f"genome.{key} contains whitespace or shell-special characters; use a simple path"
            )
        elif path and check_files:
            reference_path = Path(path)
            if not reference_path.is_file():
                errors.append(f"genome.{key} does not exist: {path}")
            elif key == "blacklist" and filtering.get("remove_blacklist", True) \
                    and reference_path.stat().st_size == 0:
                errors.append(
                    "genome.blacklist is empty while blacklist removal is enabled; "
                    "provide a BED with intervals or set filtering.remove_blacklist=false"
                )
    fasta = str(genome.get("fasta", "")).strip()
    if _is_gzip_reference(fasta, check_files):
        errors.append(
            "genome.fasta must be an uncompressed FASTA; decompress .fa.gz before use"
        )
    gtf = str(genome.get("gtf", "")).strip()
    if _is_gzip_reference(gtf, check_files):
        errors.append(
            "genome.gtf must be an uncompressed GTF; decompress .gtf.gz before use"
        )
    if build == "custom":
        custom = genome.get("custom", {})
        try:
            if int(custom.get("effective_genome_size", 0)) < 1:
                raise ValueError
        except (TypeError, ValueError):
            errors.append(
                "genome.custom.effective_genome_size must be a positive integer"
            )
        if not str(custom.get("macs_gsize", "")).strip():
            errors.append("genome.custom.macs_gsize is required for a custom genome")
        ensembl_fields = ("ensembl_species", "ensembl_release", "ensembl_assembly")
        complete_ensembl = all(str(custom.get(key, "")).strip() for key in ensembl_fields)
        try:
            complete_ensembl = complete_ensembl and int(custom.get("ensembl_release", 0)) > 0
        except (TypeError, ValueError):
            complete_ensembl = False
        for local_key in ("fasta", "gtf"):
            if not str(genome.get(local_key, "")).strip() and not complete_ensembl:
                errors.append(
                    f"custom genome needs genome.{local_key} or complete genome.custom "
                    "ensembl_species/ensembl_release/ensembl_assembly fields"
                )

        taxid = custom.get("taxid", 0)
        try:
            taxid = int(taxid)
        except (TypeError, ValueError):
            taxid = 0
        if enrichment:
            if taxid < 1:
                errors.append(
                    "functional_enrichment.enabled on a custom genome requires "
                    "a positive genome.custom.taxid"
                )
            if not str(custom.get("orgdb", "")).strip():
                errors.append(
                    "functional_enrichment.enabled on a custom genome requires "
                    "genome.custom.orgdb and its matching R package"
                )
            if config.get("functional_enrichment", {}).get("kegg", False) \
                    and taxid not in {9606, 10090, 10116}:
                errors.append(
                    "KEGG enrichment currently supports taxid 9606, 10090, or 10116; "
                    "disable functional_enrichment.kegg for this custom genome"
                )
        if chromvar:
            if taxid < 1:
                errors.append(
                    "chromvar.enabled on a custom genome requires a positive "
                    "genome.custom.taxid"
                )
            if not str(custom.get("bsgenome", "")).strip():
                errors.append(
                    "chromvar.enabled on a custom genome requires genome.custom.bsgenome "
                    "and its matching R package"
                )

    if errors:
        raise ValueError("Invalid ATAC-seq inputs:\n- " + "\n- ".join(errors))
    return selected
