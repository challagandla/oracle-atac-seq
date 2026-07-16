import gzip
import importlib.util
from pathlib import Path

import pandas as pd
import pytest


SCRIPT = Path(__file__).parents[1] / "workflow" / "scripts" / "validate_inputs.py"
SPEC = importlib.util.spec_from_file_location("validate_inputs", SCRIPT)
MOD = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MOD)


def base_config():
    return {
        "diffacc": {"enabled": True, "design": "~condition", "contrast": ["condition", "treat", "ctrl"]},
        "motif": {"enabled": False},
        "annotation": {"enabled": False},
        "functional_enrichment": {"enabled": False},
        "footprinting": {"enabled": False, "motif_db": ""},
        "peaks": {"macs3_extra": "--keep-dup all", "consensus_min_replicates": 2, "consensus_peak_width": 500},
        "genome": {"build": "human"},
        "resources": {"align_threads": 2, "sort_threads": 2, "general_threads": 2},
    }


def sample_table(tmp_path):
    rows = []
    for condition in ("ctrl", "treat"):
        for replicate in (1, 2):
            sample = f"{condition}_{replicate}"
            fq1 = tmp_path / f"{sample}_R1.fq.gz"
            fq2 = tmp_path / f"{sample}_R2.fq.gz"
            fq1.touch()
            fq2.touch()
            rows.append([sample, condition, str(replicate), str(fq1), str(fq2), ""])
    return pd.DataFrame(rows, columns=["sample", "condition", "replicate", "fq1", "fq2", "sra"])


def test_valid_table_and_reviewed_exclusion(tmp_path):
    table = sample_table(tmp_path)
    table["include"] = ["yes", "yes", "true", "1"]
    selected = MOD.validate_config_and_samples(base_config(), table)
    assert list(selected["sample"]) == ["ctrl_1", "ctrl_2", "treat_1", "treat_2"]


@pytest.mark.parametrize(
    "mutator,match",
    [
        (lambda x: x.assign(sample=["same", "same", "treat_1", "treat_2"]), "not unique"),
        (lambda x: x.assign(fq2=["", *x.fq2.iloc[1:]]), "needs both fq1 and fq2"),
        (lambda x: x.assign(sample=["bad name", *x["sample"].iloc[1:]]), "must use only"),
    ],
)
def test_invalid_sample_contract(tmp_path, mutator, match):
    with pytest.raises(ValueError, match=match):
        MOD.validate_config_and_samples(base_config(), mutator(sample_table(tmp_path)))


def test_rejects_reused_sra_accession_as_pseudoreplication(tmp_path):
    table = sample_table(tmp_path)
    table.loc[[0, 1], ["fq1", "fq2"]] = ""
    table.loc[0, "sra"] = "SRR123456"
    table.loc[1, "sra"] = "srr123456"
    with pytest.raises(ValueError, match="SRA run SRR123456 is reused"):
        MOD.validate_config_and_samples(base_config(), table, check_files=False)


def test_normalizes_valid_sra_accessions_to_canonical_uppercase(tmp_path):
    table = sample_table(tmp_path)
    table.loc[0, ["fq1", "fq2"]] = ""
    table.loc[0, "sra"] = "srr123456"
    selected = MOD.validate_config_and_samples(base_config(), table)
    assert selected.loc[selected["sample"] == "ctrl_1", "sra"].item() == "SRR123456"


def test_rejects_reused_fastq_pair_through_path_aliases(tmp_path):
    table = sample_table(tmp_path)
    table.loc[1, "fq1"] = str(Path(table.loc[0, "fq1"]).parent / "." / Path(table.loc[0, "fq1"]).name)
    table.loc[1, "fq2"] = str(Path(table.loc[0, "fq2"]).parent / "." / Path(table.loc[0, "fq2"]).name)
    with pytest.raises(ValueError, match="same paired FASTQ input is reused"):
        MOD.validate_config_and_samples(base_config(), table)


def test_rejects_one_fastq_mate_reused_with_a_different_pair(tmp_path):
    table = sample_table(tmp_path)
    table.loc[1, "fq1"] = table.loc[0, "fq1"]
    with pytest.raises(ValueError, match="FASTQ input file is reused"):
        MOD.validate_config_and_samples(base_config(), table)


def test_allows_duplicate_input_only_when_a_documented_row_is_excluded(tmp_path):
    table = sample_table(tmp_path)
    audit_row = table.iloc[[0]].copy()
    audit_row["sample"] = "ctrl_excluded"
    audit_row["replicate"] = "excluded"
    table = pd.concat([table, audit_row], ignore_index=True)
    table["include"] = [True, True, True, True, False]
    selected = MOD.validate_config_and_samples(base_config(), table)
    assert list(selected["sample"]) == ["ctrl_1", "ctrl_2", "treat_1", "treat_2"]


def test_excluded_audit_row_does_not_require_its_fastqs_to_still_exist(tmp_path):
    table = sample_table(tmp_path)
    audit_row = table.iloc[[0]].copy()
    audit_row["sample"] = "failed_library"
    audit_row["replicate"] = "failed"
    audit_row["fq1"] = str(tmp_path / "archived_R1.fastq.gz")
    audit_row["fq2"] = str(tmp_path / "archived_R2.fastq.gz")
    table = pd.concat([table, audit_row], ignore_index=True)
    table["include"] = [True, True, True, True, False]
    selected = MOD.validate_config_and_samples(base_config(), table)
    assert len(selected) == 4


def test_rejects_unreplicated_requested_contrast_factor(tmp_path):
    table = sample_table(tmp_path)
    table["batch"] = ["batch1", "batch2", "batch1", "batch1"]
    config = base_config()
    config["diffacc"]["design"] = "~batch"
    config["diffacc"]["contrast"] = ["batch", "batch2", "batch1"]
    with pytest.raises(ValueError, match="at least two independent libraries"):
        MOD.validate_config_and_samples(config, table)


def test_requires_replication_and_valid_cross_feature_config(tmp_path):
    table = sample_table(tmp_path).iloc[[0, 2]].copy()
    config = base_config()
    config["motif"]["enabled"] = True
    config["diffacc"]["enabled"] = False
    with pytest.raises(ValueError, match="motif.enabled"):
        MOD.validate_config_and_samples(config, table)


def test_rejects_single_replicate_differential_design(tmp_path):
    table = sample_table(tmp_path).iloc[[0, 2]].copy()
    with pytest.raises(ValueError, match="two biological replicates"):
        MOD.validate_config_and_samples(base_config(), table)


@pytest.mark.parametrize("name", ["ctrl_1_R1.txt", "ctrl_1_R1.FASTQ.GZ"])
def test_rejects_local_reads_without_a_supported_fastq_suffix(tmp_path, name):
    table = sample_table(tmp_path)
    table.loc[0, "fq1"] = str(tmp_path / name)

    with pytest.raises(ValueError, match=r"fq1 must end in \.fastq"):
        MOD.validate_config_and_samples(base_config(), table, check_files=False)


def test_sample_similarity_requires_two_included_libraries(tmp_path):
    table = sample_table(tmp_path).iloc[[0]].copy()
    config = base_config()
    config["diffacc"]["enabled"] = False
    config["peaks"]["consensus_min_replicates"] = 1
    config["qc"] = {"sample_similarity": True}
    with pytest.raises(ValueError, match="at least two included libraries"):
        MOD.validate_config_and_samples(config, table)

    config["qc"]["sample_similarity"] = False
    selected = MOD.validate_config_and_samples(config, table)
    assert list(selected["sample"]) == ["ctrl_1"]


def test_chromvar_requires_two_included_libraries(tmp_path):
    table = sample_table(tmp_path).iloc[[0]].copy()
    config = base_config()
    config["diffacc"]["enabled"] = False
    config["peaks"]["consensus_min_replicates"] = 1
    config["qc"] = {"sample_similarity": False}
    config["report"] = {"enabled": False}
    config["chromvar"] = {"enabled": True}
    with pytest.raises(ValueError, match="chromvar.enabled requires at least two"):
        MOD.validate_config_and_samples(config, table)


def test_categorical_reports_reject_more_than_twelve_conditions(tmp_path):
    rows = []
    for index in range(13):
        fq1 = tmp_path / f"c{index}_R1.fq.gz"
        fq2 = tmp_path / f"c{index}_R2.fq.gz"
        fq1.touch()
        fq2.touch()
        rows.append([f"sample_{index}", f"condition_{index}", "1",
                     str(fq1), str(fq2), ""])
    table = pd.DataFrame(
        rows, columns=["sample", "condition", "replicate", "fq1", "fq2", "sra"]
    )
    config = base_config()
    config["diffacc"]["enabled"] = False
    config["peaks"]["consensus_min_replicates"] = 1
    config["qc"] = {"sample_similarity": False}
    config["report"] = {"enabled": True}
    with pytest.raises(ValueError, match="at most 12 conditions"):
        MOD.validate_config_and_samples(config, table)

    config["report"]["enabled"] = False
    selected = MOD.validate_config_and_samples(config, table)
    assert len(selected) == 13


def test_rejects_single_end_peak_flags_in_bampe_mode(tmp_path):
    config = base_config()
    config["peaks"]["macs3_extra"] = "--shift -75 --extsize 150"
    with pytest.raises(ValueError, match="BAMPE"):
        MOD.validate_config_and_samples(config, sample_table(tmp_path))


def test_rejects_missing_design_covariate(tmp_path):
    config = base_config()
    config["diffacc"]["design"] = "~batch + condition"
    with pytest.raises(ValueError, match="missing sample-sheet columns: batch"):
        MOD.validate_config_and_samples(config, sample_table(tmp_path))


def test_rejects_shell_unsafe_input_path(tmp_path):
    table = sample_table(tmp_path)
    table.loc[0, "fq1"] = str(tmp_path / "unsafe path;reads_R1.fq.gz")
    with pytest.raises(ValueError, match="shell-special"):
        MOD.validate_config_and_samples(base_config(), table, check_files=False)


def test_requires_explicit_blacklist_decision_for_build_without_a_preset(tmp_path):
    config = base_config()
    config["filtering"] = {"remove_blacklist": True, "tn5_shift": True}
    config["genome"] = {"build": "mouse", "blacklist": ""}
    with pytest.raises(ValueError, match="no bundled assembly-matched blacklist"):
        MOD.validate_config_and_samples(config, sample_table(tmp_path))


def test_rejects_gzipped_fasta_and_unshifted_tss_qc(tmp_path):
    config = base_config()
    config["genome"]["fasta"] = "reference.fa.gz"
    config["filtering"] = {"remove_blacklist": True, "tn5_shift": False}
    with pytest.raises(ValueError) as error:
        MOD.validate_config_and_samples(config, sample_table(tmp_path), check_files=False)
    assert "uncompressed FASTA" in str(error.value)
    assert "filtering.tn5_shift must be true" in str(error.value)


def test_rejects_gzipped_gtf_with_an_actionable_error(tmp_path):
    config = base_config()
    config["genome"]["gtf"] = "genes.gtf.gz"
    with pytest.raises(ValueError, match="genome.gtf must be an uncompressed GTF"):
        MOD.validate_config_and_samples(config, sample_table(tmp_path), check_files=False)


@pytest.mark.parametrize("name", ["genes.gtf.GZ", "genes.gtf.bgz"])
def test_rejects_case_and_block_gzip_gtf_suffixes(tmp_path, name):
    config = base_config()
    config["genome"]["gtf"] = name
    with pytest.raises(ValueError, match="genome.gtf must be an uncompressed GTF"):
        MOD.validate_config_and_samples(config, sample_table(tmp_path), check_files=False)


def test_rejects_gzip_magic_bytes_hidden_behind_plain_gtf_suffix(tmp_path):
    compressed = tmp_path / "genes.gtf"
    with gzip.open(compressed, "wt") as handle:
        handle.write('chr1\ttest\tgene\t1\t10\t.\t+\t.\tgene_id "g1";\n')
    config = base_config()
    config["genome"]["gtf"] = str(compressed)
    with pytest.raises(ValueError, match="genome.gtf must be an uncompressed GTF"):
        MOD.validate_config_and_samples(config, sample_table(tmp_path))


def test_rejects_unshifted_run_even_when_tss_qc_is_disabled(tmp_path):
    config = base_config()
    config["filtering"] = {
        "remove_blacklist": True,
        "tn5_shift": False,
    }
    config["qc"] = {"tss": {"enabled": False}}
    with pytest.raises(ValueError, match="filtering.tn5_shift must be true"):
        MOD.validate_config_and_samples(config, sample_table(tmp_path))


def test_rejects_nonproper_pair_mode(tmp_path):
    config = base_config()
    config["filtering"] = {
        "remove_blacklist": True,
        "keep_proper_pairs": False,
        "tn5_shift": True,
    }
    with pytest.raises(ValueError, match="filtering.keep_proper_pairs must be true"):
        MOD.validate_config_and_samples(config, sample_table(tmp_path))


@pytest.mark.parametrize(
    "key",
    [
        "remove_mito", "remove_blacklist", "keep_proper_pairs",
        "remove_duplicates", "tn5_shift",
    ],
)
def test_filter_switches_reject_quoted_booleans(tmp_path, key):
    config = base_config()
    config["filtering"] = {
        "remove_mito": True,
        "remove_blacklist": True,
        "keep_proper_pairs": True,
        "remove_duplicates": True,
        "tn5_shift": True,
    }
    config["filtering"][key] = "false"
    with pytest.raises(ValueError, match=rf"filtering\.{key} must be true or false"):
        MOD.validate_config_and_samples(config, sample_table(tmp_path))


@pytest.mark.parametrize(
    "section,key",
    [
        ("trimming", "enabled"),
        ("peaks", "run_genrich"),
        ("diffacc", "enabled"),
        ("annotation", "enabled"),
        ("functional_enrichment", "kegg"),
        ("footprinting", "enabled"),
        ("chromvar", "enabled"),
    ],
)
def test_feature_switches_reject_quoted_booleans(tmp_path, section, key):
    config = base_config()
    config.setdefault(section, {})[key] = "false"
    with pytest.raises(ValueError, match=rf"{section}\.{key} must be true or false"):
        MOD.validate_config_and_samples(config, sample_table(tmp_path))


@pytest.mark.parametrize(
    "section,key",
    [
        ("functional_enrichment", "enabled"),
        ("motif", "enabled"),
        ("qc", "library_complexity"),
        ("qc", "fingerprint"),
        ("qc", "sample_similarity"),
        ("report", "enabled"),
    ],
)
def test_remaining_feature_switches_reject_quoted_booleans(
    tmp_path, section, key
):
    config = base_config()
    config.setdefault(section, {})[key] = "false"
    with pytest.raises(ValueError, match=rf"{section}\.{key} must be true or false"):
        MOD.validate_config_and_samples(config, sample_table(tmp_path))


@pytest.mark.parametrize("subsection", ["tss", "fragment_sizes"])
def test_nested_qc_switches_reject_quoted_booleans(tmp_path, subsection):
    config = base_config()
    config.setdefault("qc", {})[subsection] = {"enabled": "false"}
    with pytest.raises(
        ValueError, match=rf"qc\.{subsection}\.enabled must be true or false"
    ):
        MOD.validate_config_and_samples(config, sample_table(tmp_path))


def set_config_value(config, path, value):
    section = config
    for key in path[:-1]:
        section = section.setdefault(key, {})
    section[path[-1]] = value


@pytest.mark.parametrize(
    "path,value,label",
    [
        (("filtering", "min_mapq"), -1, "filtering.min_mapq"),
        (("filtering", "min_mapq"), 256, "filtering.min_mapq"),
        (("filtering", "min_mapq"), 30.5, "filtering.min_mapq"),
        (("filtering", "min_mapq"), 30.0, "filtering.min_mapq"),
        (("filtering", "min_mapq"), True, "filtering.min_mapq"),
        (("filtering", "tn5_shift_genome_chunk_length"), 0,
         "filtering.tn5_shift_genome_chunk_length"),
        (("sra", "max_spots"), -1, "sra.max_spots"),
        (("qc", "tss", "max_regions"), -1, "qc.tss.max_regions"),
        (("qc", "tss", "max_regions"), 100.0, "qc.tss.max_regions"),
        (("peaks", "macs3_qvalue"), 0, "peaks.macs3_qvalue"),
        (("peaks", "macs3_qvalue"), 1, "peaks.macs3_qvalue"),
        (("diffacc", "fdr"), 0, "diffacc.fdr"),
        (("diffacc", "fdr"), 1, "diffacc.fdr"),
        (("diffacc", "lfc_threshold"), -0.01, "diffacc.lfc_threshold"),
        (("annotation", "tss_upstream"), 0, "annotation.tss_upstream"),
        (("annotation", "tss_downstream"), 0, "annotation.tss_downstream"),
        (("functional_enrichment", "qvalue"), 0,
         "functional_enrichment.qvalue"),
        (("functional_enrichment", "qvalue"), 1,
         "functional_enrichment.qvalue"),
        (("functional_enrichment", "top_n"), 0,
         "functional_enrichment.top_n"),
        (("motif", "homer_size"), 0, "motif.homer_size"),
        (("chromvar", "top_n"), 0, "chromvar.top_n"),
        (("chromvar", "seed"), -1, "chromvar.seed"),
        (("chromvar", "seed"), 2_147_483_648, "chromvar.seed"),
        (("chromvar", "seed"), 1.5, "chromvar.seed"),
        (("diffacc", "fdr"), float("nan"), "diffacc.fdr"),
        (("diffacc", "lfc_threshold"), float("inf"),
         "diffacc.lfc_threshold"),
    ],
)
def test_rejects_out_of_range_or_non_finite_numeric_settings(
    tmp_path, path, value, label
):
    config = base_config()
    set_config_value(config, path, value)
    with pytest.raises(ValueError) as error:
        MOD.validate_config_and_samples(config, sample_table(tmp_path))
    assert label in str(error.value)


def test_accepts_documented_numeric_boundaries(tmp_path):
    config = base_config()
    settings = {
        ("filtering", "min_mapq"): 255,
        ("filtering", "tn5_shift_genome_chunk_length"): 1,
        ("sra", "max_spots"): 0,
        ("qc", "tss", "max_regions"): 0,
        ("peaks", "macs3_qvalue"): 1e-12,
        ("diffacc", "fdr"): 1 - 1e-12,
        ("diffacc", "lfc_threshold"): 0,
        ("annotation", "tss_upstream"): 1,
        ("annotation", "tss_downstream"): 1,
        ("functional_enrichment", "qvalue"): 0.5,
        ("functional_enrichment", "top_n"): 1,
        ("motif", "homer_size"): 1,
        ("chromvar", "top_n"): 1,
        ("chromvar", "seed"): 2_147_483_647,
    }
    for path, value in settings.items():
        set_config_value(config, path, value)

    selected = MOD.validate_config_and_samples(config, sample_table(tmp_path))
    assert len(selected) == 4


def test_rejects_consensus_replication_that_any_condition_cannot_meet(tmp_path):
    config = base_config()
    config["peaks"]["consensus_min_replicates"] = 3
    with pytest.raises(ValueError, match="every condition must satisfy"):
        MOD.validate_config_and_samples(config, sample_table(tmp_path))


@pytest.mark.parametrize(
    "value,match",
    [
        ("chrM", "must be a YAML list"),
        ([], "non-empty YAML list"),
        (["chrM", "chrM"], "duplicate"),
        (["chrM;touch_bad"], "shell-safe contig names"),
    ],
)
def test_rejects_invalid_mitochondrial_contig_contract(tmp_path, value, match):
    config = base_config()
    config["filtering"] = {
        "remove_mito": True,
        "remove_blacklist": True,
        "tn5_shift": True,
        "mitochondrial_contigs": value,
    }
    with pytest.raises(ValueError, match=match):
        MOD.validate_config_and_samples(config, sample_table(tmp_path))


def test_accepts_custom_mitochondrial_accession(tmp_path):
    config = base_config()
    config["filtering"] = {
        "remove_mito": True,
        "remove_blacklist": True,
        "tn5_shift": True,
        "mitochondrial_contigs": ["NC_012920.1"],
    }
    selected = MOD.validate_config_and_samples(config, sample_table(tmp_path))
    assert len(selected) == 4
    assert config["filtering"]["mitochondrial_contigs"] == ["NC_012920.1"]


def test_rejects_empty_blacklist_when_removal_is_enabled(tmp_path):
    blacklist = tmp_path / "blacklist.bed"
    blacklist.touch()
    config = base_config()
    config["filtering"] = {"remove_blacklist": True, "tn5_shift": True}
    config["genome"] = {"build": "human", "blacklist": str(blacklist)}

    with pytest.raises(ValueError, match="genome.blacklist is empty"):
        MOD.validate_config_and_samples(config, sample_table(tmp_path))


def custom_config():
    config = base_config()
    config["filtering"] = {"remove_blacklist": False, "tn5_shift": True}
    config["genome"] = {
        "build": "custom",
        "fasta": "",
        "gtf": "",
        "blacklist": "",
        "custom": {
            "effective_genome_size": 1000000,
            "macs_gsize": "1e6",
            "ensembl_species": "",
            "ensembl_release": 111,
            "ensembl_assembly": "",
            "taxid": 0,
            "txdb": "",
            "orgdb": "",
            "bsgenome": "",
        },
    }
    return config


def test_custom_genome_requires_local_references_or_complete_ensembl_tuple(tmp_path):
    with pytest.raises(ValueError) as error:
        MOD.validate_config_and_samples(
            custom_config(), sample_table(tmp_path), check_files=False
        )
    message = str(error.value)
    assert "genome.fasta or complete" in message
    assert "genome.gtf or complete" in message


def test_custom_optional_biology_never_defaults_to_human(tmp_path):
    config = custom_config()
    config["genome"]["custom"].update({
        "ensembl_species": "example_species",
        "ensembl_assembly": "Example1",
    })
    config["annotation"]["enabled"] = True
    config["functional_enrichment"]["enabled"] = True
    config["chromvar"] = {"enabled": True}
    with pytest.raises(ValueError) as error:
        MOD.validate_config_and_samples(config, sample_table(tmp_path), check_files=False)
    message = str(error.value)
    assert "positive genome.custom.taxid" in message
    assert "genome.custom.orgdb" in message
    assert "genome.custom.bsgenome" in message


def test_valid_custom_genome_with_local_references(tmp_path):
    config = custom_config()
    fasta = tmp_path / "genome.fa"
    gtf = tmp_path / "genes.gtf"
    fasta.touch()
    gtf.touch()
    config["genome"]["fasta"] = str(fasta)
    config["genome"]["gtf"] = str(gtf)
    selected = MOD.validate_config_and_samples(config, sample_table(tmp_path))
    assert len(selected) == 4
