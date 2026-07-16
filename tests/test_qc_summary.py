from pathlib import Path
import importlib.util

import pytest
import pysam


SCRIPT = Path(__file__).parents[1] / "workflow" / "scripts" / "qc_summary.py"
SPEC = importlib.util.spec_from_file_location("qc_summary", SCRIPT)
MOD = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MOD)


def test_qc_classification_is_explicit():
    status, notes = MOD.classify(
        {"frip": 0.15, "tss": 4.0, "mito_fraction": 0.6,
         "usable_fragments": 30_000_000, "peak_count": 50_000, "nrf": 0.9},
        "human",
    )
    assert status == "fail"
    assert "FRiP<0.20" in notes and "TSS<5" in notes


def test_nonhuman_tss_is_not_scored_against_human_thresholds():
    status, notes = MOD.classify(
        {"frip": 0.4, "tss": 2.0, "mito_fraction": 0.1,
         "usable_fragments": 20_000_000, "peak_count": 30_000, "nrf": 0.9},
        "mouse",
    )
    assert status == "pass"
    assert "TSS" not in notes


@pytest.mark.parametrize(
    "metrics,note",
    [
        ({"usable_fragments": 999_999, "peak_count": 30_000}, "usable_fragments<1M"),
        ({"usable_fragments": 20_000_000, "peak_count": 999}, "peak_count<1k"),
        ({"usable_fragments": 20_000_000, "peak_count": 30_000, "nrf": 0.49}, "NRF<0.50"),
    ],
)
def test_catastrophic_yield_or_complexity_is_a_failure(metrics, note):
    status, notes = MOD.classify(metrics, "custom")
    assert status == "fail"
    assert note in notes


def test_empty_metric_file_has_no_data_row(tmp_path):
    path = tmp_path / "empty.tsv"
    path.write_text("sample\tFRiP\n")
    assert MOD.read_data_row(path) == {}


def write_pairs(path, contigs):
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": name, "LN": 1000} for name in contigs],
    }
    with pysam.AlignmentFile(path, "wb", header=header) as bam:
        for index, contig in enumerate(contigs):
            rid = bam.get_tid(contig)
            start = 100 + 200 * index
            for is_read1 in (True, False):
                read = pysam.AlignedSegment()
                read.query_name = f"pair_{contig}"
                read.query_sequence = "A" * 50
                read.flag = 99 if is_read1 else 147
                read.reference_id = rid
                read.reference_start = start if is_read1 else start + 50
                read.mapping_quality = 60
                read.cigar = ((0, 50),)
                read.next_reference_id = rid
                read.next_reference_start = start + 50 if is_read1 else start
                read.template_length = 100 if is_read1 else -100
                read.query_qualities = pysam.qualitystring_to_array("I" * 50)
                bam.write(read)


def test_bam_counts_uses_configured_mitochondrial_contigs(tmp_path):
    raw = tmp_path / "raw.bam"
    filtered = tmp_path / "filtered.bam"
    write_pairs(raw, ["M", "NC_012920.1", "chr1"])
    write_pairs(filtered, ["chr1"])

    counts = MOD.bam_counts(
        raw, filtered, {"M", "NC_012920.1"}
    )

    assert counts["raw_fragments"] == 3
    assert counts["proper_fragments"] == 3
    assert counts["mito_fragments"] == 2
    assert counts["usable_fragments"] == 1
