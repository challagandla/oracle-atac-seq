import csv
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "workflow/scripts/merge_featurecounts.py"


def write_library(
    root: Path,
    sample: str,
    rows: list[tuple[str, str, int, int, str, int, int]],
    assigned: int | None = None,
) -> tuple[Path, Path]:
    counts = root / f"{sample}.counts.tsv"
    summary = root / f"{sample}.counts.tsv.summary"
    bam = root / f"{sample}.filtered.bam"
    with counts.open("w", encoding="utf-8", newline="") as handle:
        handle.write("# featureCounts test output\n")
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["Geneid", "Chr", "Start", "End", "Strand", "Length", bam])
        writer.writerows(rows)
    total = sum(row[-1] for row in rows) if assigned is None else assigned
    with summary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["Status", bam])
        writer.writerow(["Assigned", total])
        writer.writerow(["Unassigned_NoFeatures", 7])
    return counts, summary


def run_merge(
    tmp_path: Path,
    libraries: list[tuple[str, Path, Path]],
) -> subprocess.CompletedProcess[str]:
    out_counts = tmp_path / "merged.tsv"
    out_summary = tmp_path / "merged.tsv.summary"
    return subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--counts",
            *(str(counts) for _, counts, _ in libraries),
            "--summaries",
            *(str(summary) for _, _, summary in libraries),
            "--samples",
            *(sample for sample, _, _ in libraries),
            "--out-counts",
            str(out_counts),
            "--out-summary",
            str(out_summary),
        ],
        text=True,
        capture_output=True,
        check=False,
    )


def test_merges_valid_libraries_with_stable_sample_columns(tmp_path):
    metadata = [
        ("peak_1", "chr1", 1, 100, ".", 100),
        ("peak_2", "chr2", 201, 250, ".", 50),
    ]
    libraries = []
    for sample, values in (("control-1", (3, 4)), ("case.filtered.bam", (5, 6))):
        rows = [(*feature, value) for feature, value in zip(metadata, values, strict=True)]
        counts, summary = write_library(tmp_path, sample, rows)
        libraries.append((sample, counts, summary))

    result = run_merge(tmp_path, libraries)

    assert result.returncode == 0, result.stderr
    lines = (tmp_path / "merged.tsv").read_text(encoding="utf-8").splitlines()
    assert lines[0].startswith("# Validated merge")
    assert lines[1].split("\t") == [
        "Geneid", "Chr", "Start", "End", "Strand", "Length",
        "control-1", "case.filtered.bam",
    ]
    assert lines[2].split("\t")[-2:] == ["3", "5"]
    summary = (tmp_path / "merged.tsv.summary").read_text(encoding="utf-8")
    assert "Status\tcontrol-1\tcase.filtered.bam\n" in summary
    assert "Assigned\t7\t11\n" in summary


def test_rejects_feature_metadata_drift_without_publishing_outputs(tmp_path):
    first = [("peak_1", "chr1", 1, 100, ".", 100, 3)]
    second = [("peak_1", "chr1", 2, 101, ".", 100, 5)]
    c1, s1 = write_library(tmp_path, "control_1", first)
    c2, s2 = write_library(tmp_path, "treated_1", second)

    result = run_merge(
        tmp_path,
        [("control_1", c1, s1), ("treated_1", c2, s2)],
    )

    assert result.returncode != 0
    assert "feature metadata or row order differs" in result.stderr
    assert not (tmp_path / "merged.tsv").exists()
    assert not (tmp_path / "merged.tsv.summary").exists()


def test_rejects_count_and_assigned_summary_disagreement(tmp_path):
    rows = [("peak_1", "chr1", 1, 100, ".", 100, 3)]
    counts, summary = write_library(tmp_path, "control_1", rows, assigned=4)

    result = run_merge(tmp_path, [("control_1", counts, summary)])

    assert result.returncode != 0
    assert "matrix sum 3 differs from Assigned summary 4" in result.stderr
    assert not (tmp_path / "merged.tsv").exists()


def test_rejects_zero_assigned_library(tmp_path):
    rows = [("peak_1", "chr1", 1, 100, ".", 100, 0)]
    counts, summary = write_library(tmp_path, "control_1", rows, assigned=0)

    result = run_merge(tmp_path, [("control_1", counts, summary)])

    assert result.returncode != 0
    assert "zero fragments assigned" in result.stderr


def test_rejects_count_file_with_the_wrong_sample_identity(tmp_path):
    rows = [("peak_1", "chr1", 1, 100, ".", 100, 3)]
    counts, summary = write_library(tmp_path, "actual", rows)

    result = run_merge(tmp_path, [("expected", counts, summary)])

    assert result.returncode != 0
    assert "expected one count column for sample 'expected'" in result.stderr


def test_rejects_duplicate_feature_identifiers(tmp_path):
    rows = [
        ("peak_1", "chr1", 1, 100, ".", 100, 3),
        ("peak_1", "chr1", 101, 200, ".", 100, 4),
    ]
    counts, summary = write_library(tmp_path, "control_1", rows)

    result = run_merge(tmp_path, [("control_1", counts, summary)])

    assert result.returncode != 0
    assert "duplicate feature identifier" in result.stderr


def test_merges_summary_rows_by_status_name_not_position(tmp_path):
    rows = [("peak_1", "chr1", 1, 100, ".", 100, 3)]
    c1, s1 = write_library(tmp_path, "control_1", rows)
    c2, s2 = write_library(tmp_path, "treated_1", rows)
    lines = s2.read_text(encoding="utf-8").splitlines()
    s2.write_text("\n".join([lines[0], lines[2], lines[1]]) + "\n", encoding="utf-8")

    result = run_merge(
        tmp_path,
        [("control_1", c1, s1), ("treated_1", c2, s2)],
    )

    assert result.returncode == 0, result.stderr
    summary = (tmp_path / "merged.tsv.summary").read_text(encoding="utf-8")
    assert "Assigned\t3\t3\n" in summary


def test_rejects_unequal_feature_row_counts(tmp_path):
    first = [
        ("peak_1", "chr1", 1, 100, ".", 100, 3),
        ("peak_2", "chr1", 101, 200, ".", 100, 4),
    ]
    second = [("peak_1", "chr1", 1, 100, ".", 100, 5)]
    c1, s1 = write_library(tmp_path, "control_1", first)
    c2, s2 = write_library(tmp_path, "treated_1", second)

    result = run_merge(
        tmp_path,
        [("control_1", c1, s1), ("treated_1", c2, s2)],
    )

    assert result.returncode != 0
    assert "different row counts" in result.stderr
    assert not (tmp_path / "merged.tsv").exists()


def test_rejects_empty_feature_table(tmp_path):
    counts, summary = write_library(tmp_path, "control_1", [], assigned=1)

    result = run_merge(tmp_path, [("control_1", counts, summary)])

    assert result.returncode != 0
    assert "contain no features" in result.stderr
    assert not (tmp_path / "merged.tsv").exists()


def test_rejects_mismatched_summary_status_sets(tmp_path):
    rows = [("peak_1", "chr1", 1, 100, ".", 100, 3)]
    c1, s1 = write_library(tmp_path, "control_1", rows)
    c2, s2 = write_library(tmp_path, "treated_1", rows)
    text = s2.read_text(encoding="utf-8")
    s2.write_text(
        text.replace("Unassigned_NoFeatures", "Unassigned_MultiMapping"),
        encoding="utf-8",
    )

    result = run_merge(
        tmp_path,
        [("control_1", c1, s1), ("treated_1", c2, s2)],
    )

    assert result.returncode != 0
    assert "different status rows" in result.stderr
    assert not (tmp_path / "merged.tsv.summary").exists()


def test_rejects_duplicate_summary_status(tmp_path):
    rows = [("peak_1", "chr1", 1, 100, ".", 100, 3)]
    counts, summary = write_library(tmp_path, "control_1", rows)
    with summary.open("a", encoding="utf-8") as handle:
        handle.write("Assigned\t3\n")

    result = run_merge(tmp_path, [("control_1", counts, summary)])

    assert result.returncode != 0
    assert "duplicate status 'Assigned'" in result.stderr
    assert not (tmp_path / "merged.tsv.summary").exists()
