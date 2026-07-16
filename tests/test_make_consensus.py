"""Scientific invariants for deterministic ATAC consensus peaks."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "make_consensus.py"
CHROM_SIZES = "chr1\t1000\nchr2\t1000\nchr10\t1000\n"


def narrowpeak(rows: list[tuple[str, int, int, int, float]]) -> str:
    """Return narrowPeak text from chrom/start/end/absolute-summit/signal rows."""

    lines = []
    for index, (chrom, start, end, summit, signal) in enumerate(rows, 1):
        offset = summit - start
        lines.append(
            f"{chrom}\t{start}\t{end}\tp{index}\t100\t.\t{signal}\t3\t2\t{offset}\n"
        )
    return "".join(lines)


def run_consensus(
    tmp_path: Path,
    peak_files: list[str],
    conditions: list[str],
    *,
    min_overlap: int = 2,
    min_option: str = "--min-overlap",
    width: int = 100,
    chrom_sizes: str = CHROM_SIZES,
    blacklist: str = "",
) -> tuple[str, str, str]:
    tmp_path.mkdir(parents=True, exist_ok=True)
    paths = []
    for index, content in enumerate(peak_files):
        path = tmp_path / f"rep{index}.narrowPeak"
        path.write_text(content)
        paths.append(str(path))

    chrom_path = tmp_path / "chrom.sizes"
    chrom_path.write_text(chrom_sizes)
    bed_path = tmp_path / "consensus.bed"
    saf_path = tmp_path / "consensus.saf"
    command = [
            sys.executable,
            str(SCRIPT),
            "--peaks",
            *paths,
            "--conditions",
            *conditions,
            "--chrom",
            str(chrom_path),
            min_option,
            str(min_overlap),
            "--width",
            str(width),
            "--bed",
            str(bed_path),
            "--saf",
            str(saf_path),
        ]
    if blacklist:
        blacklist_path = tmp_path / "blacklist.bed"
        blacklist_path.write_text(blacklist)
        command.extend(["--blacklist", str(blacklist_path)])
    proc = subprocess.run(
        command,
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 0, proc.stderr
    return bed_path.read_text(), saf_path.read_text(), proc.stderr


def test_uses_stronger_summit_and_fixed_width(tmp_path: Path) -> None:
    rep1 = narrowpeak([("chr1", 130, 260, 200, 5.0)])
    rep2 = narrowpeak([("chr1", 150, 280, 220, 10.0)])

    bed, _, _ = run_consensus(tmp_path, [rep1, rep2], ["control", "control"])

    assert bed.splitlines() == ["chr1\t170\t270\tpeak_1"]


def test_support_must_be_within_one_condition(tmp_path: Path) -> None:
    control = narrowpeak([("chr1", 130, 260, 200, 20.0)])
    treated = narrowpeak([("chr1", 150, 280, 220, 20.0)])

    bed, saf, _ = run_consensus(tmp_path, [control, treated], ["control", "treated"])

    assert bed == ""
    assert saf == "GeneID\tChr\tStart\tEnd\tStrand\n"


def test_either_condition_can_independently_support_a_peak(tmp_path: Path) -> None:
    files = [
        narrowpeak([("chr1", 130, 260, 200, 5.0)]),
        narrowpeak([("chr1", 150, 280, 220, 10.0)]),
        narrowpeak([("chr1", 630, 760, 700, 5.0)]),
        narrowpeak([("chr1", 650, 780, 720, 10.0)]),
    ]
    bed, _, _ = run_consensus(
        tmp_path,
        files,
        ["control", "control", "treated", "treated"],
    )

    assert bed.splitlines() == [
        "chr1\t170\t270\tpeak_1",
        "chr1\t670\t770\tpeak_2",
    ]


def test_transitive_chain_never_creates_a_long_union(tmp_path: Path) -> None:
    # Adjacent fixed windows overlap in a chain, but the first and third do not:
    # [50,150), [140,240), [230,330). Interval-component merging would produce
    # a misleading 280 bp feature. Direct support and NMS retain one 100 bp site.
    files = [
        narrowpeak([("chr1", 60, 140, 100, 5.0)]),
        narrowpeak([("chr1", 150, 230, 190, 10.0)]),
        narrowpeak([("chr1", 240, 320, 280, 5.0)]),
    ]

    bed, _, _ = run_consensus(tmp_path, files, ["A", "A", "A"])

    assert bed.splitlines() == ["chr1\t140\t240\tpeak_1"]
    start, end = map(int, bed.splitlines()[0].split("\t")[1:3])
    assert end - start == 100


def test_many_peaks_in_one_file_count_as_one_replicate(tmp_path: Path) -> None:
    one_replicate = narrowpeak(
        [
            ("chr1", 120, 240, 180, 10.0),
            ("chr1", 140, 260, 200, 20.0),
            ("chr1", 160, 280, 220, 30.0),
        ]
    )

    bed, _, _ = run_consensus(tmp_path, [one_replicate], ["control"])

    assert bed == ""


def test_output_is_identical_when_inputs_are_permuted(tmp_path: Path) -> None:
    first = narrowpeak([("chr1", 130, 260, 200, 5.0)])
    second = narrowpeak([("chr1", 150, 280, 220, 5.0)])

    original, _, _ = run_consensus(tmp_path / "original", [first, second], ["A", "A"])
    permuted, _, _ = run_consensus(tmp_path / "permuted", [second, first], ["A", "A"])

    assert original == permuted == "chr1\t150\t250\tpeak_1\n"


def test_windows_shift_inside_chromosome_bounds(tmp_path: Path) -> None:
    chrom_sizes = "chr1\t100\n"
    files = [
        narrowpeak(
            [
                ("chr1", 0, 20, 5, 10.0),
                ("chr1", 80, 100, 95, 10.0),
            ]
        ),
        narrowpeak(
            [
                ("chr1", 0, 20, 5, 10.0),
                ("chr1", 80, 100, 95, 10.0),
            ]
        ),
    ]

    bed, _, _ = run_consensus(
        tmp_path,
        files,
        ["A", "A"],
        width=40,
        chrom_sizes=chrom_sizes,
    )

    assert bed.splitlines() == [
        "chr1\t0\t40\tpeak_1",
        "chr1\t60\t100\tpeak_2",
    ]


def test_edge_shifted_windows_find_support_beyond_one_width(tmp_path: Path) -> None:
    chrom_sizes = "chr1\t10000\n"
    left_files = [
        narrowpeak([("chr1", 0, 20, 0, 10.0)]),
        narrowpeak([("chr1", 590, 620, 600, 10.0)]),
    ]
    right_files = [
        narrowpeak([("chr1", 9980, 10000, 9999, 10.0)]),
        narrowpeak([("chr1", 9380, 9410, 9400, 10.0)]),
    ]

    left, _, _ = run_consensus(
        tmp_path / "left", left_files, ["A", "A"],
        width=500, chrom_sizes=chrom_sizes,
    )
    right, _, _ = run_consensus(
        tmp_path / "right", right_files, ["A", "A"],
        width=500, chrom_sizes=chrom_sizes,
    )

    assert len(left.splitlines()) == 1
    assert len(right.splitlines()) == 1
    assert int(left.splitlines()[0].split("\t")[2]) - int(left.splitlines()[0].split("\t")[1]) == 500
    assert int(right.splitlines()[0].split("\t")[2]) - int(right.splitlines()[0].split("\t")[1]) == 500


def test_blacklist_removes_whole_fixed_width_candidates(tmp_path: Path) -> None:
    peaks = narrowpeak([
        ("chr1", 100, 220, 150, 10.0),
        ("chr1", 700, 820, 750, 10.0),
    ])
    bed, saf, stderr = run_consensus(
        tmp_path,
        [peaks, peaks],
        ["A", "A"],
        blacklist="chr1\t120\t130\n",
    )

    assert bed == "chr1\t700\t800\tpeak_1\n"
    assert saf.splitlines()[1] == "peak_1\tchr1\t701\t800\t+"
    assert "blacklisted replicate window(s) removed" in stderr


def test_blacklisted_replicate_cannot_support_a_clean_boundary_window(tmp_path: Path) -> None:
    clean = narrowpeak([("chr1", 0, 20, 0, 10.0)])       # window [0,500)
    blocked = narrowpeak([("chr1", 490, 520, 500, 10.0)])  # window [250,750)

    bed, _, _ = run_consensus(
        tmp_path,
        [clean, blocked],
        ["A", "A"],
        width=500,
        chrom_sizes="chr1\t10000\n",
        blacklist="chr1\t500\t600\n",
    )

    assert bed == ""


def test_saf_is_one_based_and_bed_is_zero_based(tmp_path: Path) -> None:
    peak = narrowpeak([("chr1", 100, 220, 150, 10.0)])

    bed, saf, _ = run_consensus(tmp_path, [peak, peak], ["A", "A"])

    assert bed == "chr1\t100\t200\tpeak_1\n"
    assert saf.splitlines()[1] == "peak_1\tchr1\t101\t200\t+"


def test_output_order_follows_chrom_sizes(tmp_path: Path) -> None:
    one = narrowpeak(
        [
            ("chr10", 100, 220, 150, 10.0),
            ("chr2", 100, 220, 150, 10.0),
        ]
    )

    bed, _, _ = run_consensus(tmp_path, [one, one], ["A", "A"])

    assert [line.split("\t")[0] for line in bed.splitlines()] == ["chr2", "chr10"]


def test_condition_count_must_match_peak_file_count(tmp_path: Path) -> None:
    peak = tmp_path / "rep.narrowPeak"
    peak.write_text(narrowpeak([("chr1", 100, 220, 150, 10.0)]))
    chrom = tmp_path / "chrom.sizes"
    chrom.write_text(CHROM_SIZES)

    proc = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--peaks",
            str(peak),
            str(peak),
            "--conditions",
            "A",
            "--chrom",
            str(chrom),
            "--bed",
            str(tmp_path / "out.bed"),
            "--saf",
            str(tmp_path / "out.saf"),
        ],
        capture_output=True,
        text=True,
    )

    assert proc.returncode != 0
    assert "--conditions has 1 values but --peaks has 2 files" in proc.stderr


def test_min_replicates_is_a_supported_cli_alias(tmp_path: Path) -> None:
    peak = narrowpeak([("chr1", 100, 220, 150, 10.0)])

    bed, _, _ = run_consensus(
        tmp_path,
        [peak, peak],
        ["A", "A"],
        min_option="--min-replicates",
    )

    assert bed == "chr1\t100\t200\tpeak_1\n"
