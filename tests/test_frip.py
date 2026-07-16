"""Fragment-level FRiP semantics."""

from __future__ import annotations

import gzip
import importlib.util
import sys
from pathlib import Path

import pysam


SCRIPT = Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "frip.py"
SPEC = importlib.util.spec_from_file_location("oracle_frip", SCRIPT)
assert SPEC and SPEC.loader
FRIP = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = FRIP
SPEC.loader.exec_module(FRIP)


def _alignment(name: str, start: int, mate_start: int, template_length: int, read1: bool):
    read = pysam.AlignedSegment()
    read.query_name = name
    # Canonical inward-facing proper pair: 99/147.
    read.flag = 99 if read1 else 147
    read.reference_id = 0
    read.reference_start = start
    read.mapping_quality = 60
    read.cigar = ((0, 50),)
    read.next_reference_id = 0
    read.next_reference_start = mate_start
    read.template_length = template_length
    read.query_sequence = "A" * 50
    read.query_qualities = pysam.qualitystring_to_array("I" * 50)
    return read


def write_bam(tmp_path: Path, pairs: list[tuple[str, int, int]]) -> Path:
    """Write name/start1/start2 pairs as a coordinate-sorted paired BAM."""

    path = tmp_path / "reads.bam"
    reads = []
    for name, start1, start2 in pairs:
        fragment_length = start2 + 50 - start1
        reads.extend(
            [
                _alignment(name, start1, start2, fragment_length, True),
                _alignment(name, start2, start1, -fragment_length, False),
            ]
        )
    reads.sort(key=lambda read: (read.reference_start, read.query_name, not read.is_read1))
    with pysam.AlignmentFile(
        str(path), "wb", header={"HD": {"VN": "1.6", "SO": "coordinate"},
                                 "SQ": [{"SN": "chr1", "LN": 1000}]}
    ) as bam:
        for read in reads:
            bam.write(read)
    return path


def write_bed(tmp_path: Path, intervals: list[tuple[int, int]], gzipped: bool = False) -> Path:
    path = tmp_path / ("peaks.bed.gz" if gzipped else "peaks.bed")
    text = "".join(f"chr1\t{start}\t{end}\n" for start, end in intervals)
    if gzipped:
        with gzip.open(path, "wt") as handle:
            handle.write(text)
    else:
        path.write_text(text)
    return path


def count(tmp_path: Path, pairs: list[tuple[str, int, int]], peaks: list[tuple[int, int]]):
    bam = write_bam(tmp_path, pairs)
    bed = write_bed(tmp_path, peaks)
    return FRIP.count_fragments(bam, FRIP.load_peak_index(bed))


def test_peak_overlapped_only_by_mate_counts_the_fragment(tmp_path):
    # read1 is 100-150, read2 is 200-250. The peak touches only read2, but the
    # biological fragment spans 100-250 and therefore belongs in the numerator.
    assert count(tmp_path, [("pair1", 100, 200)], [(225, 230)]) == (1, 1)


def test_pair_is_counted_once_when_both_mates_overlap(tmp_path):
    # Both 300-350 and 400-450 overlap this peak. BAM-record counting would yield
    # two hits; fragment FRiP must yield one numerator and one denominator count.
    assert count(tmp_path, [("pair1", 300, 400)], [(325, 425)]) == (1, 1)


def test_fragment_with_no_peak_overlap_is_not_in_numerator(tmp_path):
    assert count(tmp_path, [("pair1", 600, 700)], [(10, 20)]) == (1, 0)


def test_multiple_fragments_have_bounded_fragment_frip(tmp_path):
    assert count(
        tmp_path,
        [("mate_only", 100, 200), ("both", 300, 400), ("none", 600, 700)],
        [(225, 230), (325, 425)],
    ) == (3, 2)


def test_gzipped_bed_and_overlapping_peaks_are_supported(tmp_path):
    bam = write_bam(tmp_path, [("pair1", 100, 200)])
    bed = write_bed(tmp_path, [(120, 180), (150, 230)], gzipped=True)
    index = FRIP.load_peak_index(bed)
    assert index["chr1"].starts == (120,)
    assert index["chr1"].ends == (230,)
    assert FRIP.count_fragments(bam, index) == (1, 1)


def test_marked_duplicate_follows_the_final_bam_contract(tmp_path):
    # When remove_duplicates=false, Picard marks but retains duplicates. FRiP
    # must then use the same final-fragment universe as peaks, tracks, and counts.
    bam = write_bam(tmp_path, [("marked", 100, 200)])
    marked = tmp_path / "marked.bam"
    with pysam.AlignmentFile(str(bam), "rb") as source, pysam.AlignmentFile(
        str(marked), "wb", header=source.header
    ) as target:
        for read in source.fetch(until_eof=True):
            read.is_duplicate = True
            target.write(read)
    bed = write_bed(tmp_path, [(120, 180)])
    assert FRIP.count_fragments(marked, FRIP.load_peak_index(bed)) == (1, 1)
