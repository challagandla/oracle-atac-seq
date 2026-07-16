"""Blacklist chromosome-name harmonisation invariants."""

from __future__ import annotations

import subprocess
import sys
import gzip
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "harmonize_blacklist.py"

HUMAN_ENSEMBL = "1\t248956422\n2\t242193529\nMT\t16569\n"
HUMAN_UCSC = "chr1\t248956422\nchr2\t242193529\nchrM\t16569\n"


def run(tmp_path: Path, blacklist: str, chrom_sizes: str, extra: list[str] | None = None):
    bl = tmp_path / "blacklist.bed"
    cs = tmp_path / "chrom.sizes"
    out = tmp_path / "out.bed"
    bl.write_text(blacklist)
    cs.write_text(chrom_sizes)
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), "--blacklist", str(bl),
         "--chrom-sizes", str(cs), "--out", str(out), *(extra or [])],
        capture_output=True, text=True,
    )
    return proc, out


def test_ucsc_blacklist_is_renamed_for_an_ensembl_genome(tmp_path):
    proc, out = run(tmp_path, "chr1\t100\t200\tHigh Signal Region\n", HUMAN_ENSEMBL)
    assert proc.returncode == 0, proc.stderr
    assert out.read_text().splitlines() == ["1\t100\t200\tHigh Signal Region"]


def test_ensembl_blacklist_is_renamed_for_a_ucsc_genome(tmp_path):
    proc, out = run(tmp_path, "1\t100\t200\tHigh Signal Region\n", HUMAN_UCSC)
    assert proc.returncode == 0, proc.stderr
    assert out.read_text().startswith("chr1\t100\t200")


def test_matching_names_pass_through_unchanged(tmp_path):
    proc, out = run(tmp_path, "chr1\t100\t200\tX\n", HUMAN_UCSC)
    assert proc.returncode == 0
    assert out.read_text().splitlines() == ["chr1\t100\t200\tX"]


def test_mitochondrion_aliases_map(tmp_path):
    # chrM and MT are the same sequence under different names.
    proc, out = run(tmp_path, "chrM\t1\t500\tX\n", HUMAN_ENSEMBL)
    assert proc.returncode == 0, proc.stderr
    assert out.read_text().startswith("MT\t1\t500")


def test_unreconcilable_blacklist_is_refused_not_silently_emptied(tmp_path):
    # Fly chromosome names against a human genome. The old behaviour was to
    # intersect nothing and continue; the only safe answer is to stop.
    proc, _ = run(tmp_path, "chr2L\t100\t200\tX\nchr3R\t1\t2\tX\n", HUMAN_ENSEMBL)
    assert proc.returncode != 0
    assert "match a chromosome" in proc.stderr


def test_contigs_absent_from_the_genome_are_dropped_not_fatal(tmp_path):
    # chrEBV ships in the ENCODE hg38 blacklist and is absent from most genomes.
    # One unknown contig among known ones must not fail the run.
    bed = "chr1\t100\t200\tX\nchr2\t100\t200\tX\nchrEBV\t1\t50\tX\n"
    proc, out = run(tmp_path, bed, HUMAN_ENSEMBL)
    assert proc.returncode == 0, proc.stderr
    lines = out.read_text().splitlines()
    assert len(lines) == 2
    assert all(line.split("\t")[0] in {"1", "2"} for line in lines)
    assert "chrEBV" in proc.stderr


def test_mostly_unmappable_blacklist_is_refused(tmp_path):
    # One in three maps: below the default 0.5 floor, so this is a genome/blacklist
    # mismatch rather than a couple of stray contigs.
    bed = "chr1\t1\t2\tX\nchr2L\t1\t2\tX\nchr3R\t1\t2\tX\n"
    proc, _ = run(tmp_path, bed, HUMAN_ENSEMBL)
    assert proc.returncode != 0


def test_output_is_coordinate_sorted(tmp_path):
    bed = "chr2\t500\t600\tX\nchr1\t900\t1000\tX\nchr1\t100\t200\tX\n"
    proc, out = run(tmp_path, bed, HUMAN_ENSEMBL)
    assert proc.returncode == 0
    rows = [line.split("\t") for line in out.read_text().splitlines()]
    keys = [(r[0], int(r[1])) for r in rows]
    assert keys == sorted(keys)


def test_empty_blacklist_is_refused_when_harmonization_was_requested(tmp_path):
    # Missing-blacklist assemblies must disable filtering explicitly; otherwise
    # the run record would claim a filtering step that removed nothing.
    proc, out = run(tmp_path, "", HUMAN_ENSEMBL)
    assert proc.returncode != 0
    assert "empty and contains no intervals" in proc.stderr
    assert not out.exists()


def test_track_and_comment_lines_are_ignored(tmp_path):
    bed = '#comment\ntrack name="bl"\nchr1\t100\t200\tX\n'
    proc, out = run(tmp_path, bed, HUMAN_ENSEMBL)
    assert proc.returncode == 0, proc.stderr
    assert out.read_text().splitlines() == ["1\t100\t200\tX"]


def test_gzipped_blacklist_is_supported(tmp_path):
    bl = tmp_path / "blacklist.bed.gz"
    cs = tmp_path / "chrom.sizes"
    out = tmp_path / "out.bed"
    with gzip.open(bl, "wt") as handle:
        handle.write("chr1\t100\t200\tX\n")
    cs.write_text(HUMAN_ENSEMBL)
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), "--blacklist", str(bl),
         "--chrom-sizes", str(cs), "--out", str(out)],
        capture_output=True, text=True,
    )
    assert proc.returncode == 0, proc.stderr
    assert out.read_text() == "1\t100\t200\tX\n"


@pytest.mark.parametrize("fraction,expect_ok", [("0.2", True), ("0.9", False)])
def test_min_mapped_fraction_is_configurable(tmp_path, fraction, expect_ok):
    bed = "chr1\t1\t2\tX\nchr2L\t1\t2\tX\nchr3R\t1\t2\tX\n"   # 1/3 maps
    proc, _ = run(tmp_path, bed, HUMAN_ENSEMBL, ["--min-mapped-fraction", fraction])
    assert (proc.returncode == 0) is expect_ok
