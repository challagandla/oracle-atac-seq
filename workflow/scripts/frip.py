#!/usr/bin/env python3
"""Compute paired-end fragment FRiP for one ATAC-seq library.

FRiP is the fraction of usable fragments that overlap at least one peak.  A
paired-end fragment is counted once, from its read-1 alignment; counting BAM
records instead counts the two mates separately and can miss a peak overlapped
only by read 2.  This implementation reconstructs the complete fragment from
the template length and tests that interval against a merged BED index.
"""

from __future__ import annotations

import argparse
import bisect
import gzip
import sys
from dataclasses import dataclass
from pathlib import Path

import pysam


@dataclass(frozen=True)
class IntervalIndex:
    """Merged, half-open BED intervals for one chromosome."""

    starts: tuple[int, ...]
    ends: tuple[int, ...]

    def overlaps(self, start: int, end: int) -> bool:
        if start >= end or not self.starts:
            return False
        # The last interval starting before ``end`` is the only possible hit:
        # intervals are sorted and merged, so every earlier interval ends no
        # later than it does.
        idx = bisect.bisect_left(self.starts, end) - 1
        return idx >= 0 and self.ends[idx] > start


def _open_text(path: str | Path):
    path = Path(path)
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def load_peak_index(path: str | Path) -> dict[str, IntervalIndex]:
    """Load BED/narrowPeak intervals, merging overlaps per chromosome."""

    raw: dict[str, list[tuple[int, int]]] = {}
    with _open_text(path) as handle:
        for lineno, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"{path}:{lineno}: expected at least 3 BED columns")
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as exc:
                raise ValueError(f"{path}:{lineno}: non-integer BED coordinates") from exc
            if start < 0 or end <= start:
                raise ValueError(
                    f"{path}:{lineno}: invalid half-open interval {fields[0]}:{start}-{end}"
                )
            raw.setdefault(fields[0], []).append((start, end))

    indexed: dict[str, IntervalIndex] = {}
    for chrom, intervals in raw.items():
        merged: list[list[int]] = []
        for start, end in sorted(intervals):
            if merged and start <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], end)
            else:
                merged.append([start, end])
        indexed[chrom] = IntervalIndex(
            tuple(interval[0] for interval in merged),
            tuple(interval[1] for interval in merged),
        )
    return indexed


def usable_read1(read: pysam.AlignedSegment) -> bool:
    """Return whether ``read`` represents one usable paired-end fragment."""

    return (
        read.is_paired
        and read.is_read1
        and read.is_proper_pair
        and not read.is_unmapped
        and not read.mate_is_unmapped
        and not read.is_secondary
        and not read.is_supplementary
        and not read.is_qcfail
        and read.reference_id == read.next_reference_id
        and read.template_length != 0
    )


def fragment_interval(read: pysam.AlignedSegment) -> tuple[int, int]:
    """Return the zero-based, half-open interval spanned by a paired fragment."""

    start = min(read.reference_start, read.next_reference_start)
    end = start + abs(read.template_length)
    return start, end


def count_fragments(
    bam_path: str | Path, peak_index: dict[str, IntervalIndex]
) -> tuple[int, int]:
    """Return ``(usable_fragments, fragments_overlapping_peaks)``."""

    total = in_peaks = 0
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if not usable_read1(read):
                continue
            start, end = fragment_interval(read)
            if start < 0 or end <= start:
                continue
            total += 1
            chrom = bam.get_reference_name(read.reference_id)
            index = peak_index.get(chrom)
            if index is not None and index.overlaps(start, end):
                in_peaks += 1
    return total, in_peaks


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--peaks", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    peak_index = load_peak_index(args.peaks)
    total, in_peaks = count_fragments(args.bam, peak_index)
    if total == 0:
        raise SystemExit(
            f"{args.bam} has no usable proper primary paired-end fragments; FRiP is undefined."
        )
    frip = in_peaks / total

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as handle:
        handle.write("sample\ttotal_fragments\tfragments_in_peaks\tFRiP\n")
        handle.write(f"{args.sample}\t{total}\t{in_peaks}\t{frip:.4f}\n")
    print(
        f"{args.sample}: FRiP={frip:.4f} ({in_peaks}/{total} fragments)",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
