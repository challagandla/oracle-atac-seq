#!/usr/bin/env python3
"""Build a replicate-supported, fixed-width ATAC consensus peak set.

Each input narrowPeak file represents one biological replicate and has one
aligned ``--conditions`` value.  Peaks are reduced to fixed-width windows around
their narrowPeak summits.  A candidate is eligible only when its window directly
overlaps peaks from at least ``--min-overlap`` distinct replicates of the same
condition.  Deterministic non-maximum suppression then selects non-overlapping
representatives without ever taking a transitive union of neighbouring peaks.
"""

from __future__ import annotations

import argparse
import bisect
import math
import sys
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Peak:
    chrom: str
    summit: int
    start: int
    end: int
    signal: float
    replicate: int
    condition: str
    source: str
    line: int


@dataclass(frozen=True)
class Candidate:
    peak: Peak
    support: int
    support_signal: float


def read_chrom_sizes(path: str | Path) -> tuple[dict[str, int], dict[str, int]]:
    sizes: dict[str, int] = {}
    order: dict[str, int] = {}
    with open(path) as handle:
        for lineno, line in enumerate(handle, 1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError(f"{path}:{lineno}: expected chromosome and size")
            try:
                size = int(fields[1])
            except ValueError as exc:
                raise ValueError(f"{path}:{lineno}: chromosome size is not an integer") from exc
            if size <= 0:
                raise ValueError(f"{path}:{lineno}: chromosome size must be positive")
            if fields[0] in sizes:
                raise ValueError(f"{path}:{lineno}: duplicate chromosome {fields[0]!r}")
            order[fields[0]] = len(order)
            sizes[fields[0]] = size
    if not sizes:
        raise ValueError(f"{path}: no chromosome sizes found")
    return sizes, order


def fixed_window(summit: int, width: int, chrom_size: int) -> tuple[int, int]:
    """Center a fixed-width half-open window, shifting it inside chromosome bounds."""

    if width <= 0:
        raise ValueError("--width must be a positive integer")
    width = min(width, chrom_size)
    start = summit - width // 2
    start = max(0, min(start, chrom_size - width))
    return start, start + width


def _number(fields: list[str], index: int, fallback: float) -> float:
    if len(fields) <= index:
        return fallback
    try:
        value = float(fields[index])
    except ValueError:
        return fallback
    return value if math.isfinite(value) and value >= 0 else fallback


def load_peaks(
    path: str | Path,
    replicate: int,
    condition: str,
    width: int,
    chrom_sizes: dict[str, int],
) -> list[Peak]:
    """Parse narrowPeak summit/signal fields and create bounded fixed windows."""

    source = str(path)
    loaded: list[Peak] = []
    with open(path) as handle:
        for lineno, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"{path}:{lineno}: expected at least 3 peak columns")
            chrom = fields[0]
            if chrom not in chrom_sizes:
                raise ValueError(f"{path}:{lineno}: chromosome {chrom!r} is absent from chrom sizes")
            try:
                peak_start, peak_end = int(fields[1]), int(fields[2])
            except ValueError as exc:
                raise ValueError(f"{path}:{lineno}: non-integer peak coordinates") from exc
            if peak_start < 0 or peak_end <= peak_start:
                raise ValueError(
                    f"{path}:{lineno}: invalid peak interval {chrom}:{peak_start}-{peak_end}"
                )

            # narrowPeak column 10 is the zero-based summit offset. MACS may emit
            # -1 when no summit is available; use the interval midpoint then.
            summit = (peak_start + peak_end) // 2
            if len(fields) > 9:
                try:
                    offset = int(fields[9])
                except ValueError:
                    offset = -1
                if 0 <= offset < peak_end - peak_start:
                    summit = peak_start + offset
            summit = max(0, min(summit, chrom_sizes[chrom] - 1))

            # Prefer narrowPeak signalValue (column 7); score (column 5) is a
            # deterministic fallback for BED-like fixtures or missing values.
            score = _number(fields, 4, 0.0)
            signal = _number(fields, 6, score)
            start, end = fixed_window(summit, width, chrom_sizes[chrom])
            loaded.append(
                Peak(chrom, summit, start, end, signal, replicate, condition,
                     source, lineno)
            )
    return loaded


def _support(candidate: Peak, neighbours: list[Peak]) -> tuple[int, float]:
    """Count distinct overlapping replicate files and sum their best signals."""

    best_by_replicate: dict[int, float] = {}
    for peak in neighbours:
        if peak.start < candidate.end and candidate.start < peak.end:
            best_by_replicate[peak.replicate] = max(
                best_by_replicate.get(peak.replicate, float("-inf")), peak.signal
            )
    return len(best_by_replicate), sum(best_by_replicate.values())


def eligible_candidates(peaks: list[Peak], min_overlap: int, width: int) -> list[Candidate]:
    if min_overlap <= 0:
        raise ValueError("--min-overlap must be a positive integer")

    by_group: dict[tuple[str, str], list[Peak]] = {}
    for peak in peaks:
        by_group.setdefault((peak.chrom, peak.condition), []).append(peak)

    sorted_groups: dict[tuple[str, str], tuple[list[int], list[Peak]]] = {}
    for key, group in by_group.items():
        ordered = sorted(
            group,
            key=lambda p: (p.summit, p.start, p.end, -p.signal, p.source, p.line),
        )
        sorted_groups[key] = ([peak.summit for peak in ordered], ordered)

    eligible: list[Candidate] = []
    for peak in peaks:
        summits, group = sorted_groups[(peak.chrom, peak.condition)]
        # Bounds-shifted windows near chromosome ends can overlap even when their
        # summits differ by more than one width. A two-width prefilter is a safe
        # upper bound; _support performs the exact half-open interval test.
        left = bisect.bisect_right(summits, peak.summit - 2 * width)
        right = bisect.bisect_left(summits, peak.summit + 2 * width)
        support, support_signal = _support(peak, group[left:right])
        if support >= min_overlap:
            eligible.append(Candidate(peak, support, support_signal))
    return eligible


def select_non_overlapping(
    candidates: list[Candidate], chrom_order: dict[str, int]
) -> list[Peak]:
    """Deterministic non-maximum suppression over fixed-width candidates."""

    ranked = sorted(
        candidates,
        key=lambda c: (
            -c.support,
            -c.support_signal,
            -c.peak.signal,
            chrom_order[c.peak.chrom],
            c.peak.start,
            c.peak.end,
            c.peak.summit,
            c.peak.condition,
            c.peak.source,
            c.peak.line,
        ),
    )

    selected: dict[str, list[Peak]] = {}
    starts: dict[str, list[int]] = {}
    for candidate in ranked:
        peak = candidate.peak
        chrom_peaks = selected.setdefault(peak.chrom, [])
        chrom_starts = starts.setdefault(peak.chrom, [])
        idx = bisect.bisect_left(chrom_starts, peak.start)
        overlaps_left = idx > 0 and chrom_peaks[idx - 1].end > peak.start
        overlaps_right = idx < len(chrom_peaks) and chrom_peaks[idx].start < peak.end
        if overlaps_left or overlaps_right:
            continue
        chrom_starts.insert(idx, peak.start)
        chrom_peaks.insert(idx, peak)

    out = [peak for chrom_peaks in selected.values() for peak in chrom_peaks]
    return sorted(out, key=lambda p: (chrom_order[p.chrom], p.start, p.end, p.summit))


def read_blacklist(path: str | Path) -> dict[str, tuple[list[int], list[tuple[int, int]]]]:
    """Read and merge half-open blacklist intervals for fast overlap checks."""

    by_chrom: dict[str, list[tuple[int, int]]] = {}
    with open(path) as handle:
        for lineno, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"{path}:{lineno}: expected at least 3 BED columns")
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as exc:
                raise ValueError(f"{path}:{lineno}: non-integer blacklist coordinates") from exc
            if start < 0 or end <= start:
                raise ValueError(f"{path}:{lineno}: invalid blacklist interval")
            by_chrom.setdefault(fields[0], []).append((start, end))

    indexed: dict[str, tuple[list[int], list[tuple[int, int]]]] = {}
    for chrom, intervals in by_chrom.items():
        merged: list[list[int]] = []
        for start, end in sorted(intervals):
            if merged and start <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], end)
            else:
                merged.append([start, end])
        frozen = [(start, end) for start, end in merged]
        indexed[chrom] = ([start for start, _ in frozen], frozen)
    return indexed


def remove_blacklisted(
    peaks: list[Peak],
    blacklist: dict[str, tuple[list[int], list[tuple[int, int]]]],
) -> list[Peak]:
    """Drop whole replicate windows that overlap any blacklist interval."""

    retained = []
    for peak in peaks:
        indexed = blacklist.get(peak.chrom)
        if indexed is None:
            retained.append(peak)
            continue
        starts, intervals = indexed
        pos = bisect.bisect_left(starts, peak.end)
        overlaps = pos > 0 and intervals[pos - 1][1] > peak.start
        if not overlaps:
            retained.append(peak)
    return retained


def write_outputs(peaks: list[Peak], bed_path: str | Path, saf_path: str | Path) -> None:
    Path(bed_path).parent.mkdir(parents=True, exist_ok=True)
    Path(saf_path).parent.mkdir(parents=True, exist_ok=True)
    with open(bed_path, "w") as bed, open(saf_path, "w") as saf:
        saf.write("GeneID\tChr\tStart\tEnd\tStrand\n")
        for index, peak in enumerate(peaks, 1):
            name = f"peak_{index}"
            bed.write(f"{peak.chrom}\t{peak.start}\t{peak.end}\t{name}\n")
            saf.write(f"{name}\t{peak.chrom}\t{peak.start + 1}\t{peak.end}\t+\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--peaks", nargs="+", required=True)
    parser.add_argument(
        "--conditions", nargs="+", required=True,
        help="one biological condition per --peaks file, in the same order",
    )
    parser.add_argument("--chrom", required=True, help="chrom.sizes")
    parser.add_argument(
        "--min-overlap",
        "--min-replicates",
        dest="min_overlap",
        type=int,
        default=2,
        help="minimum distinct replicate files supporting a locus within one condition",
    )
    parser.add_argument("--width", type=int, default=500)
    parser.add_argument(
        "--blacklist", default="",
        help="optional harmonized BED; overlapping candidate windows are removed",
    )
    parser.add_argument("--bed", required=True)
    parser.add_argument("--saf", required=True)
    args = parser.parse_args()

    if len(args.conditions) != len(args.peaks):
        parser.error(
            f"--conditions has {len(args.conditions)} values but --peaks has "
            f"{len(args.peaks)} files"
        )
    if args.min_overlap <= 0:
        parser.error("--min-overlap must be positive")
    if args.width <= 0:
        parser.error("--width must be positive")

    chrom_sizes, chrom_order = read_chrom_sizes(args.chrom)
    all_peaks: list[Peak] = []
    for replicate, (path, condition) in enumerate(zip(args.peaks, args.conditions)):
        all_peaks.extend(load_peaks(path, replicate, condition, args.width, chrom_sizes))

    before_blacklist = len(all_peaks)
    if args.blacklist:
        all_peaks = remove_blacklisted(all_peaks, read_blacklist(args.blacklist))
    eligible = eligible_candidates(all_peaks, args.min_overlap, args.width)
    consensus = select_non_overlapping(eligible, chrom_order)
    write_outputs(consensus, args.bed, args.saf)
    print(
        f"Consensus: {len(consensus)} fixed-width peaks from {len(all_peaks)} inputs; "
        f"each supported by >= {args.min_overlap} replicate(s) within a condition; "
        f"{before_blacklist - len(all_peaks)} blacklisted replicate window(s) removed.",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
