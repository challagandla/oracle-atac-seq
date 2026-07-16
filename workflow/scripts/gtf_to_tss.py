#!/usr/bin/env python3
"""Extract a genome-validated, one-base TSS BED from a gene GTF.

GTF coordinates are 1-based inclusive; BED coordinates are 0-based half-open.
Only ``gene`` features are used. Chromosome names are reconciled conservatively
against the alignment genome (direct match, ``chr`` prefix, mitochondrial
aliases), and severe partial reference mismatches fail before TSS QC can draw a
biased profile from the surviving subset.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


MITOCHONDRIAL_ALIASES = ("chrM", "MT", "chrMT", "Mito", "M")
ALIASES = {
    name: tuple(alias for alias in MITOCHONDRIAL_ALIASES if alias != name)
    for name in MITOCHONDRIAL_ALIASES
}


def read_chrom_sizes(path: str | Path) -> tuple[dict[str, int], dict[str, int]]:
    sizes = {}
    order = {}
    with open(path) as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError(f"{path}:{line_no}: expected chromosome and size")
            try:
                size = int(fields[1])
            except ValueError as exc:
                raise ValueError(f"{path}:{line_no}: invalid chromosome size") from exc
            if size <= 0 or fields[0] in sizes:
                raise ValueError(
                    f"{path}:{line_no}: chromosome size must be positive and unique"
                )
            order[fields[0]] = len(order)
            sizes[fields[0]] = size
    if not sizes:
        raise ValueError(f"{path}: no chromosome sizes found")
    return sizes, order


def contig_candidates(name: str) -> list[str]:
    candidates = [name]
    candidates.append(name[3:] if name.startswith("chr") else "chr" + name)
    candidates.extend(ALIASES.get(name, ()))
    return list(dict.fromkeys(candidates))


def extract_tss(
    gtf_path: str | Path,
    chrom_sizes_path: str | Path,
    *,
    min_retained_fraction: float = 0.95,
    warn_retained_fraction: float = 0.99,
) -> tuple[list[tuple[str, int, int, str, str, str]], dict[str, float | int]]:
    """Return valid BED6 records and mapping statistics."""

    sizes, order = read_chrom_sizes(chrom_sizes_path)
    records = []
    total = mapped = renamed = dropped_contig = dropped_bounds = 0

    with open(gtf_path) as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            total += 1
            chrom, strand, attrs = fields[0], fields[6], fields[8]
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError as exc:
                raise ValueError(f"{gtf_path}:{line_no}: non-integer gene coordinates") from exc
            if start < 1 or end < start or strand not in {"+", "-"}:
                raise ValueError(
                    f"{gtf_path}:{line_no}: invalid gene interval or strand"
                )

            target = next((name for name in contig_candidates(chrom) if name in sizes), None)
            if target is None:
                dropped_contig += 1
                continue
            tss = start - 1 if strand == "+" else end - 1
            if tss < 0 or tss >= sizes[target]:
                dropped_bounds += 1
                continue

            gene_match = re.search(r'gene_id "([^"]+)"', attrs)
            gene_id = gene_match.group(1) if gene_match else f"gene_line_{line_no}"
            records.append((target, tss, tss + 1, gene_id, ".", strand))
            mapped += 1
            renamed += int(target != chrom)

    if total == 0:
        raise ValueError(f"{gtf_path}: no gene features were found")
    retained_fraction = mapped / total
    if mapped == 0 or retained_fraction < min_retained_fraction:
        raise ValueError(
            f"only {mapped}/{total} gene TSSs ({retained_fraction:.1%}) map inside "
            f"the alignment genome; refusing biased TSS QC"
        )

    # Multiple gene records can legitimately share one TSS. DeepTools needs the
    # genomic position once; keep the first stable identifier at each strand.
    unique = {}
    for record in records:
        unique.setdefault((record[0], record[1], record[5]), record)
    records = sorted(
        unique.values(), key=lambda row: (order[row[0]], row[1], row[5], row[3])
    )
    if not records:
        raise ValueError("no unique, genome-valid TSS records remain")

    stats = {
        "gene_features": total,
        "mapped": mapped,
        "renamed": renamed,
        "dropped_contig": dropped_contig,
        "dropped_bounds": dropped_bounds,
        "unique_tss": len(records),
        "retained_fraction": retained_fraction,
        "warn": int(retained_fraction < warn_retained_fraction),
    }
    return records, stats


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gtf")
    parser.add_argument("out")
    parser.add_argument("--chrom-sizes", required=True)
    args = parser.parse_args(argv)

    records, stats = extract_tss(args.gtf, args.chrom_sizes)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as handle:
        for record in records:
            handle.write("\t".join(map(str, record)) + "\n")

    print(f"GTF gene features       : {stats['gene_features']:,}", file=sys.stderr)
    print(
        f"  mapped inside genome  : {stats['mapped']:,} "
        f"({stats['retained_fraction']:.1%})",
        file=sys.stderr,
    )
    print(f"  chromosome renamed    : {stats['renamed']:,}", file=sys.stderr)
    print(f"  absent contig dropped : {stats['dropped_contig']:,}", file=sys.stderr)
    print(f"  out-of-bounds dropped : {stats['dropped_bounds']:,}", file=sys.stderr)
    if stats["warn"]:
        print("warning: fewer than 99% of gene TSSs matched the genome", file=sys.stderr)
    print(f"Wrote {stats['unique_tss']:,} unique TSS records to {out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ValueError as exc:
        raise SystemExit(f"error: {exc}") from exc
