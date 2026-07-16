#!/usr/bin/env python3
"""Put a blacklist BED into the reference genome's chromosome naming.

The ENCODE blacklists are distributed with UCSC names (`chr1`). Ensembl genomes
name the same sequence `1`. `bedtools intersect` matches on the name, so a
UCSC blacklist applied to an Ensembl BAM removes exactly nothing -- silently,
with a zero exit status and a log that looks like success.

So: rename, verify, and refuse to continue if the two files cannot be reconciled.
A blacklist that matches no chromosome is a configuration error, not an empty
blacklist.
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

# Sequences that carry a different name rather than a different prefix.
_MITOCHONDRIAL_ALIASES = ("chrM", "MT", "chrMT", "Mito", "M")
_ALIASES = {
    name: tuple(alias for alias in _MITOCHONDRIAL_ALIASES if alias != name)
    for name in _MITOCHONDRIAL_ALIASES
}


def _read_genome_names(chrom_sizes: Path) -> set[str]:
    names = set()
    with open(chrom_sizes) as handle:
        for line in handle:
            if line.strip():
                names.add(line.split("\t", 1)[0])
    if not names:
        sys.exit(f"error: {chrom_sizes} lists no chromosomes")
    return names


def _candidates(name: str) -> list[str]:
    """Every plausible spelling of one chromosome, best first."""
    out = [name]
    if name.startswith("chr"):
        out.append(name[3:])
    else:
        out.append("chr" + name)
    out.extend(_ALIASES.get(name, ()))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--blacklist", required=True, type=Path)
    ap.add_argument("--chrom-sizes", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--min-mapped-fraction", type=float, default=0.5,
                    help="fail if fewer than this fraction of intervals map (default 0.5)")
    args = ap.parse_args()

    args.out.parent.mkdir(parents=True, exist_ok=True)

    # Assemblies without a suitable blacklist must explicitly disable
    # blacklist removal. An enabled but missing/empty source is ambiguous and
    # must never turn a filtering promise into a silent no-op.
    if not args.blacklist.is_file():
        sys.exit(f"error: blacklist does not exist: {args.blacklist}")
    if args.blacklist.stat().st_size == 0:
        sys.exit(f"error: {args.blacklist} is empty and contains no intervals")

    genome = _read_genome_names(args.chrom_sizes)

    rows: list[tuple[str, str]] = []          # (original chrom, raw line)
    opener = gzip.open if args.blacklist.suffix == ".gz" else open
    with opener(args.blacklist, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            rows.append((line.split("\t", 1)[0], line.rstrip("\n")))

    if not rows:
        sys.exit(f"error: {args.blacklist} contains no intervals")

    kept: list[str] = []
    unmapped: dict[str, int] = {}
    renamed = 0
    for chrom, line in rows:
        target = next((c for c in _candidates(chrom) if c in genome), None)
        if target is None:
            unmapped[chrom] = unmapped.get(chrom, 0) + 1
            continue
        if target != chrom:
            renamed += 1
            line = target + line[len(chrom):]
        kept.append(line)

    fraction = len(kept) / len(rows)
    print(f"blacklist intervals      : {len(rows):,}", file=sys.stderr)
    print(f"  mapped to the genome   : {len(kept):,} ({fraction:.1%})", file=sys.stderr)
    print(f"  renamed                : {renamed:,}", file=sys.stderr)
    if unmapped:
        listed = ", ".join(f"{c} ({n})" for c, n in sorted(unmapped.items())[:8])
        print(f"  dropped, not in genome : {sum(unmapped.values()):,}  [{listed}]",
              file=sys.stderr)

    if fraction < args.min_mapped_fraction:
        sys.exit(
            f"error: only {fraction:.1%} of blacklist intervals match a chromosome in\n"
            f"       {args.chrom_sizes}. The blacklist and the genome disagree about\n"
            f"       what a chromosome is called, so filtering would remove nothing.\n"
            f"       blacklist names: {sorted({c for c, _ in rows})[:5]}\n"
            f"       genome names   : {sorted(genome)[:5]}\n"
            f"       Supply a blacklist matching genome.build, or set\n"
            f"       filtering.remove_blacklist: false and say so in your methods."
        )

    # Sorted output so bedtools does not have to re-sort downstream.
    def key(line: str):
        f = line.split("\t")
        return (f[0], int(f[1]))

    args.out.write_text("\n".join(sorted(kept, key=key)) + "\n")
    print(f"wrote {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
