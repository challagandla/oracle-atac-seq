#!/usr/bin/env python3
"""Build a reproducible consensus peak set from per-sample MACS3 narrowPeaks.

Algorithm (Yan et al. 2020; ENCODE consensus approach):
  1. Concatenate all narrowPeak intervals, sort, and merge overlapping peaks.
  2. For each merged interval, count how many distinct samples contributed a
     peak; keep intervals supported by >= min_overlap samples.
  3. Emit a BED (chrom,start,end,name) and a featureCounts SAF.

Pure-Python merge (no bedtools dependency) so it is easy to audit.
"""
import argparse
import sys


def load_peaks(path):
    iv = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.split("\t")
            iv.append((f[0], int(f[1]), int(f[2])))
    return iv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--peaks", nargs="+", required=True)
    ap.add_argument("--chrom", required=True, help="chrom.sizes (for bounds)")
    ap.add_argument("--min-overlap", type=int, default=2)
    ap.add_argument("--bed", required=True)
    ap.add_argument("--saf", required=True)
    a = ap.parse_args()

    chrom_order = {}
    with open(a.chrom) as fh:
        for i, line in enumerate(fh):
            chrom_order[line.split("\t")[0]] = i

    # tag each interval with its source sample index
    tagged = []
    for sidx, p in enumerate(a.peaks):
        for chrom, start, end in load_peaks(p):
            tagged.append((chrom, start, end, sidx))

    # sort by chrom (by chrom.sizes order, then lexical), then start
    def chrom_key(c):
        return chrom_order.get(c, 10_000 + hash(c) % 1000)

    tagged.sort(key=lambda x: (chrom_key(x[0]), x[1], x[2]))

    consensus = []
    cur = None
    for chrom, start, end, sidx in tagged:
        if cur and chrom == cur[0] and start <= cur[2]:
            cur[2] = max(cur[2], end)
            cur[3].add(sidx)
        else:
            if cur:
                consensus.append(cur)
            cur = [chrom, start, end, {sidx}]
    if cur:
        consensus.append(cur)

    kept = [c for c in consensus if len(c[3]) >= a.min_overlap]

    with open(a.bed, "w") as bed, open(a.saf, "w") as saf:
        saf.write("GeneID\tChr\tStart\tEnd\tStrand\n")
        for i, (chrom, start, end, samples) in enumerate(kept, 1):
            name = f"peak_{i}"
            bed.write(f"{chrom}\t{start}\t{end}\t{name}\n")
            # SAF is 1-based inclusive
            saf.write(f"{name}\t{chrom}\t{start + 1}\t{end}\t+\n")

    sys.stderr.write(
        f"Consensus: {len(kept)} peaks kept "
        f"(>= {a.min_overlap} samples) of {len(consensus)} merged.\n"
    )


if __name__ == "__main__":
    main()
