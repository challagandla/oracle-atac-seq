#!/usr/bin/env python3
"""Compute FRiP (Fraction of Reads in Peaks) for one sample.

FRiP is a core ATAC-seq QC metric (ENCODE recommends > 0.2-0.3 for a good
library). Counts properly-paired primary reads total and those overlapping
peaks, then writes a small TSV that MultiQC can pick up via a custom table.
"""
import argparse
import subprocess
import sys


def count_reads(bam):
    out = subprocess.check_output(
        ["samtools", "view", "-c", "-f", "2", "-F", "1804", bam]
    )
    return int(out.strip())


def count_reads_in_peaks(bam, peaks):
    # bedtools intersect then count unique read pairs overlapping peaks
    p1 = subprocess.Popen(
        ["bedtools", "intersect", "-u", "-abam", bam, "-b", peaks],
        stdout=subprocess.PIPE,
    )
    p2 = subprocess.check_output(
        ["samtools", "view", "-c", "-f", "2", "-F", "1804", "-"], stdin=p1.stdout
    )
    p1.stdout.close()
    return int(p2.strip())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--peaks", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    total = count_reads(a.bam)
    in_peaks = count_reads_in_peaks(a.bam, a.peaks)
    frip = (in_peaks / total) if total else 0.0

    with open(a.out, "w") as fh:
        fh.write("sample\ttotal_reads\treads_in_peaks\tFRiP\n")
        fh.write(f"{a.sample}\t{total}\t{in_peaks}\t{frip:.4f}\n")
    sys.stderr.write(f"{a.sample}: FRiP={frip:.4f} ({in_peaks}/{total})\n")


if __name__ == "__main__":
    main()
