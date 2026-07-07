#!/usr/bin/env python3
"""Compute a numeric TSS-enrichment score from a deepTools computeMatrix file.

TSS enrichment is the headline ATAC signal-to-noise metric (ENCODE): the mean
signal in a window around the TSS divided by the background signal in the distal
flanks. A score >6-7 is a good human ATAC library; <5 is marginal.

Handles both single- and multi-sample matrices (reference-point, TSS-centred):
the JSON header's `sample_labels` / `sample_boundaries` split the value block
per sample. Emits a MultiQC custom-content table.
"""
import argparse
import gzip
import json

import numpy as np


def parse_matrix(path):
    with gzip.open(path, "rt") as fh:
        header = fh.readline()
        hdr = json.loads(header[1:] if header.startswith("@") else header)
        labels = hdr["sample_labels"]
        bounds = hdr["sample_boundaries"]     # value-block boundaries per sample
        # Column means (nan-aware) across all regions, streamed to bound memory.
        ncols = bounds[-1]
        ssum = np.zeros(ncols)
        scount = np.zeros(ncols)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            vals = np.array([np.nan if v in ("nan", "NA", "") else float(v)
                             for v in f[6:6 + ncols]])
            mask = ~np.isnan(vals)
            ssum[mask] += vals[mask]
            scount[mask] += 1
    profile = np.where(scount > 0, ssum / np.maximum(scount, 1), np.nan)
    return labels, bounds, profile


def enrichment(profile, flank_bins=10):
    """max(profile) / mean(flank background), background from distal edges."""
    prof = np.nan_to_num(profile, nan=0.0)
    if prof.size == 0:
        return 0.0
    bg = np.concatenate([prof[:flank_bins], prof[-flank_bins:]])
    bg = bg[bg > 0]
    background = bg.mean() if bg.size else np.nan
    if not background or np.isnan(background):
        return 0.0
    return float(np.nanmax(prof) / background)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", required=True, help="computeMatrix .gz")
    ap.add_argument("--out", required=True, help="MultiQC _mqc.tsv table")
    ap.add_argument("--flank-bins", type=int, default=10)
    a = ap.parse_args()

    labels, bounds, profile = parse_matrix(a.matrix)
    rows = []
    for i, lab in enumerate(labels):
        seg = profile[bounds[i]:bounds[i + 1]]
        rows.append((lab, enrichment(seg, a.flank_bins)))

    with open(a.out, "w") as out:
        out.write("# id: 'tss_enrichment'\n")
        out.write("# section_name: 'TSS enrichment score'\n")
        out.write("# description: 'Mean ATAC signal at the TSS divided by distal "
                  "background (deepTools). ENCODE: >6-7 is a good human library, "
                  "<5 is marginal.'\n")
        out.write("# plot_type: 'bargraph'\n")
        out.write("# pconfig:\n")
        out.write("#     id: 'tss_enrichment_plot'\n")
        out.write("#     namespace: 'ATAC'\n")
        out.write("#     ylab: 'TSS enrichment'\n")
        out.write("Sample\tTSS_enrichment\n")
        for lab, score in rows:
            out.write(f"{lab}\t{score:.3f}\n")
    for lab, score in rows:
        print(f"{lab}: TSS enrichment = {score:.3f}")


if __name__ == "__main__":
    main()
