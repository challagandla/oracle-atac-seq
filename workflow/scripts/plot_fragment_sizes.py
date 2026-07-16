#!/usr/bin/env python3
"""Aggregate ATAC fragment-size distributions into one publication figure.

Reads per-sample `bamPEFragmentSize --outRawFragmentLengths` tables and draws
the canonical ATAC-seq fragment-length distribution: a nucleosome-free peak
(<100 bp) followed by mono-, di- and tri-nucleosomal humps with ~200 bp
periodicity (Buenrostro et al. 2013). Lines are coloured by condition so
within-group consistency is visible at a glance.

Also emits a MultiQC custom-content table of nucleosome-partition fractions
(NFR / mono-nucleosome) — a quantitative companion to the visual check.
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from palette import condition_colours  # noqa: E402

# Nucleosome partitions (bp), by convention.
NFR = (0, 100)          # nucleosome-free
MONO = (180, 247)       # mono-nucleosome
DI = (315, 473)         # di-nucleosome


def sample_name(path):
    return os.path.basename(path).replace(".fragsize.txt", "")


def load_hist(path, max_size):
    sizes = np.zeros(max_size + 1, dtype=float)
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or line.lower().startswith("size"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 2:
                continue
            try:
                s, c = int(float(f[0])), float(f[1])
            except ValueError:
                continue
            if 0 <= s <= max_size:
                sizes[s] += c
    return sizes


def frac(hist, lo, hi):
    total = hist.sum()
    return float(hist[lo:hi + 1].sum() / total) if total > 0 else 0.0


def smooth(y, w):
    """Simple moving average so nucleosome periodicity reads cleanly."""
    if w <= 1:
        return y
    return np.convolve(y, np.ones(w) / w, mode="same")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--frag", nargs="+", required=True,
                    help="per-sample *.fragsize.txt files")
    ap.add_argument("--samples", help="sample sheet TSV (for condition colours)")
    ap.add_argument("--out-plot", required=True)
    ap.add_argument("--out-table", required=True)
    ap.add_argument("--max-size", type=int, default=1000)
    ap.add_argument("--smooth", type=int, default=5,
                    help="moving-average window (bp) for display")
    a = ap.parse_args()

    # Map sample -> condition for colouring, if a sample sheet is given.
    cond_of = {}
    if a.samples and os.path.exists(a.samples):
        with open(a.samples) as fh:
            header = None
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                cols = line.rstrip("\n").split("\t")
                if header is None:
                    header = cols
                    continue
                row = dict(zip(header, cols))
                if "sample" in row and "condition" in row:
                    cond_of[row["sample"]] = row["condition"]

    samples = [sample_name(p) for p in a.frag]
    conds = [cond_of.get(s, s) for s in samples]
    # Sorted condition order gives a colour mapping that is stable across every
    # figure in the pipeline (matches atac_theme.R's sorted assignment).
    uniq_conds = sorted(set(conds))
    colour = condition_colours(conds)

    fig, ax = plt.subplots(figsize=(7.5, 5.0), dpi=200)
    rows = []
    for path, s, c in zip(a.frag, samples, conds):
        hist = load_hist(path, a.max_size)
        total = hist.sum()
        if total == 0:
            continue
        dens = smooth(hist / total, a.smooth)
        x = np.arange(len(dens))
        ax.plot(x, dens, lw=1.1, alpha=0.8, color=colour[c])
        rows.append((s, c, frac(hist, *NFR), frac(hist, *MONO), frac(hist, *DI)))

    # Shade + annotate nucleosome partitions.
    ymax = ax.get_ylim()[1]
    for (lo, hi), lab in ((NFR, "NFR"), (MONO, "mono"), (DI, "di")):
        ax.axvspan(lo, hi, color="#9a988f", alpha=0.06, lw=0)
        ax.text((lo + hi) / 2, ymax * 0.97, lab, ha="center", va="top",
                fontsize=8, color="#52514e")

    ax.set_xlim(0, a.max_size)
    ax.set_xlabel("Fragment length (bp)")
    ax.set_ylabel("Fraction of fragments")
    ax.set_title("ATAC fragment-size distribution", fontsize=13, fontweight="bold",
                 loc="left")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#e6e5df", lw=0.4)
    handles = [plt.Line2D([0], [0], color=colour[c], lw=2) for c in uniq_conds]
    ax.legend(handles, uniq_conds, title="condition", frameon=False,
              fontsize=8, title_fontsize=9, loc="upper right")
    fig.tight_layout()
    fig.savefig(a.out_plot)
    if a.out_plot.endswith(".png"):
        fig.savefig(a.out_plot[:-4] + ".pdf")   # vector copy for publication
        if not a.out_plot.endswith("_mqc.png"):
            # MultiQC auto-embeds *_mqc.png images as their own section.
            fig.savefig(a.out_plot.replace(".png", "_mqc.png"))
    plt.close(fig)

    with open(a.out_table, "w") as out:
        out.write("# id: 'atac_fragsize'\n")
        out.write("# section_name: 'ATAC fragment partitions'\n")
        out.write("# description: 'Fraction of fragments that are nucleosome-free "
                  "(<100 bp), mono- (180-247 bp) and di-nucleosomal (315-473 bp). "
                  "A strong NFR peak and visible nucleosomal periodicity indicate "
                  "a good ATAC library.'\n")
        out.write("# plot_type: 'table'\n")
        out.write("# pconfig:\n")
        out.write("#     id: 'atac_fragsize_table'\n")
        out.write("#     namespace: 'ATAC'\n")
        out.write("Sample\tcondition\tNFR_frac\tmono_frac\tdi_frac\n")
        for s, c, nfr, mono, di in rows:
            out.write(f"{s}\t{c}\t{nfr:.4f}\t{mono:.4f}\t{di:.4f}\n")

    print(f"Wrote {a.out_plot} and {a.out_table} for {len(rows)} samples")


if __name__ == "__main__":
    main()
