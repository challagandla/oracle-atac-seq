#!/usr/bin/env python3
"""Extract a TSS BED (one position per gene) from an Ensembl/GENCODE GTF.

Usage: gtf_to_tss.py <in.gtf> <out.bed>
Output: BED6 (chrom, tss, tss+1, gene_id, ., strand), TSS = gene start on +
strand / gene end on - strand. Uses 'gene' features.
"""
import sys
import re


def main(gtf_path, out_path):
    seen = set()
    with open(gtf_path) as fh, open(out_path, "w") as out:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            chrom, start, end, strand, attrs = f[0], int(f[3]), int(f[4]), f[6], f[8]
            m = re.search(r'gene_id "([^"]+)"', attrs)
            gid = m.group(1) if m else "NA"
            tss = start - 1 if strand == "+" else end - 1
            if tss < 0:
                tss = 0
            key = (chrom, tss, strand)
            if key in seen:
                continue
            seen.add(key)
            out.write(f"{chrom}\t{tss}\t{tss + 1}\t{gid}\t.\t{strand}\n")
    sys.stderr.write(f"Wrote {len(seen)} TSS records to {out_path}\n")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("usage: gtf_to_tss.py <in.gtf> <out.bed>")
    main(sys.argv[1], sys.argv[2])
