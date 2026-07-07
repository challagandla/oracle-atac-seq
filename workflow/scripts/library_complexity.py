#!/usr/bin/env python3
"""Library-complexity metrics (ENCODE ATAC): NRF, PBC1, PBC2.

Computed from the pre-deduplication alignment, using the genomic location of
each properly-paired fragment (read1 5' + mate 5' + strand):

  M_distinct = distinct fragment locations
  M1         = locations covered by exactly one fragment
  M2         = locations covered by exactly two fragments
  total      = total fragments considered

  NRF  = M_distinct / total          (non-redundant fraction; ENCODE >0.9 ideal)
  PBC1 = M1 / M_distinct             (PCR bottleneck 1; >0.9 none/mild)
  PBC2 = M1 / M2                     (PCR bottleneck 2; >3 none/mild)

Memory is bounded by processing one reference sequence at a time.
"""
import argparse

import pysam


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-mapq", type=int, default=1)
    a = ap.parse_args()

    bam = pysam.AlignmentFile(a.bam, "rb")
    total = m_distinct = m1 = m2 = 0

    for tid in range(bam.nreferences):
        loc = {}
        for r in bam.fetch(bam.references[tid]):
            # count each fragment once, from read1 of a proper primary pair
            if not r.is_proper_pair or not r.is_read1:
                continue
            if r.is_secondary or r.is_supplementary or r.is_unmapped or r.mate_is_unmapped:
                continue
            if r.mapping_quality < a.min_mapq:
                continue
            key = (r.reference_start, r.next_reference_start, r.is_reverse)
            loc[key] = loc.get(key, 0) + 1
        for n in loc.values():
            total += n
            m_distinct += 1
            if n == 1:
                m1 += 1
            elif n == 2:
                m2 += 1
    bam.close()

    nrf = m_distinct / total if total else 0.0
    pbc1 = m1 / m_distinct if m_distinct else 0.0
    pbc2 = m1 / m2 if m2 else float("nan")

    with open(a.out, "w") as fh:
        fh.write("sample\ttotal_frags\tdistinct\tNRF\tPBC1\tPBC2\n")
        pbc2_str = f"{pbc2:.4f}" if pbc2 == pbc2 else "NA"
        fh.write(f"{a.sample}\t{total}\t{m_distinct}\t{nrf:.4f}\t{pbc1:.4f}\t{pbc2_str}\n")
    print(f"{a.sample}: NRF={nrf:.4f} PBC1={pbc1:.4f} PBC2={pbc2:.4f} "
          f"(total={total}, distinct={m_distinct})")


if __name__ == "__main__":
    main()
