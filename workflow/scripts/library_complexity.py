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


def complexity_metrics(total, distinct, singletons, doubletons):
    """Return NRF/PBC values while keeping their integer components auditable."""
    nrf = distinct / total if total else 0.0
    pbc1 = singletons / distinct if distinct else 0.0
    pbc2 = singletons / doubletons if doubletons else float("nan")
    return nrf, pbc1, pbc2


def fragment_key(read):
    """Return true paired-fragment bounds and orientation, or ``None``.

    SAM ``reference_start`` is the left edge of an alignment, not the 5' end
    of a reverse-strand mate. Template length instead spans the complete
    paired fragment and is stable when a reverse read is adapter-trimmed.
    """
    if (
        not read.is_paired
        or not read.is_proper_pair
        or not read.is_read1
        or read.is_unmapped
        or read.mate_is_unmapped
        or read.is_secondary
        or read.is_supplementary
        or read.reference_id != read.next_reference_id
        or read.template_length == 0
    ):
        return None
    start = min(read.reference_start, read.next_reference_start)
    end = start + abs(read.template_length)
    if start < 0 or end <= start:
        return None
    return start, end, read.is_reverse


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
            if r.mapping_quality < a.min_mapq:
                continue
            key = fragment_key(r)
            if key is None:
                continue
            loc[key] = loc.get(key, 0) + 1
        for n in loc.values():
            total += n
            m_distinct += 1
            if n == 1:
                m1 += 1
            elif n == 2:
                m2 += 1
    bam.close()

    nrf, pbc1, pbc2 = complexity_metrics(total, m_distinct, m1, m2)

    with open(a.out, "w") as fh:
        fh.write(
            "sample\ttotal_frags\tdistinct\tM1_singletons\tM2_doubletons\t"
            "NRF\tPBC1\tPBC2\n"
        )
        pbc2_str = f"{pbc2:.4f}" if pbc2 == pbc2 else "NA"
        fh.write(
            f"{a.sample}\t{total}\t{m_distinct}\t{m1}\t{m2}\t"
            f"{nrf:.4f}\t{pbc1:.4f}\t{pbc2_str}\n"
        )
    print(f"{a.sample}: NRF={nrf:.4f} PBC1={pbc1:.4f} PBC2={pbc2:.4f} "
          f"(total={total}, distinct={m_distinct}, M1={m1}, M2={m2})")


if __name__ == "__main__":
    main()
