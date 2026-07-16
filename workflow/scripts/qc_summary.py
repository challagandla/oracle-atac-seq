#!/usr/bin/env python3
"""Summarise sample-level ATAC-seq QC without hiding review decisions."""
from __future__ import annotations

import argparse
import csv

import pysam


def bam_counts(raw_bam, filtered_bam, mitochondrial_contigs):
    mitochondrial_contigs = set(mitochondrial_contigs)
    counts = {"raw_fragments": 0, "mapped_fragments": 0,
              "proper_fragments": 0, "mito_fragments": 0,
              "usable_fragments": 0}
    with pysam.AlignmentFile(raw_bam, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if not read.is_read1 or read.is_secondary or read.is_supplementary:
                continue
            counts["raw_fragments"] += 1
            if read.is_unmapped or read.mate_is_unmapped:
                continue
            counts["mapped_fragments"] += 1
            if not read.is_proper_pair:
                continue
            counts["proper_fragments"] += 1
            ref = bam.get_reference_name(read.reference_id)
            mate = bam.get_reference_name(read.next_reference_id)
            if ref in mitochondrial_contigs or mate in mitochondrial_contigs:
                counts["mito_fragments"] += 1
    with pysam.AlignmentFile(filtered_bam, "rb") as bam:
        counts["usable_fragments"] = sum(
            1 for read in bam.fetch(until_eof=True)
            if read.is_read1 and read.is_proper_pair and not read.is_secondary
            and not read.is_supplementary and not read.is_unmapped
        )
    return counts


def read_data_row(path):
    if not path:
        return {}
    with open(path) as handle:
        rows = [line for line in handle if line.strip() and not line.startswith("#")]
    if len(rows) < 2:
        return {}
    return next(csv.DictReader(rows, delimiter="\t"))


def read_tss(path, sample):
    if not path:
        return None
    with open(path) as handle:
        rows = [line for line in handle if line.strip() and not line.startswith("#")]
    for row in csv.DictReader(rows, delimiter="\t"):
        if row.get("Sample") == sample:
            return float(row["TSS_enrichment"])
    return None


def classify(metrics, genome_build):
    """Return a transparent review label and the rules that triggered it."""
    severity = 0
    notes = []

    def flag(level, note):
        nonlocal severity
        severity = max(severity, level)
        notes.append(note)

    frip = metrics.get("frip")
    if frip is not None:
        if frip < 0.20:
            flag(2, "FRiP<0.20")
        elif frip < 0.30:
            flag(1, "FRiP<0.30")
    tss = metrics.get("tss")
    if tss is not None and genome_build == "human":
        if tss < 5:
            flag(2, "TSS<5")
        elif tss < 7:
            flag(1, "TSS<7")
    mito = metrics.get("mito_fraction")
    if mito is not None:
        if mito > 0.50:
            flag(2, "mitochondrial_fraction>0.50")
        elif mito > 0.30:
            flag(1, "mitochondrial_fraction>0.30")
    usable = metrics.get("usable_fragments", 0)
    if usable < 1_000_000:
        flag(2, "usable_fragments<1M")
    elif usable < 10_000_000:
        flag(1, "usable_fragments<10M")
    peaks = metrics.get("peak_count", 0)
    if peaks < 1_000:
        flag(2, "peak_count<1k")
    elif peaks < 20_000:
        flag(1, "peak_count<20k")
    nrf = metrics.get("nrf")
    if nrf is not None:
        if nrf < 0.50:
            flag(2, "NRF<0.50")
        elif nrf < 0.80:
            flag(1, "NRF<0.80")

    return ("fail" if severity == 2 else "review" if severity == 1 else "pass",
            ";".join(notes) if notes else "none")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--raw-bam", required=True)
    ap.add_argument("--filtered-bam", required=True)
    ap.add_argument("--peaks", required=True)
    ap.add_argument("--frip", required=True)
    ap.add_argument("--complexity", default="")
    ap.add_argument("--tss", default="")
    ap.add_argument("--genome-build", default="custom")
    ap.add_argument(
        "--mitochondrial-contig", action="append", required=True,
        dest="mitochondrial_contigs",
        help="Exact BAM contig name; repeat for every mitochondrial alias/accession",
    )
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    counts = bam_counts(
        args.raw_bam, args.filtered_bam, args.mitochondrial_contigs
    )
    frip_row = read_data_row(args.frip)
    complexity = read_data_row(args.complexity)
    tss = read_tss(args.tss, args.sample)
    if not frip_row or "FRiP" not in frip_row:
        raise ValueError(f"FRiP file has no usable data row: {args.frip}")
    if args.complexity and (not complexity or "NRF" not in complexity):
        raise ValueError(f"complexity file has no usable data row: {args.complexity}")
    if args.tss and tss is None:
        raise ValueError(
            f"TSS table {args.tss} has no usable row for sample {args.sample}"
        )
    proper = counts["proper_fragments"]
    mito_fraction = counts["mito_fragments"] / proper if proper else 0.0
    peak_count = sum(1 for line in open(args.peaks) if line.strip() and not line.startswith("#"))
    metrics = {
        **counts,
        "mito_fraction": mito_fraction,
        "frip": float(frip_row["FRiP"]) if frip_row else None,
        "tss": tss,
        "nrf": float(complexity["NRF"]) if complexity else None,
        "pbc1": float(complexity["PBC1"]) if complexity else None,
        "pbc2": complexity.get("PBC2") if complexity else None,
        "peak_count": peak_count,
    }
    status, notes = classify(metrics, args.genome_build)
    fields = ["sample", "raw_fragments", "mapped_fragments", "proper_fragments",
              "mito_fraction", "usable_fragments", "retained_fraction", "peak_count",
              "FRiP", "TSS_enrichment", "NRF", "PBC1", "PBC2", "QC_status", "QC_notes"]
    retained = counts["usable_fragments"] / proper if proper else 0.0
    row = [args.sample, counts["raw_fragments"], counts["mapped_fragments"], proper,
           f"{mito_fraction:.4f}", counts["usable_fragments"], f"{retained:.4f}", peak_count,
           "NA" if metrics["frip"] is None else f"{metrics['frip']:.4f}",
           "NA" if tss is None else f"{tss:.3f}",
           "NA" if metrics["nrf"] is None else f"{metrics['nrf']:.4f}",
           "NA" if metrics["pbc1"] is None else f"{metrics['pbc1']:.4f}",
           metrics["pbc2"] or "NA", status, notes]
    with open(args.out, "w") as out:
        out.write("\t".join(fields) + "\n")
        out.write("\t".join(map(str, row)) + "\n")


if __name__ == "__main__":
    main()
