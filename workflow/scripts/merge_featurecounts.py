#!/usr/bin/env python3
"""Validate and merge single-library featureCounts outputs.

The final matrix uses stable sample identifiers instead of filesystem paths.
All feature metadata, sample identities, counts, status rows, and assigned-count
totals must agree before either output is published.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import tempfile
from contextlib import ExitStack
from itertools import zip_longest
from pathlib import Path
from typing import TextIO


FEATURE_COLUMNS = ("Geneid", "Chr", "Start", "End", "Strand", "Length")
SAMPLE_SUFFIX = re.compile(r"\.filtered\.bam$")


class MergeError(ValueError):
    """A per-sample count file violates the merge contract."""


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts", nargs="+", type=Path, required=True)
    parser.add_argument("--summaries", nargs="+", type=Path, required=True)
    parser.add_argument("--samples", nargs="+", required=True)
    parser.add_argument("--out-counts", type=Path, required=True)
    parser.add_argument("--out-summary", type=Path, required=True)
    return parser.parse_args(argv)


def sample_from_column(label: str) -> str:
    """Recover the sample identifier from featureCounts' BAM column."""
    return SAMPLE_SUFFIX.sub("", Path(label).name)


def column_belongs_to_sample(label: str, sample: str) -> bool:
    """Accept a stable sample ID or a native featureCounts BAM-path label."""
    return label == sample or sample_from_column(label) == sample


def count_reader(handle: TextIO, path: Path, sample: str) -> csv.reader:
    lines = (line for line in handle if not line.startswith("#"))
    reader = csv.reader(lines, delimiter="\t")
    try:
        header = next(reader)
    except StopIteration as exc:
        raise MergeError(f"{path}: no tabular header") from exc
    expected = [*FEATURE_COLUMNS, sample]
    observed_sample = (
        sample
        if len(header) == 7 and column_belongs_to_sample(header[6], sample)
        else ""
    )
    observed = [*header[:6], observed_sample]
    if observed != expected or len(header) != 7:
        raise MergeError(
            f"{path}: expected one count column for sample {sample!r}; "
            f"observed header {header!r}"
        )
    return reader


def validate_count_row(
    row: list[str], path: Path, line_number: int
) -> tuple[tuple[str, ...], int]:
    if len(row) != 7:
        raise MergeError(f"{path}:{line_number}: expected 7 columns, found {len(row)}")
    metadata = tuple(row[:6])
    gene_id, chrom, start_text, end_text, strand, length_text = metadata
    if not gene_id or not chrom:
        raise MergeError(f"{path}:{line_number}: empty feature identifier or chromosome")
    if strand not in {"+", "-", "."}:
        raise MergeError(f"{path}:{line_number}: invalid strand {strand!r}")
    try:
        start, end, length = map(int, (start_text, end_text, length_text))
    except ValueError as exc:
        raise MergeError(f"{path}:{line_number}: non-integer feature coordinates") from exc
    if start < 1 or end < start or length != end - start + 1:
        raise MergeError(f"{path}:{line_number}: invalid SAF-derived feature coordinates")
    if not row[6].isdigit():
        raise MergeError(f"{path}:{line_number}: count must be a non-negative integer")
    return metadata, int(row[6])


def parse_summary(path: Path, sample: str) -> tuple[list[str], dict[str, int]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as exc:
            raise MergeError(f"{path}: empty summary") from exc
        if (
            len(header) != 2
            or header[0] != "Status"
            or not column_belongs_to_sample(header[1], sample)
        ):
            raise MergeError(f"{path}: summary does not belong to sample {sample!r}")

        order: list[str] = []
        values: dict[str, int] = {}
        for line_number, row in enumerate(reader, start=2):
            if len(row) != 2 or not row[0] or not row[1].isdigit():
                raise MergeError(f"{path}:{line_number}: invalid summary row")
            if row[0] in values:
                raise MergeError(f"{path}:{line_number}: duplicate status {row[0]!r}")
            order.append(row[0])
            values[row[0]] = int(row[1])

    if "Assigned" not in values:
        raise MergeError(f"{path}: no Assigned row")
    if values["Assigned"] <= 0:
        raise MergeError(f"{path}: zero fragments assigned for sample {sample!r}")
    return order, values


def temporary_output(path: Path) -> tuple[TextIO, Path]:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent, text=True
    )
    return os.fdopen(descriptor, "w", encoding="utf-8", newline=""), Path(temporary)


def merge(args: argparse.Namespace) -> None:
    samples = args.samples
    if not samples or any(not sample for sample in samples) or len(set(samples)) != len(samples):
        raise MergeError("sample identifiers must be non-empty and unique")
    if len(args.counts) != len(samples) or len(args.summaries) != len(samples):
        raise MergeError("counts, summaries, and samples must have the same length")
    if args.out_counts == args.out_summary:
        raise MergeError("count and summary output paths must be different")

    count_handle, count_temporary = temporary_output(args.out_counts)
    summary_handle, summary_temporary = temporary_output(args.out_summary)
    temporary_paths = (count_temporary, summary_temporary)
    try:
        with count_handle, summary_handle, ExitStack() as stack:
            readers = []
            for path, sample in zip(args.counts, samples, strict=True):
                handle = stack.enter_context(path.open(encoding="utf-8", newline=""))
                readers.append(count_reader(handle, path, sample))

            writer = csv.writer(count_handle, delimiter="\t", lineterminator="\n")
            count_handle.write(
                "# Validated merge of independent single-library featureCounts outputs\n"
            )
            writer.writerow([*FEATURE_COLUMNS, *samples])

            seen_ids: set[str] = set()
            count_sums = [0] * len(samples)
            row_count = 0
            for line_number, grouped in enumerate(zip_longest(*readers), start=2):
                if any(row is None for row in grouped):
                    raise MergeError("per-sample count files contain different row counts")
                parsed = [
                    validate_count_row(row, path, line_number)
                    for row, path in zip(grouped, args.counts, strict=True)
                ]
                metadata = parsed[0][0]
                if any(item[0] != metadata for item in parsed[1:]):
                    raise MergeError(
                        f"feature metadata or row order differs at feature {metadata[0]!r}"
                    )
                if metadata[0] in seen_ids:
                    raise MergeError(f"duplicate feature identifier {metadata[0]!r}")
                seen_ids.add(metadata[0])
                values = [item[1] for item in parsed]
                for index, value in enumerate(values):
                    count_sums[index] += value
                writer.writerow([*metadata, *values])
                row_count += 1
            if row_count == 0:
                raise MergeError("per-sample count files contain no features")

            summary_orders: list[list[str]] = []
            summary_values: list[dict[str, int]] = []
            for path, sample in zip(args.summaries, samples, strict=True):
                order, values = parse_summary(path, sample)
                summary_orders.append(order)
                summary_values.append(values)
            status_order = summary_orders[0]
            if any(set(order) != set(status_order) for order in summary_orders[1:]):
                raise MergeError("per-sample summaries contain different status rows")
            for sample, count_sum, values in zip(
                samples, count_sums, summary_values, strict=True
            ):
                if count_sum != values["Assigned"]:
                    raise MergeError(
                        f"sample {sample!r}: matrix sum {count_sum} differs from "
                        f"Assigned summary {values['Assigned']}"
                    )

            summary_writer = csv.writer(
                summary_handle, delimiter="\t", lineterminator="\n"
            )
            summary_writer.writerow(["Status", *samples])
            for status in status_order:
                summary_writer.writerow(
                    [status, *(values[status] for values in summary_values)]
                )

        os.replace(count_temporary, args.out_counts)
        os.replace(summary_temporary, args.out_summary)
    except Exception:
        for path in temporary_paths:
            path.unlink(missing_ok=True)
        raise


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        merge(args)
    except (MergeError, OSError) as exc:
        raise SystemExit(f"error: {exc}") from exc
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
