#!/usr/bin/env python3
"""Validate narrowPeak records and normalize the UCSC display score."""
from __future__ import annotations

import argparse
import math
import os
import re
import stat
import sys
import tempfile
from pathlib import Path


INTEGER = re.compile(r"-?[0-9]+")


def finite_number(value: str, *, field: str, line_number: int) -> float:
    try:
        number = float(value)
    except ValueError as exc:
        raise ValueError(
            f"line {line_number}: {field} is not numeric: {value!r}"
        ) from exc
    if not math.isfinite(number):
        raise ValueError(f"line {line_number}: {field} must be finite")
    return number


def integer(value: str, *, field: str, line_number: int) -> int:
    if not INTEGER.fullmatch(value):
        raise ValueError(
            f"line {line_number}: {field} must be an integer: {value!r}"
        )
    return int(value)


def normalize_record(line: str, line_number: int) -> str:
    fields = line.rstrip("\n").split("\t")
    if len(fields) != 10:
        raise ValueError(
            f"line {line_number}: narrowPeak requires exactly 10 tab-separated "
            f"fields, found {len(fields)}"
        )

    chrom, start_text, end_text, name, score_text, strand = fields[:6]
    if not chrom or any(char.isspace() for char in chrom):
        raise ValueError(f"line {line_number}: invalid chromosome name {chrom!r}")
    if not name or any(char.isspace() for char in name):
        raise ValueError(f"line {line_number}: invalid peak name {name!r}")

    start = integer(start_text, field="chromStart", line_number=line_number)
    end = integer(end_text, field="chromEnd", line_number=line_number)
    if start < 0 or end <= start:
        raise ValueError(
            f"line {line_number}: require 0 <= chromStart < chromEnd, got "
            f"{start} and {end}"
        )
    if strand not in {".", "+", "-"}:
        raise ValueError(f"line {line_number}: invalid strand {strand!r}")

    score = integer(score_text, field="score", line_number=line_number)
    fields[4] = str(min(1000, max(0, score)))
    signal = finite_number(fields[6], field="signalValue", line_number=line_number)
    if signal < 0:
        raise ValueError(f"line {line_number}: signalValue must be non-negative")
    for index, field in ((7, "pValue"), (8, "qValue")):
        value = finite_number(fields[index], field=field, line_number=line_number)
        if value < 0 and value != -1:
            raise ValueError(
                f"line {line_number}: {field} must be -1 or non-negative"
            )
    summit = integer(fields[9], field="peak", line_number=line_number)
    if summit != -1 and not 0 <= summit < end - start:
        raise ValueError(
            f"line {line_number}: peak offset {summit} lies outside a "
            f"{end - start}-base interval"
        )
    return "\t".join(fields) + "\n"


def normalize_file(input_path: Path, output_path: Path) -> int:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    source_mode = stat.S_IMODE(input_path.stat().st_mode)
    temporary = None
    count = 0
    try:
        with input_path.open(encoding="utf-8") as source, tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=output_path.parent,
            prefix=f".{output_path.name}.",
            suffix=".tmp",
            delete=False,
        ) as destination:
            temporary = Path(destination.name)
            for line_number, line in enumerate(source, start=1):
                if not line.strip():
                    continue
                destination.write(normalize_record(line, line_number))
                count += 1
        if count == 0:
            raise ValueError("narrowPeak contains no records")
        os.chmod(temporary, source_mode)
        os.replace(temporary, output_path)
        temporary = None
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    return count


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    try:
        count = normalize_file(args.input, args.output)
    except (OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
    print(f"normalized {count} narrowPeak records")


if __name__ == "__main__":
    main()
