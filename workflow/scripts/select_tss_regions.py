#!/usr/bin/env python3
"""Create one reproducible genome-wide TSS subset for every QC consumer."""
from __future__ import annotations

import argparse
import os
import stat
import tempfile
from pathlib import Path


def read_regions(path: Path) -> list[str]:
    regions = []
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(
                    f"line {line_number}: TSS BED requires at least 3 fields"
                )
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as exc:
                raise ValueError(
                    f"line {line_number}: TSS coordinates must be integers"
                ) from exc
            if not fields[0] or start < 0 or end <= start:
                raise ValueError(
                    f"line {line_number}: require chromosome and "
                    "0 <= start < end"
                )
            regions.append(line if line.endswith("\n") else line + "\n")
    if not regions:
        raise ValueError("TSS BED contains no valid regions")
    return regions


def evenly_spaced_indices(total: int, requested: int) -> list[int]:
    """Choose the midpoint of equal-width partitions across the full BED."""
    if requested <= 0 or requested >= total:
        return list(range(total))
    return [((2 * index + 1) * total) // (2 * requested)
            for index in range(requested)]


def select_regions(input_path: Path, output_path: Path, max_regions: int) -> tuple[int, int]:
    if max_regions < 0:
        raise ValueError("max_regions must be zero or a positive integer")
    regions = read_regions(input_path)
    indices = evenly_spaced_indices(len(regions), max_regions)
    source_mode = stat.S_IMODE(input_path.stat().st_mode)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=output_path.parent,
            prefix=f".{output_path.name}.",
            suffix=".tmp",
            delete=False,
        ) as destination:
            temporary = Path(destination.name)
            for index in indices:
                destination.write(regions[index])
        os.chmod(temporary, source_mode)
        os.replace(temporary, output_path)
        temporary = None
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    return len(indices), len(regions)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--max-regions", required=True, type=int)
    args = parser.parse_args()
    try:
        selected, total = select_regions(
            args.input, args.output, args.max_regions
        )
    except (OSError, ValueError) as exc:
        parser.exit(1, f"error: {exc}\n")
    print(f"selected {selected} of {total} TSS regions")


if __name__ == "__main__":
    main()
