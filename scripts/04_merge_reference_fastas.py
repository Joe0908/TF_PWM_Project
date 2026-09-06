#!/usr/bin/env python3
"""Merge standard direct-PWM FASTA files while preserving source provenance."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterator


VALID_SOURCES = {"JASPAR", "CISBP", "HOCOMOCO"}


def read_fasta(path: Path) -> Iterator[tuple[str, str]]:
    header: str | None = None
    sequence: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence).upper()
                header, sequence = line[1:].strip(), []
            else:
                sequence.append(line)
    if header is not None:
        yield header, "".join(sequence).upper()


def validate_header(header: str) -> None:
    parts = header.split()[0].split("|", maxsplit=2)
    if len(parts) != 3 or parts[0].upper() not in VALID_SOURCES or not all(parts[1:]):
        raise ValueError(
            f"Invalid reference header {header!r}; expected SOURCE|ACCESSION|MOTIF_ID"
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", type=Path, help="Standard source FASTA files")
    parser.add_argument("--output", required=True, type=Path, help="Merged FASTA")
    args = parser.parse_args()

    seen: set[tuple[str, str]] = set()
    written = 0
    with args.output.open("w", encoding="utf-8") as out:
        for path in args.inputs:
            for header, sequence in read_fasta(path):
                validate_header(header)
                key = (header.split()[0], sequence)
                if not sequence or key in seen:
                    continue
                seen.add(key)
                out.write(f">{header}\n")
                for start in range(0, len(sequence), 60):
                    out.write(sequence[start : start + 60] + "\n")
                written += 1

    print(f"Wrote {written} distinct direct-PWM reference records to {args.output}")


if __name__ == "__main__":
    main()

