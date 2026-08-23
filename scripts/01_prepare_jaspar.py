#!/usr/bin/env python3
"""Convert a JASPAR UniProt JSON export to a sequence–PWM FASTA file."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def motif_ids(raw: object) -> list[str]:
    if isinstance(raw, (list, tuple)):
        values: list[str] = []
        for item in raw:
            values.extend(motif_ids(item))
        return values
    if raw is None:
        return []
    value = str(raw).strip()
    return [value] if value else []


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="JASPAR UniProt JSON export")
    parser.add_argument("output", type=Path, help="Output FASTA")
    args = parser.parse_args()

    with args.input.open(encoding="utf-8") as handle:
        data = json.load(handle)

    written = 0
    with args.output.open("w", encoding="utf-8") as out:
        for accession, info in data.items():
            if not isinstance(info, (list, tuple)) or len(info) < 2:
                continue
            sequence = "".join(str(info[1]).split()).upper()
            if not sequence:
                continue
            for motif in dict.fromkeys(motif_ids(info[0])):
                out.write(f">{str(accession).strip()}|{motif}\n{sequence}\n")
                written += 1

    print(f"Wrote {written} JASPAR sequence–PWM records to {args.output}")


if __name__ == "__main__":
    main()
