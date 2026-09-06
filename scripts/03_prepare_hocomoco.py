#!/usr/bin/env python3
"""Convert HOCOMOCO v11 human CORE annotations to standard sequence–PWM FASTA."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterator

import pandas as pd


REQUIRED_COLUMNS = ["Model", "UniProt AC"]


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


def fasta_accession(header: str) -> str:
    token = header.split()[0]
    parts = token.split("|")
    if len(parts) >= 3 and parts[0].lower() in {"sp", "tr"}:
        return parts[1]
    return parts[0]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--annotation", required=True, type=Path, help="HOCOMOCO annotation TSV")
    parser.add_argument("--sequences", required=True, type=Path, help="Human TF protein FASTA")
    parser.add_argument("--output", required=True, type=Path, help="Output FASTA")
    args = parser.parse_args()

    annotation = pd.read_csv(args.annotation, sep="\t", dtype=str)
    missing = [column for column in REQUIRED_COLUMNS if column not in annotation.columns]
    if missing:
        raise ValueError(f"Missing HOCOMOCO columns: {missing}")

    sequences: dict[str, str] = {}
    for header, sequence in read_fasta(args.sequences):
        accession = fasta_accession(header)
        if accession and sequence:
            if accession in sequences and sequences[accession] != sequence:
                raise ValueError(f"Conflicting sequences for accession {accession}")
            sequences[accession] = sequence

    records: list[tuple[str, str, str]] = []
    seen: set[tuple[str, str]] = set()
    missing_sequences: set[str] = set()
    for row in annotation[REQUIRED_COLUMNS].itertuples(index=False, name=None):
        model, accession = (str(value).strip() for value in row)
        if not model or not accession or model.lower() == "nan" or accession.lower() == "nan":
            continue
        key = (accession, model)
        if key in seen:
            continue
        seen.add(key)
        sequence = sequences.get(accession)
        if sequence is None:
            missing_sequences.add(accession)
            continue
        records.append((accession, model, sequence))

    with args.output.open("w", encoding="utf-8") as out:
        for accession, model, sequence in records:
            out.write(f">HOCOMOCO|{accession}|{model}\n{sequence}\n")

    print(f"Wrote {len(records)} HOCOMOCO sequence–PWM records to {args.output}")
    if missing_sequences:
        print(f"Skipped {len(missing_sequences)} accessions without a supplied protein sequence")


if __name__ == "__main__":
    main()

