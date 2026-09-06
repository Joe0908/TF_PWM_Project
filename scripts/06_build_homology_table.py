#!/usr/bin/env python3
"""Build sequence-evidence tiers from MMseqs2 results against the combined PWM reference."""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path
from typing import Iterator

import pandas as pd


MMSEQS_COLUMNS = [
    "query", "target", "pident", "alnlen", "mismatch", "gapopen",
    "qstart", "qend", "tstart", "tend", "evalue", "bits",
]
VALID_SOURCES = {"JASPAR", "CISBP", "HOCOMOCO"}


def fasta_ids(path: Path) -> Iterator[str]:
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            if raw.startswith(">"):
                yield normalize_query(raw[1:].strip().split()[0])


def normalize_query(value: object) -> str:
    token = str(value).strip().split()[0]
    parts = token.split("|")
    return parts[1] if len(parts) >= 3 and parts[0].lower() in {"sp", "tr"} else token


def parse_target(value: object) -> tuple[str, str, str]:
    token = str(value).strip().split()[0]
    parts = token.split("|", maxsplit=2)
    if len(parts) != 3 or parts[0].upper() not in VALID_SOURCES or not all(parts[1:]):
        raise ValueError(
            f"Invalid target identifier {token!r}; expected SOURCE|ACCESSION|MOTIF_ID"
        )
    return parts[0].upper(), parts[1], parts[2]


def add_unique(store: dict[str, list[str]], tf: str, pwm: str) -> None:
    if pwm not in store[tf]:
        store[tf].append(pwm)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--queries", required=True, type=Path, help="Query TF FASTA")
    parser.add_argument("--comparison", required=True, type=Path, help="MMseqs2 .comparison file")
    parser.add_argument("--output", required=True, type=Path, help="Wide annotation TSV")
    parser.add_argument("--provenance-output", type=Path, help="Optional long-format provenance TSV")
    parser.add_argument("--high-threshold", type=float, default=70.0)
    parser.add_argument("--lower-threshold", type=float, default=50.0)
    args = parser.parse_args()

    if not 0 <= args.lower_threshold < args.high_threshold <= 100:
        raise ValueError("Expected 0 <= lower threshold < high threshold <= 100")

    matches = pd.read_csv(
        args.comparison,
        sep="\t",
        header=None,
        names=MMSEQS_COLUMNS,
        dtype={"query": str, "target": str},
    )
    matches["pident"] = pd.to_numeric(matches["pident"], errors="raise")
    matches["bits"] = pd.to_numeric(matches["bits"], errors="coerce")
    if not matches.empty and matches["pident"].max() <= 1.0:
        raise ValueError(
            "The comparison file appears to contain fractional identity. "
            "Rerun 05_run_mmseqs.sh, which exports pident on a 0–100 scale."
        )

    matches["TF_name"] = matches["query"].map(normalize_query)
    parsed = matches["target"].map(parse_target)
    matches[["PWM_source", "reference_accession", "PWM"]] = pd.DataFrame(
        parsed.tolist(), index=matches.index
    )

    exact_tfs = set(matches.loc[matches["pident"].eq(100.0), "TF_name"])
    high_tfs = set(
        matches.loc[
            ~matches["TF_name"].isin(exact_tfs)
            & matches["pident"].ge(args.high_threshold),
            "TF_name",
        ]
    )

    identical: dict[str, list[str]] = defaultdict(list)
    homologous: dict[str, list[str]] = defaultdict(list)
    relative: dict[str, list[str]] = defaultdict(list)
    provenance_rows: list[dict[str, object]] = []

    for row in matches.itertuples(index=False):
        level: str | None = None
        if row.pident == 100.0:
            level = "Identical_PWM"
            add_unique(identical, row.TF_name, row.PWM)
        elif row.TF_name in high_tfs and row.pident >= args.high_threshold:
            level = "Homologous_PWM"
            add_unique(homologous, row.TF_name, row.PWM)
        elif (
            row.TF_name not in exact_tfs
            and row.TF_name not in high_tfs
            and row.pident >= args.lower_threshold
        ):
            level = "Relatively_Homologous_PWM"
            add_unique(relative, row.TF_name, row.PWM)

        if level is not None:
            provenance_rows.append(
                {
                    "TF_name": row.TF_name,
                    "annotation_level": level,
                    "PWM_source": row.PWM_source,
                    "reference_accession": row.reference_accession,
                    "PWM": row.PWM,
                    "pident": row.pident,
                    "bits": row.bits,
                }
            )

    rows = []
    for tf in sorted(dict.fromkeys(fasta_ids(args.queries))):
        rows.append(
            {
                "TF_name": tf,
                "TF_family": "",
                "Identical_PWM": ";".join(identical.get(tf, [])),
                "Homologous_PWM": ";".join(homologous.get(tf, [])),
                "Relatively_Homologous_PWM": ";".join(relative.get(tf, [])),
                "ModCRE": "",
                "AlphaFold": "",
            }
        )

    pd.DataFrame(rows).to_csv(args.output, sep="\t", index=False)
    if args.provenance_output:
        pd.DataFrame(
            provenance_rows,
            columns=[
                "TF_name", "annotation_level", "PWM_source", "reference_accession",
                "PWM", "pident", "bits",
            ],
        ).drop_duplicates().to_csv(args.provenance_output, sep="\t", index=False)
    print(f"Wrote {len(rows)} TF rows to {args.output}")


if __name__ == "__main__":
    main()

