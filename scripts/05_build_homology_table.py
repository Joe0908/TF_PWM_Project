#!/usr/bin/env python3
"""Build the identity/homology annotation table from explicit MMseqs2 pident output."""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path
from typing import Iterator

import pandas as pd


MMSEQS_COLUMNS = [
    "query",
    "target",
    "pident",
    "alnlen",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "tstart",
    "tend",
    "evalue",
    "bits",
]


def fasta_ids(path: Path) -> Iterator[str]:
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            if raw.startswith(">"):
                yield normalize_query(raw[1:].strip().split()[0])


def normalize_query(value: object) -> str:
    token = str(value).strip().split()[0]
    parts = token.split("|")
    return parts[1] if len(parts) >= 3 and parts[1] else token


def motif_id(value: object) -> str | None:
    token = str(value).strip().split()[0]
    parts = token.split("|", maxsplit=1)
    return parts[1] if len(parts) == 2 and parts[1] else None


def add_unique(store: dict[str, list[str]], tf: str, pwm: str | None) -> None:
    if pwm and pwm not in store[tf]:
        store[tf].append(pwm)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--queries", required=True, type=Path, help="Query TF FASTA")
    parser.add_argument("--comparison", required=True, type=Path, help="MMseqs2 .comparison file")
    parser.add_argument("--output", required=True, type=Path, help="Output TSV")
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
    if matches["pident"].max() <= 1.0:
        raise ValueError(
            "The comparison file appears to contain fractional identity. "
            "Rerun 04_run_mmseqs.sh, which exports pident on a 0–100 scale."
        )
    matches["TF_name"] = matches["query"].map(normalize_query)
    matches["PWM"] = matches["target"].map(motif_id)

    identical: dict[str, list[str]] = defaultdict(list)
    homologous: dict[str, list[str]] = defaultdict(list)
    relative: dict[str, list[str]] = defaultdict(list)

    exact_tfs = set(matches.loc[matches["pident"] == 100.0, "TF_name"])
    for row in matches.itertuples(index=False):
        tf, pwm, identity = row.TF_name, row.PWM, row.pident
        if identity == 100.0:
            add_unique(identical, tf, pwm)
        elif tf not in exact_tfs and identity >= args.high_threshold:
            add_unique(homologous, tf, pwm)
        elif tf not in exact_tfs and identity >= args.lower_threshold:
            add_unique(relative, tf, pwm)

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
    print(f"Wrote {len(rows)} TF rows to {args.output}")


if __name__ == "__main__":
    main()
