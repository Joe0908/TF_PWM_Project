#!/usr/bin/env python3
"""Attach filtered structural models only where sequence-based PWM evidence is absent."""

from __future__ import annotations

import argparse
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd


SEQUENCE_COLUMNS = ["Identical_PWM", "Homologous_PWM", "Relatively_Homologous_PWM"]


def populated(value: object) -> bool:
    return not pd.isna(value) and str(value).strip() not in {"", "nan", "None"}


def load_exclusions(path: Path) -> set[str]:
    excluded: set[str] = set()
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            for token in re.split(r"\s+", raw.strip()):
                if token:
                    excluded.add(token if token.endswith(".meme") else token + ".meme")
    return excluded


def done_models(path: Path, excluded: set[str]) -> list[str]:
    models: list[str] = []
    seen: set[str] = set()
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 2 or not fields[1].strip().startswith("DONE"):
                continue
            model = Path(fields[0].strip()).name
            if model in excluded or model in seen:
                continue
            seen.add(model)
            models.append(model)
    return models


def append_unique(store: dict[str, list[str]], key: str, value: str) -> None:
    if value not in store[key]:
        store[key].append(value)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--modcre-log", required=True, type=Path)
    parser.add_argument("--af3-log", required=True, type=Path)
    parser.add_argument("--modcre-exclusions", required=True, type=Path)
    parser.add_argument("--af3-exclusions", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    table = pd.read_csv(args.input, sep="\t", dtype=str)
    missing = [column for column in SEQUENCE_COLUMNS if column not in table.columns]
    if missing:
        raise ValueError(f"Missing sequence-evidence columns: {missing}")

    accessions = set(table["TF_name"].astype(str).str.strip())
    accessions_upper = {value.upper(): value for value in accessions}

    modcre_map: dict[str, list[str]] = defaultdict(list)
    for model in done_models(args.modcre_log, load_exclusions(args.modcre_exclusions)):
        for part in model.split("_"):
            accession = accessions_upper.get(part.upper())
            if accession is not None:
                append_unique(modcre_map, accession, model)
                break

    af3_map: dict[str, list[str]] = defaultdict(list)
    for model in done_models(args.af3_log, load_exclusions(args.af3_exclusions)):
        accession = accessions_upper.get(model.split(".", maxsplit=1)[0].upper())
        if accession is not None:
            append_unique(af3_map, accession, model)

    def has_sequence_evidence(row: pd.Series) -> bool:
        return any(populated(row[column]) for column in SEQUENCE_COLUMNS)

    modcre_values: list[str] = []
    af3_values: list[str] = []
    for _, row in table.iterrows():
        tf = row["TF_name"]
        if has_sequence_evidence(row):
            modcre_values.append("")
            af3_values.append("")
            continue
        modcre = ";".join(modcre_map.get(tf, []))
        modcre_values.append(modcre)
        af3_values.append("" if modcre else ";".join(af3_map.get(tf, [])))

    table["ModCRE"] = modcre_values
    table["AlphaFold"] = af3_values
    table.to_csv(args.output, sep="\t", index=False)
    print(f"Retained ModCRE models for {sum(bool(value) for value in modcre_values)} TFs")
    print(f"Retained AlphaFold-derived models for {sum(bool(value) for value in af3_values)} TFs")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()

