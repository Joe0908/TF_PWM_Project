#!/usr/bin/env python3
"""Attach the supplied accession-to-family mapping to a TF annotation table."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--mapping", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    table = pd.read_csv(args.input, sep="\t", dtype=str)
    mapping = pd.read_csv(
        args.mapping,
        header=None,
        names=["TF_name", "TF_family_mapped"],
        dtype=str,
    )
    mapping["TF_name"] = mapping["TF_name"].str.strip()
    if mapping["TF_name"].duplicated().any():
        duplicates = mapping.loc[mapping["TF_name"].duplicated(), "TF_name"].tolist()
        raise ValueError(f"Duplicate TF accessions in family mapping: {duplicates[:5]}")

    table = table.drop(columns=["TF_family"], errors="ignore")
    table = table.merge(mapping, on="TF_name", how="left", validate="one_to_one")
    table = table.rename(columns={"TF_family_mapped": "TF_family"})
    columns = table.columns.tolist()
    columns.insert(1, columns.pop(columns.index("TF_family")))
    table[columns].to_csv(args.output, sep="\t", index=False)
    print(f"Attached family annotations to {table['TF_family'].notna().sum()} of {len(table)} TFs")


if __name__ == "__main__":
    main()

