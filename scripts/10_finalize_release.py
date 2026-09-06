#!/usr/bin/env python3
"""Enforce the evidence hierarchy and add release summary fields."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


EVIDENCE_COLUMNS = [
    "Direct_PWM",
    "Homologous_PWM",
    "Relatively_Homologous_PWM",
    "ModCRE",
    "AlphaFold",
]
BASE_COLUMNS = ["TF_name", "TF_family", *EVIDENCE_COLUMNS]


def populated(value: object) -> bool:
    if pd.isna(value):
        return False
    return str(value).strip().lower() not in {"", "nan", "none", "na", "n/a", "-"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    table = pd.read_csv(args.input, sep="\t", dtype=str)
    missing = [column for column in BASE_COLUMNS if column not in table.columns]
    if missing:
        raise ValueError(f"Missing input columns: {missing}")
    if table["TF_name"].isna().any() or table["TF_name"].duplicated().any():
        raise ValueError("TF_name must be populated and unique")

    table = table[BASE_COLUMNS].copy()
    best_levels: list[str] = []
    best_values: list[object] = []
    nonempty_counts: list[int] = []

    for index, row in table.iterrows():
        retained_level = "Unannotated"
        retained_value: object = pd.NA
        for column in EVIDENCE_COLUMNS:
            if retained_level == "Unannotated" and populated(row[column]):
                retained_level = column
                retained_value = row[column]
            elif retained_level != "Unannotated":
                table.at[index, column] = pd.NA
        best_levels.append(retained_level)
        best_values.append(retained_value)
        nonempty_counts.append(0 if retained_level == "Unannotated" else 1)

    table["Best_annotation_level"] = best_levels
    table["Best_PWM_or_model"] = best_values
    table["N_nonempty_annotation_columns"] = nonempty_counts
    table.to_csv(args.output, sep="\t", index=False)
    print(f"Wrote {len(table)} TF rows to {args.output}")
    print(table["Best_annotation_level"].value_counts().to_string())


if __name__ == "__main__":
    main()
