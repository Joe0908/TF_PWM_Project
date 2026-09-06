#!/usr/bin/env python3
"""Validate the schema, hierarchy and expected counts of the 5,417-TF release."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


EVIDENCE_COLUMNS = [
    "Identical_PWM",
    "Homologous_PWM",
    "Relatively_Homologous_PWM",
    "ModCRE",
    "AlphaFold",
]
REQUIRED_COLUMNS = [
    "TF_name", "TF_family", *EVIDENCE_COLUMNS,
    "Best_annotation_level", "Best_PWM_or_model", "N_nonempty_annotation_columns",
]
EXPECTED_COUNTS = {
    "Identical_PWM": 2160,
    "Homologous_PWM": 1786,
    "Relatively_Homologous_PWM": 700,
    "ModCRE": 400,
    "AlphaFold": 242,
    "Unannotated": 129,
}


def populated(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip().ne("")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("table", type=Path)
    parser.add_argument(
        "--structural-exclusions",
        nargs=2,
        metavar=("MODCRE_LIST", "AF3_LIST"),
        type=Path,
        help="Optionally confirm excluded structural models are absent",
    )
    args = parser.parse_args()

    table = pd.read_csv(args.table, sep="\t", dtype=str)
    errors: list[str] = []
    missing = [column for column in REQUIRED_COLUMNS if column not in table.columns]
    extra = [column for column in table.columns if column not in REQUIRED_COLUMNS]
    if missing:
        errors.append(f"missing columns: {missing}")
    if extra:
        errors.append(f"unexpected columns: {extra}")
    if errors:
        raise SystemExit("Validation failed:\n- " + "\n- ".join(errors))

    if len(table) != 5417:
        errors.append(f"expected 5,417 rows, found {len(table):,}")
    if table["TF_name"].isna().any() or populated(table["TF_name"]).sum() != len(table):
        errors.append("TF_name contains missing or blank values")
    duplicated = table["TF_name"].duplicated().sum()
    if duplicated:
        errors.append(f"found {duplicated} duplicated TF_name values")

    evidence_count = sum(populated(table[column]).astype(int) for column in EVIDENCE_COLUMNS)
    stored_count = pd.to_numeric(table["N_nonempty_annotation_columns"], errors="coerce")
    if not evidence_count.equals(stored_count):
        errors.append("N_nonempty_annotation_columns does not match the evidence columns")
    if not evidence_count.isin([0, 1]).all():
        errors.append("at least one TF retains more than one evidence level")

    expected_level = pd.Series("Unannotated", index=table.index)
    for column in reversed(EVIDENCE_COLUMNS):
        expected_level.loc[populated(table[column])] = column
    if not expected_level.equals(table["Best_annotation_level"]):
        errors.append("Best_annotation_level is inconsistent with retained evidence")

    expected_model = pd.Series(pd.NA, index=table.index, dtype="object")
    for column in EVIDENCE_COLUMNS:
        mask = table["Best_annotation_level"].eq(column)
        expected_model.loc[mask] = table.loc[mask, column]
    actual_model = table["Best_PWM_or_model"].replace("", pd.NA)
    if not expected_model.fillna("").astype(str).eq(actual_model.fillna("").astype(str)).all():
        errors.append("Best_PWM_or_model is inconsistent with Best_annotation_level")

    observed_counts = table["Best_annotation_level"].value_counts().to_dict()
    if observed_counts != EXPECTED_COUNTS:
        errors.append(f"unexpected annotation-level counts: {observed_counts}")

    if args.structural_exclusions:
        invalid: set[str] = set()
        for path in args.structural_exclusions:
            for raw in path.read_text(encoding="utf-8").splitlines():
                token = raw.strip()
                if token:
                    invalid.add(token if token.endswith(".meme") else token + ".meme")
        present = set()
        for column in ("ModCRE", "AlphaFold"):
            for value in table[column].dropna():
                present.update(item.strip() for item in value.split(";") if item.strip())
        overlap = sorted(invalid & present)
        if overlap:
            errors.append(f"found {len(overlap)} excluded structural models in release")

    if errors:
        raise SystemExit("Validation failed:\n- " + "\n- ".join(errors))

    print(f"PASS: {args.table}")
    print(f"Rows: {len(table):,}; unique TFs: {table['TF_name'].nunique():,}")
    for level, count in EXPECTED_COUNTS.items():
        print(f"{level}: {count:,}")


if __name__ == "__main__":
    main()

