#!/usr/bin/env python3
"""Integrate HOCOMOCO v11 by UniProt accession and enforce evidence hierarchy."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


FINAL_COLUMNS = [
    "TF_name",
    "TF_family",
    "Identical_PWM",
    "Homologous_PWM",
    "Relatively_Homologous_PWM",
    "ModCRE",
    "AlphaFold",
]
HOCO_COLUMNS = ["Model", "Transcription factor", "TF family", "UniProt AC", "UniProt ID"]
EVIDENCE_COLUMNS = FINAL_COLUMNS[2:]


def is_missing(value: object) -> bool:
    if pd.isna(value):
        return True
    text = str(value).strip()
    return text == "" or text.lower() in {"nan", "none", "na", "n/a", "-"}


def clean_join_unique(values: Iterable[object]) -> object:
    seen: set[str] = set()
    output: list[str] = []
    for value in values:
        if is_missing(value):
            continue
        for item in str(value).split(";"):
            item = item.strip()
            if item and item not in seen:
                seen.add(item)
                output.append(item)
    return ";".join(output) if output else np.nan


def best_level(row: pd.Series) -> str:
    for column in EVIDENCE_COLUMNS:
        if not is_missing(row[column]):
            return column
    return "Unannotated"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path, help="Pre-HOCOMOCO TF chart")
    parser.add_argument("--hocomoco", required=True, type=Path, help="HOCOMOCO annotation TSV")
    parser.add_argument("--output", required=True, type=Path, help="Integrated output TSV")
    parser.add_argument("--audit", required=True, type=Path, help="Integration audit TSV")
    args = parser.parse_args()

    original = pd.read_csv(args.input, sep="\t", dtype=str)
    hocomoco = pd.read_csv(args.hocomoco, sep="\t", dtype=str)
    missing_final = [column for column in FINAL_COLUMNS if column not in original.columns]
    missing_hoco = [column for column in HOCO_COLUMNS if column not in hocomoco.columns]
    if missing_final:
        raise ValueError(f"Missing input columns: {missing_final}")
    if missing_hoco:
        raise ValueError(f"Missing HOCOMOCO columns: {missing_hoco}")

    original["TF_name"] = original["TF_name"].astype(str).str.strip()
    hocomoco["UniProt AC"] = hocomoco["UniProt AC"].fillna("").astype(str).str.strip()
    hocomoco = hocomoco[hocomoco["UniProt AC"] != ""].copy()

    grouped = (
        hocomoco.groupby("UniProt AC", as_index=False)
        .agg(
            {
                "Model": clean_join_unique,
                "Transcription factor": clean_join_unique,
                "TF family": clean_join_unique,
                "UniProt ID": clean_join_unique,
            }
        )
        .rename(
            columns={
                "UniProt AC": "TF_name",
                "Model": "HOCOMOCO_PWM",
                "Transcription factor": "HOCOMOCO_TF_name",
                "TF family": "HOCOMOCO_family",
                "UniProt ID": "HOCOMOCO_UniProt_ID",
            }
        )
    )
    lookup = grouped.set_index("TF_name").to_dict(orient="index")
    integrated = original.copy()
    audit_rows: list[dict[str, object]] = []

    for index, row in integrated.iterrows():
        tf = row["TF_name"]
        if tf not in lookup:
            continue
        hoco = lookup[tf]
        old_identical = row["Identical_PWM"]
        new_identical = clean_join_unique([old_identical, hoco["HOCOMOCO_PWM"]])
        integrated.at[index, "Identical_PWM"] = new_identical
        audit_rows.append(
            {
                "TF_name": tf,
                "action": (
                    "upgrade_to_identical_pwm"
                    if is_missing(old_identical)
                    else "add_hocomoco_to_existing_identical_pwm"
                ),
                "old_Identical_PWM": old_identical,
                "new_Identical_PWM": new_identical,
                "old_Homologous_PWM": row["Homologous_PWM"],
                "old_Relatively_Homologous_PWM": row["Relatively_Homologous_PWM"],
                "old_ModCRE": row["ModCRE"],
                "old_AlphaFold": row["AlphaFold"],
                "HOCOMOCO_TF_name": hoco["HOCOMOCO_TF_name"],
                "HOCOMOCO_family": hoco["HOCOMOCO_family"],
                "HOCOMOCO_UniProt_ID": hoco["HOCOMOCO_UniProt_ID"],
            }
        )

    existing = set(integrated["TF_name"])
    new_rows: list[dict[str, object]] = []
    for _, hoco in grouped.loc[~grouped["TF_name"].isin(existing)].iterrows():
        new_rows.append(
            {
                "TF_name": hoco["TF_name"],
                "TF_family": hoco["HOCOMOCO_family"],
                "Identical_PWM": hoco["HOCOMOCO_PWM"],
                "Homologous_PWM": np.nan,
                "Relatively_Homologous_PWM": np.nan,
                "ModCRE": np.nan,
                "AlphaFold": np.nan,
            }
        )
        audit_rows.append(
            {
                "TF_name": hoco["TF_name"],
                "action": "new_tf_from_hocomoco",
                "old_Identical_PWM": np.nan,
                "new_Identical_PWM": hoco["HOCOMOCO_PWM"],
                "old_Homologous_PWM": np.nan,
                "old_Relatively_Homologous_PWM": np.nan,
                "old_ModCRE": np.nan,
                "old_AlphaFold": np.nan,
                "HOCOMOCO_TF_name": hoco["HOCOMOCO_TF_name"],
                "HOCOMOCO_family": hoco["HOCOMOCO_family"],
                "HOCOMOCO_UniProt_ID": hoco["HOCOMOCO_UniProt_ID"],
            }
        )
    if new_rows:
        integrated = pd.concat([integrated, pd.DataFrame(new_rows)], ignore_index=True)

    for index, row in integrated.iterrows():
        retained = False
        for column in EVIDENCE_COLUMNS:
            if retained:
                integrated.at[index, column] = np.nan
            elif not is_missing(row[column]):
                retained = True

    integrated["Best_annotation_level"] = integrated.apply(best_level, axis=1)
    integrated["Best_PWM_or_model"] = integrated.apply(
        lambda row: np.nan
        if row["Best_annotation_level"] == "Unannotated"
        else row[row["Best_annotation_level"]],
        axis=1,
    )
    integrated["N_nonempty_annotation_columns"] = integrated.apply(
        lambda row: sum(not is_missing(row[column]) for column in EVIDENCE_COLUMNS), axis=1
    )

    audit = pd.DataFrame(audit_rows)
    integrated.to_csv(args.output, sep="\t", index=False)
    audit.to_csv(args.audit, sep="\t", index=False)

    print(f"Original rows: {len(original)}")
    print(f"HOCOMOCO unique accessions: {len(grouped)}")
    print(f"New HOCOMOCO TFs: {len(new_rows)}")
    print(f"Integrated rows: {len(integrated)}")
    print(integrated["Best_annotation_level"].value_counts().to_string())


if __name__ == "__main__":
    main()
