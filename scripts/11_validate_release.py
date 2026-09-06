#!/usr/bin/env python3
"""Validate the schema, hierarchy, Codebook provenance and release snapshot."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd


EVIDENCE_COLUMNS = [
    "Direct_PWM",
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
    "Direct_PWM": 2395,
    "Homologous_PWM": 1741,
    "Relatively_Homologous_PWM": 602,
    "ModCRE": 372,
    "AlphaFold": 234,
    "Unannotated": 121,
}
CODEBOOK_MOTIF_COLUMNS = [
    "HGNC_symbol", "Ensembl_gene_id", "DBD_family", "DBD_count", "PWM",
    "experiment_id", "assay", "derivation_method", "plasmid_id", "insert_name",
    "insert_composition", "amino_acid_sequence", "study_accession", "study_doi",
]
CODEBOOK_MAP_COLUMNS = [
    "TF_name", "PWM_source", "reference_accession", "PWM", "match_basis",
    "HGNC_symbol", "Ensembl_gene_id", "DBD_family", "DBD_count",
    "sequence_relation", "canonical_accession",
]
CODEBOOK_CENSUS_COUNTS = {
    "Motif (measured) - Lambert 2018": 1107,
    "Motif (inferred) - Lambert 2018": 104,
    "New motif - Codebook": 177,
    "New motif - literature": 33,
    "No motif - Codebook tested": 146,
    "No data": 72,
}


def populated(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip().ne("")


def motif_set(value: object) -> set[str]:
    if pd.isna(value):
        return set()
    return {item.strip() for item in str(value).split(";") if item.strip()}


def validate_codebook(
    release: pd.DataFrame,
    motifs_path: Path,
    map_path: Path,
    census_path: Path,
    errors: list[str],
) -> None:
    motifs = pd.read_csv(motifs_path, sep="\t", dtype=str).fillna("")
    assignments = pd.read_csv(map_path, sep="\t", dtype=str).fillna("")
    census = pd.read_csv(census_path, sep="\t", dtype=str).fillna("")

    if list(motifs.columns) != CODEBOOK_MOTIF_COLUMNS:
        errors.append("Codebook motif metadata has an unexpected schema")
        return
    if list(assignments.columns) != CODEBOOK_MAP_COLUMNS:
        errors.append("Codebook accession map has an unexpected schema")
        return
    if len(motifs) != 177 or motifs["HGNC_symbol"].nunique() != 177:
        errors.append("Codebook motif metadata must contain 177 distinct Codebook TFs")
    if motifs["PWM"].nunique() != 177 or motifs["experiment_id"].nunique() != 177:
        errors.append("Codebook PWM and experiment IDs must each be unique")
    if set(motifs["study_accession"]) != {"Jolma2026a"}:
        errors.append("Codebook study accession must be Jolma2026a")
    if set(motifs["study_doi"]) != {"10.1038/s41586-026-10798-9"}:
        errors.append("Codebook DOI is inconsistent")
    valid_aa = motifs["amino_acid_sequence"].str.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+")
    if not valid_aa.all():
        errors.append("Codebook metadata contains an invalid amino-acid sequence")

    motif_pairs = set(zip(motifs["Ensembl_gene_id"], motifs["PWM"]))
    map_pairs = set(zip(assignments["reference_accession"], assignments["PWM"]))
    if not map_pairs <= motif_pairs:
        errors.append("Codebook accession map refers to unknown Ensembl/PWM pairs")
    if set(assignments["PWM_source"]) != {"CODEBOOK"}:
        errors.append("Codebook accession map has an invalid PWM_source")
    if set(assignments["match_basis"]) - {"exact_construct", "stable_gene_id"}:
        errors.append("Codebook accession map has an invalid match_basis")
    canonical = assignments.loc[assignments["canonical_accession"].eq("yes")]
    if len(canonical) != 177 or canonical["HGNC_symbol"].nunique() != 177:
        errors.append("Codebook accession map must designate one canonical accession per motif TF")
    if canonical["HGNC_symbol"].duplicated().any():
        errors.append("A Codebook TF has more than one canonical accession")

    indexed = release.set_index("TF_name")
    missing_accessions = sorted(set(assignments["TF_name"]) - set(indexed.index))
    if missing_accessions:
        errors.append(f"Codebook assignments absent from release: {missing_accessions[:5]}")
    else:
        for row in assignments.itertuples(index=False):
            if row.PWM not in motif_set(indexed.at[row.TF_name, "Direct_PWM"]):
                errors.append(f"Codebook PWM {row.PWM} is not direct for {row.TF_name}")
                break
            if indexed.at[row.TF_name, "Best_annotation_level"] != "Direct_PWM":
                errors.append(f"Codebook assignment is not retained as direct for {row.TF_name}")
                break

    if len(census) != 1639 or census["Ensembl_gene_id"].nunique() != 1639:
        errors.append("Codebook TF census must contain 1,639 distinct Ensembl records")
    if census["is_TF"].value_counts().to_dict() != {"Yes": 1638, "No": 1}:
        errors.append("Codebook TF census assessment counts are inconsistent")
    if census["motif_status"].value_counts().to_dict() != CODEBOOK_CENSUS_COUNTS:
        errors.append("Codebook TF census motif-status counts are inconsistent")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("table", type=Path)
    parser.add_argument(
        "--structural-exclusions", nargs=2, metavar=("MODCRE_LIST", "AF3_LIST"),
        type=Path, help="Optionally confirm excluded structural models are absent",
    )
    parser.add_argument(
        "--codebook-reference", nargs=3,
        metavar=("MOTIFS_TSV", "ACCESSION_MAP_TSV", "TF_CENSUS_TSV"), type=Path,
        help="Validate Codebook metadata, direct assignments and census",
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

    expected_rows = sum(EXPECTED_COUNTS.values())
    if len(table) != expected_rows:
        errors.append(f"expected {expected_rows:,} rows, found {len(table):,}")
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

    identity_pattern = re.compile(r"^.+ \(([0-9]+(?:\.[0-9]+)?)%\)$")
    for column in ("Homologous_PWM", "Relatively_Homologous_PWM"):
        for value in table[column].dropna():
            identities = []
            for item in value.split(";"):
                match = identity_pattern.fullmatch(item.strip())
                if not match:
                    errors.append(f"invalid homology item in {column}: {item}")
                    break
                identities.append(float(match.group(1)))
            if identities and max(identities) - min(identities) > 0.05:
                errors.append(f"non-nearest homology motifs retained in {column}")
                break

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

    if args.codebook_reference:
        validate_codebook(table, *args.codebook_reference, errors)

    if errors:
        raise SystemExit("Validation failed:\n- " + "\n- ".join(errors))

    print(f"PASS: {args.table}")
    print(f"Rows: {len(table):,}; unique TFs: {table['TF_name'].nunique():,}")
    for level, count in EXPECTED_COUNTS.items():
        print(f"{level}: {count:,}")


if __name__ == "__main__":
    main()
