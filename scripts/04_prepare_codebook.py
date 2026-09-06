#!/usr/bin/env python3
"""Convert the final Codebook supplementary tables to sequence–PWM records."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


PROTEIN_COLUMNS = [
    "TF",
    "TF category",
    "Ensembl ID",
    "DBD(s)",
    "DBD count",
    "Expert curated PWM ID",
    "Expert curated PWM experiment source",
    "Expert curated PWM experiment source type",
    "Expert curated PWM derivation method",
]
EXPERIMENT_COLUMNS = ["Experiment ID", "TF", "Plasmid ID"]
PLASMID_COLUMNS = ["Plasmid ID", "Insert Name", "TF"]
INSERT_COLUMNS = ["Insert Name", "TF", "Insert composition", "Amino acid sequence"]
METADATA_COLUMNS = [
    "HGNC_symbol",
    "Ensembl_gene_id",
    "DBD_family",
    "DBD_count",
    "PWM",
    "experiment_id",
    "assay",
    "derivation_method",
    "plasmid_id",
    "insert_name",
    "insert_composition",
    "amino_acid_sequence",
    "study_accession",
    "study_doi",
]
CENSUS_SOURCE_COLUMNS = [
    "TF (Ensembl ID)",
    "TF (HGNC)",
    "TF DBD(s) (Lambert 2018)",
    "TF DBD(s) (Codebook, Jolma 2026)",
    "Is TF?",
    "TF assessment (Lambert 2018)",
    "Conservative TF re-assessment",
    "TF re-assessment",
    "Motif status",
    "Number of representative motifs",
    "Representative motif ID(s) (see Supplementary Table 10, Supplementary Data 1)",
]
CENSUS_COLUMNS = [
    "Ensembl_gene_id",
    "HGNC_symbol",
    "DBD_Lambert2018",
    "DBD_Codebook2026",
    "is_TF",
    "assessment_Lambert2018",
    "conservative_reassessment",
    "reassessment",
    "motif_status",
    "representative_motif_count",
    "representative_motif_ids",
]


def read_sheet(path: Path, sheet: str, required: list[str]) -> pd.DataFrame:
    table = pd.read_excel(path, sheet_name=sheet, dtype=str)
    missing = [column for column in required if column not in table.columns]
    if missing:
        raise ValueError(f"{path.name}:{sheet} is missing columns: {missing}")
    return table


def indexed_rows(table: pd.DataFrame, key: str, values: set[str]) -> dict[str, pd.Series]:
    selected = table.loc[table[key].isin(values)].copy()
    duplicated = selected.loc[selected[key].duplicated(keep=False), key].dropna().unique()
    if len(duplicated):
        raise ValueError(f"Duplicate {key} values in selected records: {duplicated[:5].tolist()}")
    rows = {str(row[key]).strip(): row for _, row in selected.iterrows()}
    missing = sorted(values - set(rows))
    if missing:
        raise ValueError(f"Missing {key} values: {missing[:5]}")
    return rows


def clean(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--proteins", required=True, type=Path, help="Supplementary Table 1 XLSX")
    parser.add_argument("--inserts", required=True, type=Path, help="Supplementary Table 2 XLSX")
    parser.add_argument("--plasmids", required=True, type=Path, help="Supplementary Table 3 XLSX")
    parser.add_argument("--experiments", required=True, type=Path, help="Supplementary Table 4 XLSX")
    parser.add_argument("--output", required=True, type=Path, help="Output sequence–PWM FASTA")
    parser.add_argument("--metadata-output", type=Path, help="Optional normalized provenance TSV")
    parser.add_argument("--census", type=Path, help="Optional Supplementary Table 14 XLSX")
    parser.add_argument("--census-output", type=Path, help="Optional normalized TF-census TSV")
    parser.add_argument("--tf-category", default="Codebook TF")
    args = parser.parse_args()

    proteins = read_sheet(args.proteins, "TableS1", PROTEIN_COLUMNS)
    experiments = read_sheet(args.experiments, "TableS4", EXPERIMENT_COLUMNS)
    plasmids = read_sheet(args.plasmids, "TableS3", PLASMID_COLUMNS)
    inserts = read_sheet(args.inserts, "TableS2", INSERT_COLUMNS)

    selected = proteins.loc[
        proteins["TF category"].eq(args.tf_category)
        & proteins["Expert curated PWM ID"].notna()
        & proteins["Expert curated PWM ID"].astype(str).str.strip().ne("")
    ].copy()
    if selected["TF"].duplicated().any() or selected["Expert curated PWM ID"].duplicated().any():
        raise ValueError("Expected one distinct expert-curated PWM per selected TF")

    experiment_ids = set(selected["Expert curated PWM experiment source"].map(clean))
    experiment_rows = indexed_rows(experiments, "Experiment ID", experiment_ids)
    plasmid_ids = {clean(row["Plasmid ID"]) for row in experiment_rows.values()}
    plasmid_rows = indexed_rows(plasmids, "Plasmid ID", plasmid_ids)
    insert_names = {clean(row["Insert Name"]) for row in plasmid_rows.values()}
    insert_rows = indexed_rows(inserts, "Insert Name", insert_names)

    records: list[dict[str, str]] = []
    for _, protein in selected.sort_values(["TF", "Expert curated PWM ID"]).iterrows():
        gene = clean(protein["TF"])
        ensembl = clean(protein["Ensembl ID"])
        pwm = clean(protein["Expert curated PWM ID"])
        experiment_id = clean(protein["Expert curated PWM experiment source"])
        experiment = experiment_rows[experiment_id]
        plasmid_id = clean(experiment["Plasmid ID"])
        plasmid = plasmid_rows[plasmid_id]
        insert_name = clean(plasmid["Insert Name"])
        insert = insert_rows[insert_name]

        linked_genes = {clean(experiment["TF"]), clean(plasmid["TF"]), clean(insert["TF"])} - {""}
        if linked_genes != {gene}:
            raise ValueError(
                f"Inconsistent TF linkage for {gene}/{pwm}: {sorted(linked_genes)}"
            )
        sequence = "".join(clean(insert["Amino acid sequence"]).split()).upper()
        if not sequence or set(sequence) - set("ACDEFGHIKLMNPQRSTVWY"):
            raise ValueError(f"Invalid amino-acid sequence for {gene}/{pwm}")
        if not ensembl or not pwm or "|" in ensembl or "|" in pwm:
            raise ValueError(f"Invalid reference identifier for {gene}/{pwm}")

        records.append(
            {
                "HGNC_symbol": gene,
                "Ensembl_gene_id": ensembl,
                "DBD_family": clean(protein["DBD(s)"]),
                "DBD_count": clean(protein["DBD count"]),
                "PWM": pwm,
                "experiment_id": experiment_id,
                "assay": clean(protein["Expert curated PWM experiment source type"]),
                "derivation_method": clean(protein["Expert curated PWM derivation method"]),
                "plasmid_id": plasmid_id,
                "insert_name": insert_name,
                "insert_composition": clean(insert["Insert composition"]),
                "amino_acid_sequence": sequence,
                "study_accession": "Jolma2026a",
                "study_doi": "10.1038/s41586-026-10798-9",
            }
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(
                f">CODEBOOK|{record['Ensembl_gene_id']}|{record['PWM']}\n"
            )
            sequence = record["amino_acid_sequence"]
            for start in range(0, len(sequence), 60):
                handle.write(sequence[start : start + 60] + "\n")

    if args.metadata_output:
        args.metadata_output.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(records, columns=METADATA_COLUMNS).to_csv(
            args.metadata_output, sep="\t", index=False
        )
    if bool(args.census) != bool(args.census_output):
        raise ValueError("--census and --census-output must be supplied together")
    if args.census and args.census_output:
        census = read_sheet(args.census, "TableS14", CENSUS_SOURCE_COLUMNS)
        census = census[CENSUS_SOURCE_COLUMNS].rename(
            columns=dict(zip(CENSUS_SOURCE_COLUMNS, CENSUS_COLUMNS))
        )
        for column in CENSUS_COLUMNS:
            census[column] = census[column].map(clean)
        if census["Ensembl_gene_id"].eq("").any() or census["Ensembl_gene_id"].duplicated().any():
            raise ValueError("TF census Ensembl_gene_id values must be populated and unique")
        args.census_output.parent.mkdir(parents=True, exist_ok=True)
        census.replace("", "NA").to_csv(args.census_output, sep="\t", index=False)
    print(f"Wrote {len(records)} Codebook sequence–PWM records to {args.output}")


if __name__ == "__main__":
    main()
