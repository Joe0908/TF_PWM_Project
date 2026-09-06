#!/usr/bin/env python3
"""Convert Cis-BP 2.00 SQL exports to standard sequence–PWM FASTA."""

from __future__ import annotations

import argparse
import re
from collections import defaultdict
from pathlib import Path


MOTIF_PATTERN = re.compile(
    r"\('([^']+)',\s*'([^']+)',\s*[^,]+,\s*[^,]+,\s*[^,]+,\s*[^,]+,\s*'([^']+)'"
)
PROTEIN_PATTERN = re.compile(
    r"\('([^']+)',\s*'([^']+)',\s*[^,]+,\s*[^,]+,\s*'([^']+)'\)"
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--motifs", required=True, type=Path, help="CisBP_2.00.all.motifs.sql")
    parser.add_argument("--proteins", required=True, type=Path, help="CisBP_2.00.all.proteins.sql")
    parser.add_argument(
        "--accession-map",
        type=Path,
        help="Optional TSV with TF_ID and accession columns for stable-ID direct matching",
    )
    parser.add_argument("--output", required=True, type=Path, help="Output FASTA")
    args = parser.parse_args()

    accession_map: dict[str, str] = {}
    if args.accession_map:
        import pandas as pd

        crosswalk = pd.read_csv(args.accession_map, sep="\t", dtype=str)
        required = ["TF_ID", "accession"]
        missing = [column for column in required if column not in crosswalk.columns]
        if missing:
            raise ValueError(f"Cis-BP accession map is missing columns: {missing}")
        crosswalk = crosswalk[required].dropna()
        crosswalk["TF_ID"] = crosswalk["TF_ID"].str.strip()
        crosswalk["accession"] = crosswalk["accession"].str.strip()
        if crosswalk["TF_ID"].duplicated().any():
            raise ValueError("Cis-BP accession map contains duplicate TF_ID values")
        accession_map = dict(zip(crosswalk["TF_ID"], crosswalk["accession"]))

    motif_map: dict[str, list[str]] = defaultdict(list)
    motif_text = args.motifs.read_text(encoding="utf-8", errors="replace")
    for motif_id, tf_id, iupac in MOTIF_PATTERN.findall(motif_text):
        if iupac.upper() != "NULL" and motif_id not in motif_map[tf_id]:
            motif_map[tf_id].append(motif_id.strip())

    protein_map: dict[str, str] = {}
    protein_text = args.proteins.read_text(encoding="utf-8", errors="replace")
    for _, tf_id, sequence in PROTEIN_PATTERN.findall(protein_text):
        sequence = "".join(sequence.split()).upper()
        if sequence and sequence != "NULL":
            protein_map[tf_id] = sequence

    written = 0
    with args.output.open("w", encoding="utf-8") as out:
        for tf_id, sequence in protein_map.items():
            for motif_id in motif_map.get(tf_id, []):
                accession = accession_map.get(tf_id, tf_id)
                out.write(f">CISBP|{accession}|{motif_id}\n{sequence}\n")
                written += 1

    print(f"Wrote {written} Cis-BP sequence–PWM records to {args.output}")


if __name__ == "__main__":
    main()
