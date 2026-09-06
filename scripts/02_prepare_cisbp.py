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
    parser.add_argument("--output", required=True, type=Path, help="Output FASTA")
    args = parser.parse_args()

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
                out.write(f">CISBP|{tf_id}|{motif_id}\n{sequence}\n")
                written += 1

    print(f"Wrote {written} Cis-BP sequence–PWM records to {args.output}")


if __name__ == "__main__":
    main()

