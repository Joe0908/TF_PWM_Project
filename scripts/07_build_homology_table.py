#!/usr/bin/env python3
"""Build direct and nearest-neighbour PWM evidence from MMseqs2 results."""

from __future__ import annotations

import argparse
import re
from collections import defaultdict
from pathlib import Path
from typing import Iterator

import pandas as pd


MMSEQS_COLUMNS = [
    "query", "target", "pident", "alnlen", "mismatch", "gapopen",
    "qstart", "qend", "tstart", "tend", "evalue", "bits",
    "qcov", "tcov", "qlen", "tlen",
]
VALID_SOURCES = {"JASPAR", "CISBP", "HOCOMOCO", "CODEBOOK"}
VALID_MATCH_BASES = {
    "stable_accession", "stable_gene_id", "exact_construct", "exact_sequence"
}
DIRECT_COLUMNS = [
    "TF_name", "PWM_source", "reference_accession", "PWM", "match_basis"
]
DOMAIN_COLUMNS = ["identifier", "DBD_family", "DBD_count"]


def fasta_ids(path: Path) -> Iterator[str]:
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            if raw.startswith(">"):
                yield normalize_query(raw[1:].strip().split()[0])


def normalize_query(value: object) -> str:
    token = str(value).strip().split()[0]
    parts = token.split("|")
    return parts[1] if len(parts) >= 3 and parts[0].lower() in {"sp", "tr"} else token


def parse_target(value: object) -> tuple[str, str, str]:
    token = str(value).strip().split()[0]
    parts = token.split("|", maxsplit=2)
    if len(parts) != 3 or parts[0].upper() not in VALID_SOURCES or not all(parts[1:]):
        raise ValueError(
            f"Invalid target identifier {token!r}; expected SOURCE|ACCESSION|MOTIF_ID"
        )
    return parts[0].upper(), parts[1], parts[2]


def add_unique(store: dict[str, list[str]], tf: str, pwm: str) -> None:
    if pwm not in store[tf]:
        store[tf].append(pwm)


def populated(value: object) -> bool:
    return not pd.isna(value) and str(value).strip().lower() not in {
        "", "nan", "none", "na", "n/a", "-"
    }


def canonical_domains(value: object) -> set[str]:
    if not populated(value):
        return set()
    text = str(value).lower()
    aliases = {
        "c2h2": (r"c2h2",),
        "bed": (r"\bbed\b", r"zf-bed"),
        "myb_sant": (r"\bmyb\b", r"\bsant\b"),
        "cenpb": (r"cenp.?b",),
        "at_hook": (r"at.?hook",),
        "cxxc": (r"cxxc", r"cxc"),
        "homeodomain": (r"homeo",),
        "sand": (r"\bsand\b",),
        "gata": (r"\bgata\b",),
        "bhlh": (r"\bbhlh\b", r"helix.loop.helix", r"\bhlh\b"),
        "bzip": (r"\bbzip\b", r"brlz"),
        "methyl_cpg": (r"methyl.?cpg", r"\bmbd\b"),
        "hmg": (r"hmg", r"\bsox\b"),
        "ndt80": (r"ndt80", r"phog"),
        "grainyhead": (r"grainyhead", r"cp2"),
        "gtf2i": (r"gtf2i",),
        "flywch": (r"flywch",),
        "ccch": (r"ccch",),
    }
    found = {
        canonical
        for canonical, patterns in aliases.items()
        if any(re.search(pattern, text) for pattern in patterns)
    }
    if found:
        return found
    return {
        token
        for token in re.split(r"[^a-z0-9]+", text)
        if len(token) >= 3 and token not in {"unknown", "domain"}
    }


def load_domain_map(path: Path | None) -> dict[str, tuple[set[str], int | None]]:
    if path is None:
        return {}
    table = pd.read_csv(path, sep="\t", dtype=str)
    missing = [column for column in DOMAIN_COLUMNS if column not in table.columns]
    if missing:
        raise ValueError(f"Domain map is missing columns: {missing}")
    result: dict[str, tuple[set[str], int | None]] = {}
    for row in table[DOMAIN_COLUMNS].itertuples(index=False):
        identifier = str(row.identifier).strip()
        if not identifier:
            continue
        domains = canonical_domains(row.DBD_family)
        count = int(float(row.DBD_count)) if populated(row.DBD_count) else None
        current = result.get(identifier)
        if current and current != (domains, count):
            raise ValueError(f"Conflicting domain annotations for {identifier}")
        result[identifier] = (domains, count)
    return result


def domain_compatibility(
    query: str,
    target: str,
    domain_map: dict[str, tuple[set[str], int | None]],
) -> str:
    query_domains, query_count = domain_map.get(query, (set(), None))
    target_domains, target_count = domain_map.get(target, (set(), None))
    if not query_domains or not target_domains:
        return "unknown"
    shared = query_domains & target_domains
    if not shared:
        return "no"
    if "c2h2" in shared and query_count is not None and target_count is not None:
        return "yes" if query_count == target_count else "no"
    return "yes"


def load_direct_assignments(paths: list[Path]) -> pd.DataFrame:
    tables: list[pd.DataFrame] = []
    for path in paths:
        table = pd.read_csv(path, sep="\t", dtype=str)
        missing = [column for column in DIRECT_COLUMNS if column not in table.columns]
        if missing:
            raise ValueError(f"{path} is missing direct-assignment columns: {missing}")
        tables.append(table[DIRECT_COLUMNS].copy())
    if not tables:
        return pd.DataFrame(columns=DIRECT_COLUMNS)
    assignments = pd.concat(tables, ignore_index=True).fillna("")
    for column in DIRECT_COLUMNS:
        assignments[column] = assignments[column].astype(str).str.strip()
    invalid_sources = sorted(set(assignments["PWM_source"].str.upper()) - VALID_SOURCES)
    invalid_bases = sorted(set(assignments["match_basis"]) - VALID_MATCH_BASES)
    if invalid_sources:
        raise ValueError(f"Invalid direct-assignment sources: {invalid_sources}")
    if invalid_bases:
        raise ValueError(f"Invalid direct match_basis values: {invalid_bases}")
    if (assignments[DIRECT_COLUMNS] == "").any(axis=None):
        raise ValueError("Direct assignments contain blank required values")
    assignments["PWM_source"] = assignments["PWM_source"].str.upper()
    return assignments.drop_duplicates()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--queries", required=True, type=Path, help="Query TF FASTA")
    parser.add_argument("--comparison", required=True, type=Path, help="MMseqs2 .comparison file")
    parser.add_argument("--output", required=True, type=Path, help="Wide annotation TSV")
    parser.add_argument("--provenance-output", type=Path, help="Optional long-format provenance TSV")
    parser.add_argument(
        "--direct-assignments", nargs="*", type=Path, default=[],
        help="TSV files with explicit TF-to-PWM direct assignments",
    )
    parser.add_argument(
        "--domain-map", type=Path,
        help="Optional identifier-to-DBD map; known incompatible hits are rejected",
    )
    parser.add_argument("--high-threshold", type=float, default=70.0)
    parser.add_argument("--lower-threshold", type=float, default=50.0)
    parser.add_argument("--min-target-coverage", type=float, default=0.80)
    parser.add_argument("--min-alignment-length", type=int, default=40)
    parser.add_argument("--max-evalue", type=float, default=1e-5)
    parser.add_argument("--identity-tie-tolerance", type=float, default=0.05)
    parser.add_argument("--bits-tie-tolerance", type=float, default=0.05)
    args = parser.parse_args()

    if not 0 <= args.lower_threshold < args.high_threshold <= 100:
        raise ValueError("Expected 0 <= lower threshold < high threshold <= 100")
    if not 0 <= args.min_target_coverage <= 1:
        raise ValueError("--min-target-coverage must be between 0 and 1")
    if args.min_alignment_length < 1 or args.max_evalue < 0:
        raise ValueError("Alignment length must be positive and e-value non-negative")

    query_ids = sorted(dict.fromkeys(fasta_ids(args.queries)))
    query_set = set(query_ids)
    matches = pd.read_csv(
        args.comparison, sep="\t", header=None, names=MMSEQS_COLUMNS,
        dtype={"query": str, "target": str},
    )
    if matches[["qcov", "tcov", "qlen", "tlen"]].isna().all().any():
        raise ValueError(
            "Expected the 16-column output from scripts/06_run_mmseqs.sh, including coverage and lengths"
        )
    numeric = ["pident", "alnlen", "evalue", "bits", "qcov", "tcov", "qlen", "tlen"]
    for column in numeric:
        matches[column] = pd.to_numeric(matches[column], errors="raise")
    if not matches.empty and matches["pident"].max() <= 1.0:
        raise ValueError("The comparison file contains fractional rather than percentage identity")

    matches["TF_name"] = matches["query"].map(normalize_query)
    parsed = matches["target"].map(parse_target)
    matches[["PWM_source", "reference_accession", "PWM"]] = pd.DataFrame(
        parsed.tolist(), index=matches.index
    )
    domain_map = load_domain_map(args.domain_map)
    matches["domain_compatible"] = [
        domain_compatibility(query, target, domain_map)
        for query, target in zip(matches["TF_name"], matches["reference_accession"])
    ]

    direct = load_direct_assignments(args.direct_assignments)
    unknown_direct = sorted(set(direct["TF_name"]) - query_set)
    if unknown_direct:
        raise ValueError(
            f"Direct assignments contain TFs absent from the query FASTA: {unknown_direct[:5]}"
        )

    exact_self = matches.loc[
        matches["TF_name"].eq(matches["reference_accession"])
        & matches["pident"].eq(100.0)
        & matches["qcov"].ge(0.99)
        & matches["tcov"].ge(0.99),
        ["TF_name", "PWM_source", "reference_accession", "PWM"],
    ].copy()
    exact_self["match_basis"] = "stable_accession"
    direct = pd.concat([direct, exact_self[DIRECT_COLUMNS]], ignore_index=True).drop_duplicates()

    direct_store: dict[str, list[str]] = defaultdict(list)
    provenance_rows: list[dict[str, object]] = []
    for row in direct.itertuples(index=False):
        add_unique(direct_store, row.TF_name, row.PWM)
        hit = matches.loc[
            matches["TF_name"].eq(row.TF_name)
            & matches["PWM_source"].eq(row.PWM_source)
            & matches["reference_accession"].eq(row.reference_accession)
            & matches["PWM"].eq(row.PWM)
        ].sort_values(["pident", "bits"], ascending=False).head(1)
        stats = hit.iloc[0] if len(hit) else None
        provenance_rows.append(
            {
                "TF_name": row.TF_name,
                "annotation_level": "Direct_PWM",
                "PWM_source": row.PWM_source,
                "reference_accession": row.reference_accession,
                "PWM": row.PWM,
                "match_basis": row.match_basis,
                "pident": stats["pident"] if stats is not None else pd.NA,
                "alnlen": stats["alnlen"] if stats is not None else pd.NA,
                "qcov": stats["qcov"] if stats is not None else pd.NA,
                "tcov": stats["tcov"] if stats is not None else pd.NA,
                "evalue": stats["evalue"] if stats is not None else pd.NA,
                "bits": stats["bits"] if stats is not None else pd.NA,
                "domain_compatible": stats["domain_compatible"] if stats is not None else "not_applicable",
            }
        )

    direct_tfs = set(direct_store)
    qualified = matches.loc[
        ~matches["TF_name"].isin(direct_tfs)
        & matches["pident"].ge(args.lower_threshold)
        & matches["alnlen"].ge(args.min_alignment_length)
        & matches["tcov"].ge(args.min_target_coverage)
        & matches["evalue"].le(args.max_evalue)
        & matches["domain_compatible"].ne("no")
    ].copy()

    homologous: dict[str, list[str]] = defaultdict(list)
    relative: dict[str, list[str]] = defaultdict(list)
    for tf, group in qualified.groupby("TF_name", sort=False):
        best_identity = group["pident"].max()
        nearest = group.loc[
            group["pident"].ge(best_identity - args.identity_tie_tolerance)
        ].copy()
        best_bits = nearest["bits"].max()
        nearest = nearest.loc[
            nearest["bits"].ge(best_bits - args.bits_tie_tolerance)
        ].sort_values(["PWM_source", "reference_accession", "PWM"])
        level = (
            "Homologous_PWM"
            if best_identity >= args.high_threshold
            else "Relatively_Homologous_PWM"
        )
        store = homologous if level == "Homologous_PWM" else relative
        for row in nearest.itertuples(index=False):
            add_unique(store, tf, row.PWM)
            provenance_rows.append(
                {
                    "TF_name": tf,
                    "annotation_level": level,
                    "PWM_source": row.PWM_source,
                    "reference_accession": row.reference_accession,
                    "PWM": row.PWM,
                    "match_basis": "nearest_sequence",
                    "pident": row.pident,
                    "alnlen": row.alnlen,
                    "qcov": row.qcov,
                    "tcov": row.tcov,
                    "evalue": row.evalue,
                    "bits": row.bits,
                    "domain_compatible": row.domain_compatible,
                }
            )

    rows = [
        {
            "TF_name": tf,
            "TF_family": "",
            "Direct_PWM": ";".join(direct_store.get(tf, [])),
            "Homologous_PWM": ";".join(homologous.get(tf, [])),
            "Relatively_Homologous_PWM": ";".join(relative.get(tf, [])),
            "ModCRE": "",
            "AlphaFold": "",
        }
        for tf in query_ids
    ]
    pd.DataFrame(rows).to_csv(args.output, sep="\t", index=False)
    if args.provenance_output:
        provenance_columns = [
            "TF_name", "annotation_level", "PWM_source", "reference_accession",
            "PWM", "match_basis", "pident", "alnlen", "qcov", "tcov",
            "evalue", "bits", "domain_compatible",
        ]
        pd.DataFrame(provenance_rows, columns=provenance_columns).drop_duplicates().to_csv(
            args.provenance_output, sep="\t", index=False
        )
    print(f"Wrote {len(rows)} TF rows to {args.output}")
    print(
        f"Direct: {len(direct_store):,}; high homology: {len(homologous):,}; "
        f"moderate homology: {len(relative):,}"
    )


if __name__ == "__main__":
    main()
