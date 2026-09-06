# Data dictionary

The main release is a tab-separated table with one row per human TF accession.

| Column | Meaning |
|---|---|
| `TF_name` | Human TF identifier, generally a UniProt accession. |
| `TF_family` | TF-family or DBD annotation joined through the accession mapping. |
| `Direct_PWM` | Semicolon-separated experimental PWM IDs assigned to the same TF by stable ID, exact construct or exact sequence evidence. JASPAR, Cis-BP, HOCOMOCO and Codebook all contribute to this field. |
| `Homologous_PWM` | PWM IDs from the nearest compatible reference at the high-homology operational threshold (≥70% identity). |
| `Relatively_Homologous_PWM` | PWM IDs from the nearest compatible reference at the moderate-homology operational threshold (≥50% and <70% identity). |
| `ModCRE` | Valid retained ModCRE model names for TFs without sequence-based PWM evidence. |
| `AlphaFold` | Valid retained AlphaFold3-assisted ModCRE model names for TFs not covered by any earlier tier. |
| `Best_annotation_level` | Highest populated level in the fixed evidence hierarchy. |
| `Best_PWM_or_model` | Value copied from the retained evidence column. |
| `N_nonempty_annotation_columns` | Number of populated evidence columns after hierarchy enforcement; expected to be 0 or 1. |

Missing values are empty TSV fields. Multiple motifs within one evidence level are separated by semicolons, with duplicates removed while preserving deterministic order.

## Evidence hierarchy

```text
Direct_PWM
  > Homologous_PWM
  > Relatively_Homologous_PWM
  > ModCRE
  > AlphaFold
  > Unannotated
```

The four PWM resources differ only in source-specific parsing and identifier crosswalks before their sequence–PWM records are combined. Direct evidence is never inferred merely from a 100% local match to a different TF. When a higher-priority value exists, lower-priority evidence columns are cleared.

## Sequence provenance table

`scripts/07_build_homology_table.py` can write a long-format provenance table with one row per retained TF–motif assignment:

| Column | Meaning |
|---|---|
| `TF_name` | Query TF accession. |
| `annotation_level` | Direct, high-homology or moderate-homology tier. |
| `PWM_source` | `JASPAR`, `CISBP`, `HOCOMOCO` or `CODEBOOK`. |
| `reference_accession` | Stable accession, source TF identifier or Ensembl gene ID carried by the source record. |
| `PWM` | Source motif/model identifier. |
| `match_basis` | Direct-ID/construct basis or `nearest_sequence` for homology transfer. |
| `pident` | MMseqs2 percentage sequence identity, when applicable. |
| `alnlen` | Alignment length in amino acids. |
| `qcov` | Fraction of the query sequence covered by the alignment. |
| `tcov` | Fraction of the reference sequence or construct covered by the alignment. |
| `evalue` | MMseqs2 E-value. |
| `bits` | MMseqs2 bit score. |
| `domain_compatible` | `yes`, `no`, `unknown` or `not_applicable`. Known `no` matches are not retained. |

## Codebook reference tables

| File | Grain | Purpose |
|---|---|---|
| `CODEBOOK_2026_motifs.tsv` | 177 experimental TF–PWM records | Connects each motif to its Ensembl gene, assay, derivation method, experiment, plasmid, insert and assayed protein sequence. |
| `CODEBOOK_2026_accession_map.tsv` | 236 UniProt TF–PWM assignments | Pins the direct accession crosswalk and records whether it is supported by exact construct containment or stable gene identity. |
| `CODEBOOK_2026_tf_census.tsv` | 1,639 assessed genes | Stores the final TF/DBD reassessment and motif status used to audit the human TF query catalogue. |
| `TF_DBD_map.tsv` | identifier-level DBD records | Supplies DBD compatibility information for sequence-transfer filtering. |
