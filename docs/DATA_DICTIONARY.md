# Data dictionary

The main release is a tab-separated table with one row per human TF accession.

| Column | Meaning |
|---|---|
| `TF_name` | Human TF identifier, generally a UniProt accession. |
| `TF_family` | TF-family annotation joined through the supplied accession mapping. |
| `Identical_PWM` | Semicolon-separated PWM IDs obtained by an exact sequence match in the combined JASPAR, Cis-BP and HOCOMOCO direct-PWM reference. |
| `Homologous_PWM` | PWM IDs transferred from sequence matches at the high-homology operational threshold (≥70% and <100%). |
| `Relatively_Homologous_PWM` | PWM IDs transferred from the lower-confidence operational tier (≥50% and <70%). |
| `ModCRE` | Valid retained ModCRE model names for TFs without sequence-based PWM evidence. |
| `AlphaFold` | Valid retained AlphaFold3-assisted ModCRE model names for TFs not covered by any earlier tier. |
| `Best_annotation_level` | Highest populated level in the fixed evidence hierarchy. |
| `Best_PWM_or_model` | Value copied from the retained evidence column. |
| `N_nonempty_annotation_columns` | Number of populated evidence columns after hierarchy enforcement; expected to be 0 or 1. |

Missing values are written as empty TSV fields. Multiple models within one evidence level are separated by semicolons, with duplicates removed while preserving first occurrence.

## Evidence hierarchy

```text
Identical_PWM
  > Homologous_PWM
  > Relatively_Homologous_PWM
  > ModCRE
  > AlphaFold
  > Unannotated
```

JASPAR, Cis-BP and HOCOMOCO differ only in source-specific parsing before their sequence–PWM records are merged. All three contribute to the same `Identical_PWM` layer and to the same PWM-bearing reference used for homology transfer. When a higher-priority value exists, lower-priority evidence columns are cleared in the final release.

## Sequence provenance table

`scripts/06_build_homology_table.py` can optionally write a long-format provenance table with one row per retained MMseqs2-supported TF–motif match:

| Column | Meaning |
|---|---|
| `TF_name` | Query TF accession. |
| `annotation_level` | Exact, high-homology or moderate-homology tier. |
| `PWM_source` | `JASPAR`, `CISBP` or `HOCOMOCO`. |
| `reference_accession` | Accession or TF identifier carried by the source record. |
| `PWM` | Source motif/model identifier. |
| `pident` | MMseqs2 percentage sequence identity. |
| `bits` | MMseqs2 alignment bit score. |

