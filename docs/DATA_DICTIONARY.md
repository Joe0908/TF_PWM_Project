# Data dictionary

The main release is a tab-separated table with one row per TF accession.

| Column | Meaning |
|---|---|
| `TF_name` | Human TF identifier, generally a UniProt accession. |
| `TF_family` | Imported TF-family annotation. The original 5,384-row set uses the supplied accession-to-family mapping; HOCOMOCO-only additions use HOCOMOCO's `TF family` label. These vocabularies have not been harmonised. |
| `Identical_PWM` | Semicolon-separated PWM model IDs linked by an identical sequence or direct HOCOMOCO accession match. |
| `Homologous_PWM` | PWM IDs transferred from sequence matches at the high-homology operational threshold (≥70% and <100%). |
| `Relatively_Homologous_PWM` | PWM IDs transferred from the lower-confidence operational tier (≥50% and <70%). |
| `ModCRE` | Valid retained ModCRE model names. |
| `AlphaFold` | Valid retained AlphaFold-derived model names. |
| `Best_annotation_level` | Highest non-empty evidence level after hierarchy enforcement. |
| `Best_PWM_or_model` | Value copied from the retained evidence column. |
| `N_nonempty_annotation_columns` | Number of populated evidence columns after hierarchy enforcement; expected to be 0 or 1. |

Missing values are written as empty TSV fields. Multiple models within one evidence level are separated by semicolons, with duplicates removed while preserving their first occurrence.

## Evidence hierarchy

```text
Identical_PWM
  > Homologous_PWM
  > Relatively_Homologous_PWM
  > ModCRE
  > AlphaFold
  > Unannotated
```

When a higher-priority value exists, all lower-priority evidence columns are cleared in the final release. The pre-HOCOMOCO and interim tables retain overlapping raw evidence and should not be interpreted as hierarchy-enforced releases.
