# Hierarchical human TF–PWM annotation

A curated annotation table linking human transcription factors (TFs) to position weight matrices (PWMs) or predicted models through a fixed evidence hierarchy:

1. identical-sequence PWM;
2. high-homology PWM transfer;
3. moderate-homology PWM transfer;
4. ModCRE model;
5. AlphaFold-derived model.

The current release contains **5,417 unique TF accessions**. Of these, **5,288 (97.6%)** have one retained annotation and **129** remain unannotated.

## Current release

| Best annotation level | TFs |
|---|---:|
| Identical PWM | 2,160 |
| Homologous PWM | 1,786 |
| Relatively homologous PWM | 700 |
| ModCRE | 400 |
| AlphaFold | 242 |
| Unannotated | 129 |
| **Total** | **5,417** |

Download the main table: [`data/releases/TF_PWM_chart_final.tsv`](data/releases/TF_PWM_chart_final.tsv)

The 2025 pre-HOCOMOCO release (5,384 TFs) is retained as [`data/releases/TF_PWM_chart_pre_hocomoco.tsv`](data/releases/TF_PWM_chart_pre_hocomoco.tsv). HOCOMOCO v11 was incorporated in a subsequent 2026 update, adding 33 TFs and improving the evidence level of 29 existing TFs.

## Repository structure

```text
.
├── data/
│   ├── releases/           # versioned publication tables
│   ├── interim/            # retained homology-stage table
│   ├── reference/          # family mapping and HOCOMOCO annotation
│   ├── model_predictions/  # ModCRE/AF3 logs and exclusion lists
│   └── audit/              # HOCOMOCO update audit trail
├── scripts/                # numbered processing and validation scripts
├── docs/                   # data dictionary and reproducibility notes
└── .github/workflows/      # automated release validation
```

## Quick validation

```bash
python -m pip install -r requirements.txt
python scripts/09_validate_release.py data/releases/TF_PWM_chart_final.tsv
```

To reproduce the verified HOCOMOCO update:

```bash
python scripts/08_integrate_hocomoco.py \
  --input data/releases/TF_PWM_chart_pre_hocomoco.tsv \
  --hocomoco data/reference/HOCOMOCOv11_core_annotation_HUMAN_mono.tsv \
  --output reproduced_final.tsv \
  --audit reproduced_audit.tsv

cmp reproduced_final.tsv data/releases/TF_PWM_chart_final.tsv
cmp reproduced_audit.tsv data/audit/HOCOMOCO_integration_audit.tsv
```

## Workflow

The original 2025 workflow combined JASPAR and CisBP sequence–motif records, searched the human TF set using MMseqs2, added TF-family annotations, and then attached filtered ModCRE and AlphaFold-derived predictions. The 2026 update added HOCOMOCO v11 and enforced the hierarchy so that each TF retains only its highest-priority evidence column.

The 70% and 50% sequence-identity cutoffs are operational thresholds used in this project, not universal biological boundaries. They separate high-confidence transfer (≥70%) from a lower-confidence 50–70% tier. The motivation for using structure-based ModCRE predictions under remote homology follows [Fornes et al. (2024)](https://doi.org/10.1093/nargab/lqae068).

Full field definitions and limitations are documented in [`docs/DATA_DICTIONARY.md`](docs/DATA_DICTIONARY.md) and [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md).

## Data and software sources

- [JASPAR](https://jaspar.elixir.no/) (2024 release)
- [Cis-BP](http://cisbp.ccbr.utoronto.ca/) (version 2.00)
- [HOCOMOCO](https://hocomoco11.autosome.org/) (version 11, human core mononucleotide models)
- [MMseqs2](https://github.com/soedinglab/MMseqs2)
- [ModCRE](https://doi.org/10.1093/nargab/lqae068)

Users should also cite the original databases and methods appropriate to their analysis.

## Status

The final release and HOCOMOCO integration are exactly reproducible from the files in this repository. The historical raw MMseqs2 alignment file was not archived; consequently, the homology-search stage cannot be reproduced byte-for-byte, although its 5,384-row intermediate output is retained and validated.
