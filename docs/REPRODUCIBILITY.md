# Reproducibility and provenance

## Canonical construction

The workflow treats JASPAR, Cis-BP and HOCOMOCO as parallel direct-PWM resources:

1. convert each database export to the common `SOURCE|ACCESSION|MOTIF_ID` FASTA-header schema;
2. merge all sequence–PWM records without discarding distinct motif identifiers;
3. search the complete human TF query set against the combined reference with MMseqs2;
4. assign exact, ≥70% and 50–70% sequence tiers;
5. attach TF-family annotations;
6. attach valid ModCRE and AlphaFold3-assisted ModCRE models only to TFs not covered by an earlier tier;
7. retain one highest-priority evidence level per TF;
8. validate schema, uniqueness, hierarchy and release counts.

The optional long-format sequence provenance output preserves the database source, reference accession, motif identifier, percentage identity and bit score for every retained sequence-supported match.

## Input scope

The source manifest at `data/reference/SOURCES.tsv` records the resource version, selected scope and expected input form. HOCOMOCO is restricted to human CORE mononucleotide models. JASPAR, Cis-BP and HOCOMOCO are all processed before the direct-PWM reference is constructed.

## Reproducibility boundary

The final table is deposited with a SHA-256 checksum and is validated automatically. The repository contains the complete transformation and validation code, but does not redistribute every upstream database export, the comprehensive human TF query FASTA or the raw MMseqs2 alignment output. An end-to-end rerun therefore requires the source inputs listed in the manifest.

Database exports should be snapshotted locally because live resources may change. MMseqs2 version, command-line parameters and the generated `.comparison` file should be retained for any new release.

## Threshold interpretation

The 70% and 50% cutoffs are project-specific operating thresholds. They should not be interpreted as universal TF-family-independent biological boundaries. All qualifying unique motif IDs within a tier are retained; the table does not claim that every entry is a single nearest-neighbour assignment.

## Structural-model filtering

Model logs may contain repeated `DONE` records. Processing deduplicates model names before joining them into a table cell. Entries in `wrong_model_modcre.txt` and `wrong_model_af3.txt` are excluded. The final release contains none of those invalid model names.

