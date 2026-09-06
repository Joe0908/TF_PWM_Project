# Reproducibility and provenance

## Canonical construction

JASPAR, Cis-BP, HOCOMOCO and Codebook enter one direct-PWM workflow:

1. convert every source to the common `SOURCE|ACCESSION|MOTIF_ID` FASTA-header schema;
2. merge all sequence–PWM records without discarding distinct motif identifiers;
3. attach explicit stable-ID and construct-level direct assignments;
4. search the complete human TF query catalogue against the combined reference with MMseqs2;
5. filter sequence transfers by reference coverage, alignment length, E-value and known DBD compatibility;
6. retain only the highest-identity/highest-bit-score neighbour and exact ties;
7. attach TF-family annotations;
8. attach valid ModCRE and AlphaFold3-assisted ModCRE models only to TFs not covered by an earlier tier;
9. retain one highest-priority evidence level per TF;
10. validate schema, uniqueness, hierarchy, source provenance and release counts.

The long-format sequence provenance output preserves source, reference accession, motif ID, assignment basis, identity, alignment length, query/reference coverage, E-value, bit score and DBD-compatibility state.

## Direct evidence

A direct PWM requires a same-TF identifier relationship. Supported `match_basis` values are `stable_accession`, `stable_gene_id`, `exact_construct` and `exact_sequence`. A 100% local alignment to a differently identified TF is not direct evidence.

For Codebook, the source identifier is the Supplementary Table 1 Ensembl gene ID. `CODEBOOK_2026_accession_map.tsv` pins one canonical UniProt accession for every one of the 177 motif-bearing TFs and compatible alternative accessions when the assayed construct is contained exactly. The map retains `Jolma2026a` PWM IDs and does not relabel literature-derived motifs as Codebook experiments.

## Homology transfer

The production defaults require:

- identity ≥50%;
- reference/construct coverage ≥80%;
- alignment length ≥40 amino acids;
- E-value ≤1e-5;
- no known DBD incompatibility.

The best percentage-identity hit is selected first, bit score breaks identity ties, and exact ties retain all associated PWM IDs. Identity alone is not used to create direct evidence. DBD compatibility is enforced when both identifiers have mapped DBD annotations; missing DBD metadata remains visible as `unknown` in provenance.

## Input scope

`data/reference/SOURCES.tsv` records resource version, selected scope and expected input form. HOCOMOCO is restricted to human CORE mononucleotide models. Codebook uses the final Nature supplementary Tables 1–4 and 14; the 177 experimental motifs are distributed under Cis-BP study accession `Jolma2026a`.

The TF census snapshot records 1,638 assessed human TF genes plus one non-TF reassessment. Its motif-status totals are validated independently, including the 146 `No motif - Codebook tested` entries. A negative assay result is an audit flag, not proof that a protein cannot bind DNA, so model evidence is not deleted solely because of that status.

## Reproducibility boundary

The release, normalized Codebook provenance, accession map, DBD map and census are deposited with SHA-256 checksums. The repository contains transformation and validation code, but does not redistribute every upstream database export, the full human TF query FASTA, the raw Codebook workbooks or MMseqs2 alignments. An end-to-end rerun requires the source inputs listed in the manifest.

Database exports should be snapshotted locally because live resources may change. Retain the MMseqs2 version, command-line parameters and generated `.comparison` file for every regenerated release.

The raw MMseqs2 table underlying the deposited snapshot is not present. Its retained homology cells can therefore be checked for nearest percentage identity, but their coverage, E-value and DBD filters cannot be reconstructed from the snapshot alone. Those filters are enforced and recorded by the canonical workflow for every full regeneration. Codebook direct assignments are independently reproducible from the deposited construct sequences and pinned UniProtKB 2026_03 crosswalk.

## Structural-model filtering

Model logs may contain repeated `DONE` records. Processing deduplicates model names before joining them into a table cell. Entries in `wrong_model_modcre.txt` and `wrong_model_af3.txt` are excluded, and validation confirms that none appear in the final release.
