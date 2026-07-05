# Open Science Strategy for AESSP

AESSP should be open-science-native without being careless about licensing, provenance, or raw biological sequence data.

## Why open science helps this project

Open science turns AESSP from a private literature scrape into a reproducible research artifact. The core value is not that every raw file is public. The value is that every conclusion can be traced to an evidence source, every automatic extraction is flagged for review, and every ranking score can be regenerated.

## What should be open

- Query files and collection scripts.
- Data schemas.
- Small mock fixtures for tests.
- Processed candidate tables after manual review, when licensing allows.
- Ranking formulas and scoring weights.
- Notebooks that generate summary tables and figures.
- Limitations and manual-review notes.

## What should not be open by default

- Large raw API dumps.
- Login-protected or paywalled content.
- Supplier quotations or private culture-collection terms.
- Full protein FASTA exports unless reviewed for license and necessity.
- Company production data, NMR files, GPC files, or batch records.

## Recommended open-science artifacts

1. A versioned GitHub repository.
2. A small public demo dataset with mock records.
3. A reproducible pipeline script.
4. A candidate-table schema.
5. A manually reviewed Top 20 table with source identifiers.
6. A project report or preprint describing the workflow.
7. Optional Zenodo/OSF archive after the pipeline stabilizes.

## Evidence standard

Every candidate strain or enzyme should have at least one source field:

- PMID
- DOI
- UniProt accession
- NCBI Protein accession
- NCBI Taxonomy ID
- culture collection accession
- manual review note

Automatically extracted numeric claims must be marked `needs_manual_review = true` until checked.

## Project principle

Agents accelerate evidence gathering; they do not replace scientific judgment. AESSP rankings should always be treated as screening priorities, not biological facts.
