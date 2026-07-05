# Methods

## Overview

AESSP is an agent-driven evidence pipeline for dextran strain and enzyme prioritization.

The workflow is intentionally split into low-complexity automation and high-complexity scientific judgment.

```text
Official APIs and manual tables
→ evidence records
→ rule-based candidate extraction
→ candidate strain/enzyme tables
→ transparent multi-objective scoring
→ human review
→ pilot-screening shortlist
```

## Data sources

Initial supported source classes:

- PubMed metadata via NCBI E-utilities.
- Crossref metadata via the Crossref REST API.
- NCBI Protein/Nucleotide metadata via E-utilities.
- Optional UniProt metadata.
- Manual culture collection review tables.

## Candidate classes

AESSP initially focuses on:

- Dextran-producing strains.
- Dextransucrase/glucansucrase enzymes.
- Related GH70 glucansucrases such as alternansucrase, mutansucrase, and reuteransucrase when relevant.

## Evidence extraction

First-pass extraction should be rule-based:

- Species and strain mentions.
- Enzyme terms.
- Molecular-weight mentions.
- Branching/linkage terms.
- NMR evidence terms.
- Viscosity/processability terms.
- Yield terms.

All extracted numeric values should be marked for manual review.

## Scoring

The first scoring model is transparent and editable. It should include:

- literature confidence
- product-structure match
- target molecular-weight fit
- branching evidence
- enzyme sequence availability
- procurement availability
- safety/processability notes
- uncertainty penalties

The score is a prioritization tool, not a final scientific conclusion.

## Expert-tool interfaces

AESSP should expose adapter layers for:

- CatPred-style enzyme kinetic priors.
- Future KETCHUP/Simscape-style kinetic modeling.
- NMR/GPC/manual-review validation.

Adapters should never fabricate outputs. If a tool has not actually been run, the adapter should output a pending status.
