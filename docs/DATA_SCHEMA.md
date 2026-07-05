# Data Schema

## `literature_records.csv`

One row per literature record returned by PubMed, Crossref, or another bibliographic source.

Columns:

- `query`
- `source_api`
- `title`
- `abstract`
- `authors`
- `year`
- `journal`
- `doi`
- `pmid`
- `url`
- `fetched_at`
- `source_record_id`

## `protein_records.csv`

One row per protein record returned by NCBI Protein, UniProt, or another structured protein source.

Columns:

- `query`
- `source_api`
- `accession`
- `protein_name`
- `gene_name`
- `organism`
- `taxonomy_id`
- `sequence_length`
- `sequence_hash`
- `ec_number`
- `reviewed_status`
- `source_url`
- `fetched_at`
- `notes`

Reports should not print full sequences.

## `candidate_mentions.csv`

Rule-based candidate mentions extracted from literature metadata.

Columns:

- `record_id`
- `title`
- `species_mentioned`
- `strain_mentioned`
- `enzyme_terms`
- `product_terms`
- `mw_mentions`
- `branching_mentions`
- `nmr_mentions`
- `viscosity_mentions`
- `yield_mentions`
- `confidence`
- `needs_manual_review`

## `dextran_candidate_table.csv`

Merged strain candidate table. Automatic title/abstract extraction is stored
as evidence text and must not be treated as a verified reported value.

Columns:

- `candidate_id`
- `genus`
- `species`
- `strain`
- `organism_label`
- `enzyme_name`
- `protein_accession`
- `culture_collection`
- `culture_accession`
- `available_in_china`
- `availability_status`
- `mw_evidence_text`
- `branching_evidence_text`
- `viscosity_evidence_text`
- `yield_evidence_text`
- `nmr_evidence_text`
- `reported_yield`
- `reported_Mw`
- `reported_PDI`
- `reported_branching`
- `reported_alpha_1_6`
- `reported_alpha_1_3`
- `reported_viscosity`
- `literature_evidence_count`
- `protein_sequence_available`
- `safety_notes`
- `source_pmids`
- `source_dois`
- `evidence_confidence`
- `evidence_level`
- `structure_evidence_level`
- `enzyme_evidence_level`
- `availability_evidence_level`
- `manual_verified`
- `pilot_ready`
- `needs_manual_review`

## `enzyme_candidate_table.csv`

Protein-accession candidates that are not automatically promoted to strain
candidates unless linked to literature/manual strain evidence.

Columns:

- `enzyme_candidate_id`
- `accession`
- `protein_name`
- `gene_name`
- `organism`
- `taxonomy_id`
- `sequence_length`
- `sequence_hash`
- `ec_number`
- `source_url`
- `query`
- `literature_linked`
- `manual_relevant`
- `needs_manual_review`

## `dextran_candidate_scores.csv`

Scored candidate table.

Columns:

- all candidate-table columns
- `literature_confidence_score`
- `product_structure_score`
- `target_mw_score`
- `branching_evidence_score`
- `enzyme_sequence_score`
- `availability_score`
- `safety_score`
- `processability_score`
- `uncertainty_penalty`
- `total_score`
- `score_notes`

## Manual review tables

Manual review is expected. Suggested files:

- `data/manual/culture_collection_manual.csv`
- `data/manual/manual_curated_candidates.csv`
- `data/reports/top20_manual_review.csv`
- `data/reports/top8_pilot_screening_recommendations.csv`
