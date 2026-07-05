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

## `enzyme_validation_table.csv`

Pre-CatPred enzyme validation table. This table verifies whether a candidate
accession is technically ready for CatPred input preparation. It does not imply
that CatPred has been run.

Columns:

- `enzyme_validation_id`
- `candidate_id`
- `genus`
- `species`
- `strain`
- `organism_label`
- `enzyme_name_from_candidate`
- `protein_accession`
- `protein_name`
- `protein_organism`
- `taxonomy_id`
- `sequence_length`
- `sequence_hash`
- `source_url`
- `enzyme_class`
- `product_class_hint`
- `accession_link_confidence`
- `catpred_ready`
- `catpred_block_reason`
- `needs_manual_review`
- `manual_verified`
- `notes`

## CatPred preparation and result tables

Private CatPred inputs should live under `data/private/catpred/`.
Do not commit full FASTA sequences or private CatPred input files by default.

`catpred_input_manifest.csv` contains:

- `protein_accession`
- `enzyme_validation_id`
- `candidate_id`
- `organism_label`
- `enzyme_class`
- `sequence_length`
- `sequence_hash`
- `source_url`
- `catpred_ready`
- `needs_manual_review`

`catpred_input_table.csv` contains:

- `protein_accession`
- `enzyme_validation_id`
- `sequence_id`
- `substrate_name`
- `substrate_smiles`
- `reaction_note`
- `enzyme_class`
- `candidate_id`
- `organism_label`
- `catpred_input_status`

`enzyme_catpred_scores.csv` is created only after actual CatPred output is
ingested. It includes CatPred predictions and derived conservative score fields:

- `catpred_kcat_score`
- `catpred_Km_score`
- `catpred_efficiency_score`
- `catpred_overall_score`
- `catpred_result_available`
- `catpred_validation_level`

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
- `catpred_kinetic_score`
- `catpred_result_available`
- `catpred_validation_level`
- `evidence_priority_score`
- `catpred_validated_score`
- `production_priority_score`
- `uncertainty_penalty`
- `total_score`
- `score_notes`

## Manual review tables

Manual review is expected. Suggested files:

- `data/manual/culture_collection_manual.csv`
- `data/manual/manual_curated_candidates.csv`
- `data/reports/top20_manual_review.csv`
- `data/reports/top8_pilot_screening_recommendations.csv`
