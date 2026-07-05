from __future__ import annotations

import pandas as pd

from aessp.scoring.score import score_candidates
from aessp.validation.catpred import (
    build_catpred_input_table,
    build_enzyme_validation_table,
    catpred_manifest_from_validation,
    ingest_catpred_results,
)


WEIGHTS = {
    "weights": {
        "literature_confidence_score": 0.20,
        "product_structure_score": 0.20,
        "target_mw_score": 0.15,
        "branching_evidence_score": 0.15,
        "enzyme_sequence_score": 0.10,
        "availability_score": 0.10,
        "safety_score": 0.05,
        "processability_score": 0.05,
    },
    "uncertainty": {"missing_value_penalty": 0.10, "auto_extracted_numeric_requires_review": True},
}


def _candidate(strain: str = "NRRL B-512F") -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "candidate_id": "cand_0001",
                "genus": "Leuconostoc",
                "species": "mesenteroides",
                "strain": strain,
                "organism_label": f"Leuconostoc mesenteroides {strain}".strip(),
                "enzyme_name": "dextransucrase",
                "protein_accession": "ABC123.1",
                "manual_verified": False,
            }
        ]
    )


def _protein(organism: str) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "accession": "ABC123.1",
                "protein_name": "dextransucrase",
                "gene_name": "dsr",
                "organism": organism,
                "taxonomy_id": "1245",
                "sequence_length": 1500,
                "sequence_hash": "hash123",
                "source_url": "https://www.ncbi.nlm.nih.gov/protein/ABC123.1",
                "query": "dextransucrase",
            }
        ]
    )


def test_exact_strain_matching_is_catpred_ready():
    validation = build_enzyme_validation_table(
        _candidate(),
        _protein("Leuconostoc mesenteroides NRRL B-512F"),
    )
    row = validation.iloc[0]

    assert row["accession_link_confidence"] == "exact_strain"
    assert row["enzyme_class"] == "dextransucrase"
    assert bool(row["catpred_ready"]) is True


def test_species_only_match_is_lower_confidence_but_ready():
    validation = build_enzyme_validation_table(
        _candidate(strain="NRRL B-512F"),
        _protein("Leuconostoc mesenteroides"),
    )
    row = validation.iloc[0]

    assert row["accession_link_confidence"] == "exact_species"
    assert bool(row["catpred_ready"]) is True


def test_genus_only_does_not_become_catpred_ready():
    validation = build_enzyme_validation_table(
        _candidate(),
        _protein("Leuconostoc citreum"),
    )
    row = validation.iloc[0]

    assert row["accession_link_confidence"] == "genus_only"
    assert bool(row["catpred_ready"]) is False
    assert "insufficient_accession_link_confidence" in row["catpred_block_reason"]


def test_catpred_manifest_does_not_write_full_sequences():
    full_sequence = "M" + "A" * 80
    validation = build_enzyme_validation_table(_candidate(), _protein("Leuconostoc mesenteroides NRRL B-512F"))
    validation["sequence"] = full_sequence

    manifest = catpred_manifest_from_validation(validation, only_ready=True)
    manifest_text = manifest.to_csv(index=False)

    assert "sequence" not in manifest.columns
    assert full_sequence not in manifest_text


def test_catpred_input_table_marks_missing_substrate_smiles():
    validation = build_enzyme_validation_table(_candidate(), _protein("Leuconostoc mesenteroides NRRL B-512F"))
    manifest = catpred_manifest_from_validation(validation, only_ready=True)

    table = build_catpred_input_table(
        validation,
        manifest,
        {"substrates": {"sucrose": {"smiles": "", "source_note": "manual definition required"}}},
    )

    assert table.iloc[0]["substrate_smiles"] == ""
    assert table.iloc[0]["catpred_input_status"] == "needs_substrate_definition"


def test_catpred_results_ingestion_fails_clearly_on_missing_columns():
    validation = build_enzyme_validation_table(_candidate(), _protein("Leuconostoc mesenteroides NRRL B-512F"))
    output = pd.DataFrame([{"protein_accession": "ABC123.1", "predicted_kcat": 1.0}])

    try:
        ingest_catpred_results(output, validation)
    except ValueError as exc:
        assert "CatPred output missing required columns" in str(exc)
        assert "predicted_Km" in str(exc)
    else:  # pragma: no cover - defensive failure branch
        raise AssertionError("missing CatPred columns should raise")


def test_candidate_scoring_does_not_treat_missing_catpred_as_result():
    candidates = pd.DataFrame(
        [
            {
                "candidate_id": "cand_0001",
                "protein_accession": "ABC123.1",
                "protein_sequence_available": True,
                "literature_evidence_count": 1,
                "evidence_confidence": 0.8,
            }
        ]
    )

    scored = score_candidates(candidates, WEIGHTS)
    row = scored.iloc[0]

    assert not bool(row["catpred_result_available"])
    assert row["catpred_kinetic_score"] == ""
    assert row["catpred_validation_level"] == "not_yet_validated"


def test_genus_only_link_cannot_receive_catpred_score():
    candidates = pd.DataFrame(
        [
            {
                "candidate_id": "cand_0001",
                "protein_accession": "ABC123.1",
                "protein_sequence_available": True,
                "literature_evidence_count": 1,
                "evidence_confidence": 0.8,
            }
        ]
    )
    catpred_scores = pd.DataFrame(
        [
            {
                "candidate_id": "cand_0001",
                "protein_accession": "ABC123.1",
                "accession_link_confidence": "genus_only",
                "catpred_result_available": True,
                "catpred_overall_score": 1.0,
            }
        ]
    )

    scored = score_candidates(candidates, WEIGHTS, catpred_scores=catpred_scores)
    row = scored.iloc[0]

    assert not bool(row["catpred_result_available"])
    assert row["catpred_kinetic_score"] == ""
    assert row["catpred_validation_level"] == "blocked_link_confidence"
