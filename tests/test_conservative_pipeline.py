from __future__ import annotations

import pandas as pd

from aessp.scoring.score import score_candidates, write_score_outputs
from scripts.build_candidate_table import build_candidate_table


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


def _literature() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "source_record_id": "lit1",
                "pmid": "123",
                "doi": "10.mock/dextran",
                "abstract": "Weissella confusa produced dextran.",
            }
        ]
    )


def test_no_candidate_created_from_empty_species_mention():
    mentions = pd.DataFrame(
        [
            {
                "record_id": "lit1",
                "species_mentioned": "",
                "strain_mentioned": "",
                "enzyme_terms": "dextransucrase",
                "product_terms": "dextran",
                "confidence": 0.9,
                "needs_manual_review": True,
            }
        ]
    )

    strain_candidates, enzyme_candidates = build_candidate_table(_literature(), pd.DataFrame(), mentions)

    assert strain_candidates.empty
    assert enzyme_candidates.empty


def test_auto_extraction_fills_evidence_text_not_reported_fields():
    mentions = pd.DataFrame(
        [
            {
                "record_id": "lit1",
                "species_mentioned": "Weissella confusa",
                "strain_mentioned": "strain X",
                "enzyme_terms": "dextransucrase",
                "product_terms": "dextran",
                "mw_mentions": "molecular weight 2.0 x 10^6 Da",
                "branching_mentions": "alpha-1,6 linkage mentioned",
                "nmr_mentions": "1H NMR indicated linkage",
                "viscosity_mentions": "viscosity increased",
                "yield_mentions": "yield 12 g/L",
                "confidence": 0.9,
                "needs_manual_review": True,
            }
        ]
    )

    strain_candidates, _ = build_candidate_table(_literature(), pd.DataFrame(), mentions)
    row = strain_candidates.iloc[0]

    assert row["mw_evidence_text"]
    assert row["branching_evidence_text"]
    assert row["nmr_evidence_text"]
    assert row["reported_Mw"] == ""
    assert row["reported_branching"] == ""
    assert row["reported_yield"] == ""
    assert row["reported_viscosity"] == ""


def test_mention_only_branching_does_not_get_full_score():
    candidates = pd.DataFrame(
        [
            {
                "candidate_id": "cand_0001",
                "branching_evidence_text": "branching was mentioned in the abstract",
                "availability_status": "manual_review_needed",
                "manual_verified": False,
                "needs_manual_review": True,
            }
        ]
    )

    scored = score_candidates(candidates, WEIGHTS)

    assert scored.iloc[0]["branching_evidence_score"] <= 0.3
    assert scored.iloc[0]["product_structure_score"] < 1.0


def test_unmatched_protein_goes_to_enzyme_table_not_strain_table():
    proteins = pd.DataFrame(
        [
            {
                "query": "dextransucrase",
                "accession": "ABC123.1",
                "protein_name": "dextransucrase",
                "gene_name": "dsr",
                "organism": "Leuconostoc mesenteroides",
                "taxonomy_id": "1245",
                "sequence_length": 1500,
                "sequence_hash": "",
                "ec_number": "2.4.1.5",
                "source_url": "https://www.ncbi.nlm.nih.gov/protein/ABC123.1",
            }
        ]
    )

    strain_candidates, enzyme_candidates = build_candidate_table(pd.DataFrame(), proteins, pd.DataFrame())

    assert strain_candidates.empty
    assert len(enzyme_candidates) == 1
    assert enzyme_candidates.iloc[0]["accession"] == "ABC123.1"
    assert not bool(enzyme_candidates.iloc[0]["literature_linked"])


def test_unknown_availability_score_is_conservative():
    candidates = pd.DataFrame(
        [
            {"candidate_id": "cand_unknown", "availability_status": "", "needs_manual_review": True},
            {
                "candidate_id": "cand_manual",
                "availability_status": "manual_review_needed",
                "needs_manual_review": True,
            },
        ]
    )

    scored = score_candidates(candidates, WEIGHTS).set_index("candidate_id")

    assert scored.loc["cand_unknown", "availability_score"] == 0.0
    assert scored.loc["cand_manual", "availability_score"] <= 0.1


def test_top8_is_not_real_recommendations_without_manual_ready_candidates(tmp_path):
    scored = pd.DataFrame(
        [
            {
                "candidate_id": "cand_0001",
                "total_score": 0.8,
                "manual_verified": False,
                "pilot_ready": False,
            }
        ]
    )
    out_csv = tmp_path / "scores.csv"
    top20 = tmp_path / "top20.csv"
    top8 = tmp_path / "top8_pilot_screening_recommendations.csv"

    write_score_outputs(scored, out_csv, top20, top8)

    assert out_csv.exists()
    assert top20.exists()
    assert not top8.exists()
    not_ready = tmp_path / "top8_pilot_screening_NOT_READY.md"
    assert "not generated from auto-extracted evidence alone" in not_ready.read_text(encoding="utf-8")
