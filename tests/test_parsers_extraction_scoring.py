from __future__ import annotations

import pandas as pd

from aessp.api.crossref import parse_crossref_works
from aessp.api.ncbi import parse_protein_summary_records, parse_pubmed_xml_records
from aessp.extraction.mentions import extract_mentions_from_text
from aessp.scoring.score import score_candidates, write_score_outputs


def test_ncbi_pubmed_and_protein_parsers_work_on_mock_payloads():
    xml = """
    <PubmedArticleSet>
      <PubmedArticle>
        <MedlineCitation>
          <PMID>12345</PMID>
          <Article>
            <ArticleTitle>Dextran from Weissella confusa strain X</ArticleTitle>
            <Abstract><AbstractText>Dextransucrase produced dextran.</AbstractText></Abstract>
            <Journal><Title>Mock Journal</Title><JournalIssue><PubDate><Year>2024</Year></PubDate></JournalIssue></Journal>
            <AuthorList><Author><ForeName>Ada</ForeName><LastName>Lovelace</LastName></Author></AuthorList>
          </Article>
        </MedlineCitation>
        <PubmedData><ArticleIdList><ArticleId IdType="doi">10.123/mock</ArticleId></ArticleIdList></PubmedData>
      </PubmedArticle>
    </PubmedArticleSet>
    """
    pubmed = parse_pubmed_xml_records(xml, query="dextran", fetched_at="2026-01-01T00:00:00+00:00")

    assert pubmed[0]["pmid"] == "12345"
    assert pubmed[0]["doi"] == "10.123/mock"
    assert pubmed[0]["authors"] == "Ada Lovelace"

    protein_payload = {
        "result": {
            "uids": ["1"],
            "1": {
                "accessionversion": "ABC123.1",
                "title": "dextransucrase EC 2.4.1.5 [Leuconostoc mesenteroides]",
                "slen": 1527,
                "taxid": 1245,
            },
        }
    }
    proteins = parse_protein_summary_records(
        protein_payload, query="dextransucrase", fetched_at="2026-01-01T00:00:00+00:00"
    )

    assert proteins[0]["accession"] == "ABC123.1"
    assert proteins[0]["organism"] == "Leuconostoc mesenteroides"
    assert proteins[0]["ec_number"] == "2.4.1.5"
    assert proteins[0]["sequence_hash"] == ""


def test_crossref_parser_works_on_mock_payload():
    payload = {
        "message": {
            "items": [
                {
                    "DOI": "10.5555/dextran",
                    "title": ["Dextran branching by NMR"],
                    "container-title": ["Carbohydrate Mock Letters"],
                    "issued": {"date-parts": [[2023, 5, 1]]},
                    "author": [{"given": "Grace", "family": "Hopper"}],
                    "URL": "https://doi.org/10.5555/dextran",
                    "abstract": "<jats:p>Mock abstract.</jats:p>",
                }
            ]
        }
    }

    records = parse_crossref_works(payload, query="branching", fetched_at="2026-01-01T00:00:00+00:00")

    assert records[0]["doi"] == "10.5555/dextran"
    assert records[0]["title"] == "Dextran branching by NMR"
    assert records[0]["year"] == "2023"
    assert records[0]["authors"] == "Grace Hopper"
    assert records[0]["abstract"] == "Mock abstract."


def test_candidate_mention_extraction_flags_numeric_values_for_review():
    record = extract_mentions_from_text(
        title="Weissella confusa strain DSM 20196 dextransucrase",
        abstract=(
            "The strain produced dextran with molecular weight 2.3 x 10^6 Da. "
            "1H NMR indicated alpha-1,6 linkages with alpha-1,3 branching. "
            "Yield reached 12 g/L in a flask."
        ),
        record_id="rec1",
    )

    assert record["species_mentioned"] == "Weissella confusa"
    assert "DSM 20196" in record["strain_mentioned"]
    assert "dextransucrase" in record["enzyme_terms"]
    assert record["mw_mentions"]
    assert record["yield_mentions"]
    assert record["needs_manual_review"] is True


def test_candidate_mention_extraction_ignores_strain_stopwords():
    record = extract_mentions_from_text(
        title="Weissella confusa dextran production",
        abstract="The strain was productive and the dextran showed branching.",
        record_id="rec2",
    )

    assert "strain was" not in record["strain_mentioned"].lower()


def test_scoring_handles_missing_values():
    candidates = pd.DataFrame(
        [
            {
                "candidate_id": "cand_0001",
                "protein_accession": "ABC123.1",
                "protein_sequence_available": True,
                "literature_evidence_count": 1,
                "evidence_confidence": "",
                "needs_manual_review": True,
            }
        ]
    )
    weights = {
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

    scored = score_candidates(candidates, weights)

    assert "total_score" in scored.columns
    assert 0.0 <= scored.iloc[0]["total_score"] <= 1.0
    assert scored.iloc[0]["score_notes"]


def test_reports_do_not_contain_full_sequences(tmp_path):
    full_sequence = "M" + "A" * 80
    scored = pd.DataFrame(
        [
            {
                "candidate_id": f"cand_{index:04d}",
                "total_score": 0.9,
                "sequence": full_sequence,
                "protein_sequence_available": True,
                "manual_verified": True,
                "pilot_ready": True,
            }
            for index in range(1, 9)
        ]
    )
    out_csv = tmp_path / "scores.csv"
    top20 = tmp_path / "top20.csv"
    top8 = tmp_path / "top8.csv"

    write_score_outputs(scored, out_csv, top20, top8)

    assert full_sequence not in top20.read_text(encoding="utf-8")
    assert full_sequence not in top8.read_text(encoding="utf-8")
    assert "sequence" not in pd.read_csv(top20).columns
