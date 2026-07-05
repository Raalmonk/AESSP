#!/usr/bin/env python
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Iterable

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from aessp.io import append_unique, ensure_parent, read_csv_or_empty


RELEVANT_PRODUCT_TERMS = {"dextran", "glucan", "eps", "polysaccharide"}
RELEVANT_ENZYME_TERMS = {
    "dextransucrase",
    "glucansucrase",
    "alternansucrase",
    "mutansucrase",
    "reuteransucrase",
    "gh70",
}


CANDIDATE_COLUMNS = [
    "candidate_id",
    "genus",
    "species",
    "strain",
    "organism_label",
    "enzyme_name",
    "protein_accession",
    "culture_collection",
    "culture_accession",
    "available_in_china",
    "availability_status",
    "mw_evidence_text",
    "branching_evidence_text",
    "viscosity_evidence_text",
    "yield_evidence_text",
    "nmr_evidence_text",
    "reported_yield",
    "reported_Mw",
    "reported_PDI",
    "reported_branching",
    "reported_alpha_1_6",
    "reported_alpha_1_3",
    "reported_viscosity",
    "literature_evidence_count",
    "protein_sequence_available",
    "safety_notes",
    "source_pmids",
    "source_dois",
    "evidence_confidence",
    "evidence_level",
    "structure_evidence_level",
    "enzyme_evidence_level",
    "availability_evidence_level",
    "manual_verified",
    "pilot_ready",
    "needs_manual_review",
]


ENZYME_CANDIDATE_COLUMNS = [
    "enzyme_candidate_id",
    "accession",
    "protein_name",
    "gene_name",
    "organism",
    "taxonomy_id",
    "sequence_length",
    "sequence_hash",
    "ec_number",
    "source_url",
    "query",
    "literature_linked",
    "manual_relevant",
    "needs_manual_review",
]


def _split_values(value: object) -> list[str]:
    return [part.strip() for part in str(value or "").split(";") if part.strip()]


def _species_parts(label: str) -> tuple[str, str]:
    parts = label.split()
    if len(parts) >= 2:
        return parts[0], parts[1]
    if len(parts) == 1:
        return parts[0], ""
    return "", ""


def _organism_label(genus: str, species: str, strain: str = "") -> str:
    base = " ".join(part for part in [genus, species] if part).strip()
    return " ".join(part for part in [base, strain] if part).strip()


def _default_candidate(genus: str, species: str, strain: str, enzyme: str) -> dict[str, object]:
    return {
        "candidate_id": "",
        "genus": genus,
        "species": species,
        "strain": strain,
        "organism_label": _organism_label(genus, species, strain),
        "enzyme_name": enzyme,
        "protein_accession": "",
        "culture_collection": "",
        "culture_accession": "",
        "available_in_china": "",
        "availability_status": "manual_review_needed",
        "mw_evidence_text": "",
        "branching_evidence_text": "",
        "viscosity_evidence_text": "",
        "yield_evidence_text": "",
        "nmr_evidence_text": "",
        "reported_yield": "",
        "reported_Mw": "",
        "reported_PDI": "",
        "reported_branching": "",
        "reported_alpha_1_6": "",
        "reported_alpha_1_3": "",
        "reported_viscosity": "",
        "literature_evidence_count": 0,
        "protein_sequence_available": False,
        "safety_notes": "",
        "source_pmids": "",
        "source_dois": "",
        "evidence_confidence": 0.0,
        "evidence_level": "mention_only",
        "structure_evidence_level": "mention_only",
        "enzyme_evidence_level": "mention_only",
        "availability_evidence_level": "mention_only",
        "manual_verified": False,
        "pilot_ready": False,
        "needs_manual_review": True,
        "_confidence_values": [],
        "_has_product_literature": True,
    }


def _candidate_key(genus: str, species: str, strain: str, enzyme: str) -> tuple[str, str, str, str]:
    return (genus.lower(), species.lower(), strain.lower(), enzyme.lower())


def _literature_lookup(literature: pd.DataFrame) -> dict[str, dict[str, object]]:
    lookup: dict[str, dict[str, object]] = {}
    for index, row in literature.fillna("").iterrows():
        row_dict = row.to_dict()
        keys = {
            str(index),
            str(row_dict.get("source_record_id", "")),
            str(row_dict.get("pmid", "")),
            str(row_dict.get("doi", "")),
        }
        for key in keys:
            if key:
                lookup[key] = row_dict
    return lookup


def _append_field(candidate: dict[str, object], field: str, values: Iterable[object]) -> None:
    candidate[field] = append_unique([candidate.get(field, ""), *values])


def _as_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "y", "verified", "ready"}


def _max_evidence_level(*levels: object) -> str:
    rank = {
        "": 0,
        "mention_only": 1,
        "abstract_supported": 2,
        "accession_supported": 3,
        "manually_verified": 4,
    }
    best = "mention_only"
    for level in levels:
        level_text = str(level or "").strip()
        if rank.get(level_text, 0) > rank[best]:
            best = level_text
    return best


def _has_relevant_product(mention: dict[str, object]) -> bool:
    product_terms = {term.lower() for term in _split_values(mention.get("product_terms", ""))}
    return bool(product_terms & RELEVANT_PRODUCT_TERMS)


def _has_relevant_enzyme_text(*values: object) -> bool:
    text = " ".join(str(value or "") for value in values).lower()
    return any(term in text for term in RELEVANT_ENZYME_TERMS)


def _can_create_strain_candidate(mention: dict[str, object]) -> bool:
    if not _split_values(mention.get("species_mentioned", "")):
        return False
    if not _has_relevant_product(mention):
        return False
    try:
        confidence = float(mention.get("confidence", 0) or 0)
    except (TypeError, ValueError):
        confidence = 0.0
    return confidence >= 0.5


def _add_mention_evidence(candidate: dict[str, object], mention: dict[str, object], literature: dict[str, object]) -> None:
    candidate["literature_evidence_count"] = int(candidate.get("literature_evidence_count", 0) or 0) + 1
    try:
        confidence = float(mention.get("confidence", 0) or 0)
    except (TypeError, ValueError):
        confidence = 0.0
    candidate["_confidence_values"].append(confidence)
    if str(mention.get("needs_manual_review", "")).lower() in {"true", "1", "yes"}:
        candidate["needs_manual_review"] = True

    _append_field(candidate, "source_pmids", [literature.get("pmid", "")])
    _append_field(candidate, "source_dois", [literature.get("doi", "")])
    _append_field(candidate, "mw_evidence_text", [mention.get("mw_mentions", "")])
    _append_field(candidate, "branching_evidence_text", [mention.get("branching_mentions", "")])
    _append_field(candidate, "viscosity_evidence_text", [mention.get("viscosity_mentions", "")])
    _append_field(candidate, "yield_evidence_text", [mention.get("yield_mentions", "")])
    _append_field(candidate, "nmr_evidence_text", [mention.get("nmr_mentions", "")])

    mention_level = "abstract_supported" if str(literature.get("abstract", "")).strip() else "mention_only"
    candidate["evidence_level"] = _max_evidence_level(candidate.get("evidence_level"), mention_level)
    if any(
        str(mention.get(field, "")).strip()
        for field in ("mw_mentions", "branching_mentions", "nmr_mentions", "viscosity_mentions", "yield_mentions")
    ):
        candidate["structure_evidence_level"] = _max_evidence_level(
            candidate.get("structure_evidence_level"), mention_level
        )


def _matches_protein(candidate: dict[str, object], protein: dict[str, object]) -> bool:
    organism = str(protein.get("organism", "")).lower()
    genus_species = " ".join(
        part for part in [candidate.get("genus", ""), candidate.get("species", "")] if part
    ).lower()
    if not genus_species or genus_species not in organism:
        return False
    if not candidate.get("_has_product_literature"):
        return False
    enzyme = str(candidate.get("enzyme_name", "")).lower()
    protein_name = str(protein.get("protein_name", "")).lower()
    query = str(protein.get("query", "")).lower()
    return _has_relevant_enzyme_text(enzyme, protein_name, query)


def _add_protein(candidate: dict[str, object], protein: dict[str, object]) -> None:
    _append_field(candidate, "protein_accession", [protein.get("accession", "")])
    if not candidate.get("enzyme_name") and protein.get("protein_name"):
        candidate["enzyme_name"] = protein["protein_name"]
    if protein.get("accession"):
        candidate["protein_sequence_available"] = True
        candidate["enzyme_evidence_level"] = _max_evidence_level(
            candidate.get("enzyme_evidence_level"), "accession_supported"
        )
        candidate["evidence_level"] = _max_evidence_level(candidate.get("evidence_level"), "accession_supported")
    if protein.get("organism") and not candidate.get("organism_label"):
        candidate["organism_label"] = protein["organism"]


def _protein_to_enzyme_candidate(
    protein: dict[str, object],
    index: int,
    literature_linked: bool = False,
) -> dict[str, object]:
    return {
        "enzyme_candidate_id": f"enz_{index:04d}",
        "accession": protein.get("accession", ""),
        "protein_name": protein.get("protein_name", ""),
        "gene_name": protein.get("gene_name", ""),
        "organism": protein.get("organism", ""),
        "taxonomy_id": protein.get("taxonomy_id", ""),
        "sequence_length": protein.get("sequence_length", ""),
        "sequence_hash": protein.get("sequence_hash", ""),
        "ec_number": protein.get("ec_number", ""),
        "source_url": protein.get("source_url", ""),
        "query": protein.get("query", ""),
        "literature_linked": literature_linked,
        "manual_relevant": False,
        "needs_manual_review": True,
    }


def build_enzyme_candidate_table(
    proteins: pd.DataFrame,
    linked_accessions: set[str] | None = None,
) -> pd.DataFrame:
    linked_accessions = linked_accessions or set()
    rows = []
    for index, (_, protein_row) in enumerate(proteins.fillna("").iterrows(), start=1):
        protein = protein_row.to_dict()
        accession = str(protein.get("accession", ""))
        rows.append(
            _protein_to_enzyme_candidate(
                protein,
                index,
                literature_linked=bool(accession and accession in linked_accessions),
            )
        )
    return pd.DataFrame(rows, columns=ENZYME_CANDIDATE_COLUMNS)


def build_candidate_table(
    literature: pd.DataFrame,
    proteins: pd.DataFrame,
    mentions: pd.DataFrame,
    manual_culture: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    candidates: dict[tuple[str, str, str, str], dict[str, object]] = {}
    literature_by_id = _literature_lookup(literature)

    for _, mention_row in mentions.fillna("").iterrows():
        mention = mention_row.to_dict()
        if not _can_create_strain_candidate(mention):
            continue
        species_values = _split_values(mention.get("species_mentioned", ""))
        strain_values = _split_values(mention.get("strain_mentioned", "")) or [""]
        enzyme_values = _split_values(mention.get("enzyme_terms", "")) or [""]
        literature_record = literature_by_id.get(str(mention.get("record_id", "")), {})
        for species_label in species_values:
            genus, species = _species_parts(species_label)
            if not genus or not species:
                continue
            for strain in strain_values:
                for enzyme in enzyme_values:
                    key = _candidate_key(genus, species, strain, enzyme)
                    candidate = candidates.setdefault(key, _default_candidate(genus, species, strain, enzyme))
                    _add_mention_evidence(candidate, mention, literature_record)

    linked_accessions: set[str] = set()
    for _, protein_row in proteins.fillna("").iterrows():
        protein = protein_row.to_dict()
        matched_keys = [key for key, candidate in candidates.items() if _matches_protein(candidate, protein)]
        for key in matched_keys:
            _add_protein(candidates[key], protein)
            accession = str(protein.get("accession", ""))
            if accession:
                linked_accessions.add(accession)

    manual_culture = manual_culture if manual_culture is not None else pd.DataFrame()
    for _, manual_row in manual_culture.fillna("").iterrows():
        manual = manual_row.to_dict()
        manual_species = str(manual.get("species", ""))
        manual_label = str(manual.get("organism_label", ""))
        manual_strain = str(manual.get("strain", ""))
        for candidate in candidates.values():
            label_match = manual_label and manual_label.lower() == str(candidate.get("organism_label", "")).lower()
            species_match = manual_species and manual_species.lower() == str(candidate.get("species", "")).lower()
            strain_match = not manual_strain or manual_strain.lower() == str(candidate.get("strain", "")).lower()
            if label_match or (species_match and strain_match):
                for column in CANDIDATE_COLUMNS:
                    if column in manual and manual[column] != "":
                        candidate[column] = manual[column]
                if _as_bool(manual.get("manual_verified", "")):
                    candidate["manual_verified"] = True
                    candidate["evidence_level"] = "manually_verified"
                    candidate["structure_evidence_level"] = "manually_verified"
                    candidate["enzyme_evidence_level"] = _max_evidence_level(
                        candidate.get("enzyme_evidence_level"), "manually_verified"
                    )
                    candidate["availability_evidence_level"] = _max_evidence_level(
                        candidate.get("availability_evidence_level"), "manually_verified"
                    )
                candidate["needs_manual_review"] = True

    rows: list[dict[str, object]] = []
    for index, candidate in enumerate(candidates.values(), start=1):
        confidences = candidate.pop("_confidence_values", [])
        candidate.pop("_has_product_literature", None)
        if confidences:
            candidate["evidence_confidence"] = round(sum(confidences) / len(confidences), 3)
        candidate["candidate_id"] = f"cand_{index:04d}"
        rows.append({column: candidate.get(column, "") for column in CANDIDATE_COLUMNS})

    strain_candidates = pd.DataFrame(rows, columns=CANDIDATE_COLUMNS)
    enzyme_candidates = build_enzyme_candidate_table(proteins, linked_accessions=linked_accessions)
    return strain_candidates, enzyme_candidates


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build merged dextran strain/enzyme candidate table.")
    parser.add_argument("--literature-csv", default="data/processed/literature_records.csv")
    parser.add_argument("--protein-csv", default="data/processed/protein_records.csv")
    parser.add_argument("--mentions-csv", default="data/processed/candidate_mentions.csv")
    parser.add_argument("--manual-culture-csv", default="data/manual/culture_collection_manual.csv")
    parser.add_argument("--out-csv", default="data/processed/dextran_candidate_table.csv")
    parser.add_argument("--enzyme-out-csv", default="data/processed/enzyme_candidate_table.csv")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    literature = read_csv_or_empty(args.literature_csv)
    proteins = read_csv_or_empty(args.protein_csv)
    mentions = read_csv_or_empty(args.mentions_csv)
    manual_path = Path(args.manual_culture_csv)
    manual = read_csv_or_empty(manual_path) if manual_path.exists() else pd.DataFrame()
    candidate_table, enzyme_candidate_table = build_candidate_table(literature, proteins, mentions, manual)
    ensure_parent(args.out_csv)
    candidate_table.to_csv(args.out_csv, index=False)
    ensure_parent(args.enzyme_out_csv)
    enzyme_candidate_table.to_csv(args.enzyme_out_csv, index=False)


if __name__ == "__main__":
    main()
