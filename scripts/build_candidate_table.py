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
        "needs_manual_review": True,
        "_confidence_values": [],
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


def _add_mention_evidence(candidate: dict[str, object], mention: dict[str, object], literature: dict[str, object]) -> None:
    candidate["literature_evidence_count"] = int(candidate.get("literature_evidence_count", 0) or 0) + 1
    candidate["_confidence_values"].append(float(mention.get("confidence", 0) or 0))
    if str(mention.get("needs_manual_review", "")).lower() in {"true", "1", "yes"}:
        candidate["needs_manual_review"] = True

    _append_field(candidate, "source_pmids", [literature.get("pmid", "")])
    _append_field(candidate, "source_dois", [literature.get("doi", "")])
    _append_field(candidate, "reported_Mw", [mention.get("mw_mentions", "")])
    _append_field(candidate, "reported_branching", [mention.get("branching_mentions", "")])
    _append_field(candidate, "reported_viscosity", [mention.get("viscosity_mentions", "")])
    _append_field(candidate, "reported_yield", [mention.get("yield_mentions", "")])

    branching = str(mention.get("branching_mentions", "")).lower()
    if "1,6" in branching or "(1,6)" in branching:
        _append_field(candidate, "reported_alpha_1_6", ["mentioned"])
    if "1,3" in branching or "(1,3)" in branching:
        _append_field(candidate, "reported_alpha_1_3", ["mentioned"])


def _matches_protein(candidate: dict[str, object], protein: dict[str, object]) -> bool:
    organism = str(protein.get("organism", "")).lower()
    genus_species = " ".join(
        part for part in [candidate.get("genus", ""), candidate.get("species", "")] if part
    ).lower()
    if genus_species and genus_species not in organism:
        return False
    enzyme = str(candidate.get("enzyme_name", "")).lower()
    protein_name = str(protein.get("protein_name", "")).lower()
    return not enzyme or enzyme in protein_name or protein_name in enzyme


def _add_protein(candidate: dict[str, object], protein: dict[str, object]) -> None:
    _append_field(candidate, "protein_accession", [protein.get("accession", "")])
    if not candidate.get("enzyme_name") and protein.get("protein_name"):
        candidate["enzyme_name"] = protein["protein_name"]
    if protein.get("accession"):
        candidate["protein_sequence_available"] = True
    if protein.get("organism") and not candidate.get("organism_label"):
        candidate["organism_label"] = protein["organism"]


def build_candidate_table(
    literature: pd.DataFrame,
    proteins: pd.DataFrame,
    mentions: pd.DataFrame,
    manual_culture: pd.DataFrame | None = None,
) -> pd.DataFrame:
    candidates: dict[tuple[str, str, str, str], dict[str, object]] = {}
    literature_by_id = _literature_lookup(literature)

    for _, mention_row in mentions.fillna("").iterrows():
        mention = mention_row.to_dict()
        species_values = _split_values(mention.get("species_mentioned", "")) or [""]
        strain_values = _split_values(mention.get("strain_mentioned", "")) or [""]
        enzyme_values = _split_values(mention.get("enzyme_terms", "")) or [""]
        literature_record = literature_by_id.get(str(mention.get("record_id", "")), {})
        for species_label in species_values:
            genus, species = _species_parts(species_label)
            for strain in strain_values:
                for enzyme in enzyme_values:
                    key = _candidate_key(genus, species, strain, enzyme)
                    candidate = candidates.setdefault(key, _default_candidate(genus, species, strain, enzyme))
                    _add_mention_evidence(candidate, mention, literature_record)

    for _, protein_row in proteins.fillna("").iterrows():
        protein = protein_row.to_dict()
        matched_keys = [key for key, candidate in candidates.items() if _matches_protein(candidate, protein)]
        if not matched_keys:
            genus, species = _species_parts(str(protein.get("organism", "")))
            enzyme = str(protein.get("protein_name", ""))
            key = _candidate_key(genus, species, "", enzyme)
            candidates.setdefault(key, _default_candidate(genus, species, "", enzyme))
            matched_keys = [key]
        for key in matched_keys:
            _add_protein(candidates[key], protein)

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
                candidate["needs_manual_review"] = True

    rows: list[dict[str, object]] = []
    for index, candidate in enumerate(candidates.values(), start=1):
        confidences = candidate.pop("_confidence_values", [])
        if confidences:
            candidate["evidence_confidence"] = round(sum(confidences) / len(confidences), 3)
        candidate["candidate_id"] = f"cand_{index:04d}"
        rows.append({column: candidate.get(column, "") for column in CANDIDATE_COLUMNS})

    return pd.DataFrame(rows, columns=CANDIDATE_COLUMNS)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build merged dextran strain/enzyme candidate table.")
    parser.add_argument("--literature-csv", default="data/processed/literature_records.csv")
    parser.add_argument("--protein-csv", default="data/processed/protein_records.csv")
    parser.add_argument("--mentions-csv", default="data/processed/candidate_mentions.csv")
    parser.add_argument("--manual-culture-csv", default="data/manual/culture_collection_manual.csv")
    parser.add_argument("--out-csv", default="data/processed/dextran_candidate_table.csv")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    literature = read_csv_or_empty(args.literature_csv)
    proteins = read_csv_or_empty(args.protein_csv)
    mentions = read_csv_or_empty(args.mentions_csv)
    manual_path = Path(args.manual_culture_csv)
    manual = read_csv_or_empty(manual_path) if manual_path.exists() else pd.DataFrame()
    candidate_table = build_candidate_table(literature, proteins, mentions, manual)
    ensure_parent(args.out_csv)
    candidate_table.to_csv(args.out_csv, index=False)


if __name__ == "__main__":
    main()
