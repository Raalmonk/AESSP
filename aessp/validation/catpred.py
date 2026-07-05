from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Mapping

import pandas as pd


ENZYME_VALIDATION_COLUMNS = [
    "enzyme_validation_id",
    "candidate_id",
    "genus",
    "species",
    "strain",
    "organism_label",
    "enzyme_name_from_candidate",
    "protein_accession",
    "protein_name",
    "protein_organism",
    "taxonomy_id",
    "sequence_length",
    "sequence_hash",
    "source_url",
    "enzyme_class",
    "product_class_hint",
    "accession_link_confidence",
    "catpred_ready",
    "catpred_block_reason",
    "needs_manual_review",
    "manual_verified",
    "notes",
]

CATPRED_MANIFEST_COLUMNS = [
    "protein_accession",
    "enzyme_validation_id",
    "candidate_id",
    "organism_label",
    "enzyme_class",
    "sequence_length",
    "sequence_hash",
    "source_url",
    "catpred_ready",
    "needs_manual_review",
]

CATPRED_INPUT_COLUMNS = [
    "protein_accession",
    "enzyme_validation_id",
    "sequence_id",
    "substrate_name",
    "substrate_smiles",
    "reaction_note",
    "enzyme_class",
    "candidate_id",
    "organism_label",
    "catpred_input_status",
]

READY_ENZYME_CLASSES = {"dextransucrase", "glucansucrase", "other_GH70"}
READY_LINK_CONFIDENCE = {"exact_strain", "exact_species"}


def split_multi(value: object) -> list[str]:
    return [part.strip() for part in str(value or "").split(";") if part.strip()]


def has_value(value: object) -> bool:
    if value is None:
        return False
    if isinstance(value, float) and math.isnan(value):
        return False
    text = str(value).strip()
    return text.lower() not in {"", "nan", "none", "null", "unknown", "not_available_ncbi"}


def as_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "y", "verified", "ready", "present"}


def normalize_accession(value: object) -> str:
    return str(value or "").strip()


def accession_keys(accession: object) -> set[str]:
    text = normalize_accession(accession)
    if not text:
        return set()
    keys = {text}
    if "." in text:
        keys.add(text.split(".", 1)[0])
    return keys


def normalize_text(value: object) -> str:
    return re.sub(r"\s+", " ", str(value or "")).strip()


def normalize_match_text(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", " ", str(value or "").lower()).strip()


def classify_enzyme(*values: object) -> str:
    text = " ".join(normalize_match_text(value) for value in values)
    if "reuteransucrase" in text:
        return "reuteransucrase"
    if "alternansucrase" in text:
        return "alternansucrase"
    if "mutansucrase" in text:
        return "mutansucrase"
    if "dextransucrase" in text:
        return "dextransucrase"
    if "glucansucrase" in text:
        return "glucansucrase"
    if re.search(r"\bgh70\b", text):
        return "other_GH70"
    return "unknown"


def product_class_hint(enzyme_class: str, *values: object) -> str:
    text = " ".join(normalize_match_text(value) for value in values)
    if enzyme_class == "dextransucrase" or "dextran" in text:
        return "dextran"
    if enzyme_class == "reuteransucrase" or "reuteran" in text:
        return "reuteran"
    if enzyme_class == "alternansucrase" or "alternan" in text:
        return "alternan"
    if enzyme_class == "mutansucrase" or "mutan" in text:
        return "mutan"
    if "insoluble glucan" in text:
        return "insoluble_glucan"
    if "eps" in text or "exopolysaccharide" in text or "glucan" in text:
        return "EPS_unknown"
    return "unknown"


def _strain_token(strain: object) -> str:
    text = normalize_match_text(strain)
    text = re.sub(r"^strain\s+", "", text).strip()
    return text


def accession_link_confidence(candidate: Mapping[str, object], protein: Mapping[str, object]) -> str:
    genus = normalize_match_text(candidate.get("genus", ""))
    species = normalize_match_text(candidate.get("species", ""))
    strain = _strain_token(candidate.get("strain", ""))
    organism = normalize_match_text(protein.get("organism", ""))
    protein_text = normalize_match_text(
        " ".join(
            str(protein.get(field, ""))
            for field in ("organism", "protein_name", "gene_name", "source_url", "accession")
        )
    )

    genus_species = " ".join(part for part in [genus, species] if part)
    if strain and strain in protein_text:
        return "exact_strain"
    if genus_species and genus_species in organism:
        return "exact_species"
    if genus and re.search(rf"\b{re.escape(genus)}\b", organism):
        return "genus_only"

    candidate_enzyme = normalize_match_text(candidate.get("enzyme_name", ""))
    protein_name = normalize_match_text(protein.get("protein_name", ""))
    if candidate_enzyme and protein_name and (candidate_enzyme in protein_name or protein_name in candidate_enzyme):
        return "weak_text_match"
    return "unlinked"


def protein_lookup(proteins: pd.DataFrame) -> dict[str, dict[str, object]]:
    lookup: dict[str, dict[str, object]] = {}
    for _, row in proteins.fillna("").iterrows():
        record = row.to_dict()
        accession = record.get("accession", "")
        for key in accession_keys(accession):
            lookup[key] = record
    return lookup


def _catpred_block_reason(row: Mapping[str, object], allow_fasta_export: bool = False) -> str:
    reasons: list[str] = []
    if not has_value(row.get("protein_accession", "")):
        reasons.append("missing_protein_accession")
    if not has_value(row.get("sequence_hash", "")) and not allow_fasta_export:
        reasons.append("missing_sequence_hash_or_fasta_export")
    if row.get("enzyme_class") not in READY_ENZYME_CLASSES:
        reasons.append("unsupported_enzyme_class")
    if row.get("accession_link_confidence") not in READY_LINK_CONFIDENCE:
        reasons.append("insufficient_accession_link_confidence")
    return "; ".join(reasons)


def _manual_lookup(manual: pd.DataFrame) -> list[dict[str, object]]:
    if manual is None or manual.empty:
        return []
    return [row.to_dict() for _, row in manual.fillna("").iterrows()]


def _manual_matches(row: Mapping[str, object], manual: Mapping[str, object]) -> bool:
    manual_candidate = str(manual.get("candidate_id", "")).strip()
    manual_accession = str(manual.get("protein_accession", "") or manual.get("accession", "")).strip()
    if manual_candidate and manual_accession:
        return (
            manual_candidate == str(row.get("candidate_id", "")).strip()
            and manual_accession == str(row.get("protein_accession", "")).strip()
        )
    if manual_candidate:
        return manual_candidate == str(row.get("candidate_id", "")).strip()
    if manual_accession:
        return manual_accession == str(row.get("protein_accession", "")).strip()
    return False


def _apply_manual_overrides(row: dict[str, object], manual_rows: list[dict[str, object]]) -> dict[str, object]:
    for manual in manual_rows:
        if not _manual_matches(row, manual):
            continue
        for column in ENZYME_VALIDATION_COLUMNS:
            if column in manual and has_value(manual[column]):
                row[column] = manual[column]
        if as_bool(manual.get("manual_verified", "")):
            row["manual_verified"] = True
            row["needs_manual_review"] = False
    return row


def build_enzyme_validation_table(
    candidates: pd.DataFrame,
    proteins: pd.DataFrame,
    manual_validation: pd.DataFrame | None = None,
    allow_fasta_export: bool = False,
) -> pd.DataFrame:
    proteins_by_accession = protein_lookup(proteins)
    manual_rows = _manual_lookup(manual_validation if manual_validation is not None else pd.DataFrame())
    rows: list[dict[str, object]] = []

    for _, candidate_row in candidates.fillna("").iterrows():
        candidate = candidate_row.to_dict()
        accessions = split_multi(candidate.get("protein_accession", "")) or [""]
        for accession in accessions:
            protein = {}
            for key in accession_keys(accession):
                if key in proteins_by_accession:
                    protein = proteins_by_accession[key]
                    break

            enzyme_class = classify_enzyme(
                candidate.get("enzyme_name", ""),
                protein.get("protein_name", ""),
                protein.get("query", ""),
            )
            link_confidence = accession_link_confidence(candidate, protein) if protein else "unlinked"
            row = {
                "enzyme_validation_id": f"ev_{len(rows) + 1:04d}",
                "candidate_id": candidate.get("candidate_id", ""),
                "genus": candidate.get("genus", ""),
                "species": candidate.get("species", ""),
                "strain": candidate.get("strain", ""),
                "organism_label": candidate.get("organism_label", ""),
                "enzyme_name_from_candidate": candidate.get("enzyme_name", ""),
                "protein_accession": accession,
                "protein_name": protein.get("protein_name", ""),
                "protein_organism": protein.get("organism", ""),
                "taxonomy_id": protein.get("taxonomy_id", ""),
                "sequence_length": protein.get("sequence_length", ""),
                "sequence_hash": protein.get("sequence_hash", ""),
                "source_url": protein.get("source_url", ""),
                "enzyme_class": enzyme_class,
                "product_class_hint": product_class_hint(
                    enzyme_class,
                    candidate.get("product_terms", ""),
                    candidate.get("organism_label", ""),
                    protein.get("protein_name", ""),
                ),
                "accession_link_confidence": link_confidence,
                "catpred_ready": False,
                "catpred_block_reason": "",
                "needs_manual_review": True,
                "manual_verified": as_bool(candidate.get("manual_verified", "")),
                "notes": "",
            }
            row = _apply_manual_overrides(row, manual_rows)
            block_reason = _catpred_block_reason(row, allow_fasta_export=allow_fasta_export)
            row["catpred_block_reason"] = block_reason
            row["catpred_ready"] = not block_reason
            if row["manual_verified"]:
                row["needs_manual_review"] = False
            else:
                row["needs_manual_review"] = True
            if not protein and accession:
                row["notes"] = "protein accession was present on candidate but not found in protein_records"
            elif not accession:
                row["notes"] = "candidate has no linked protein accession"
            rows.append(row)

    return pd.DataFrame(rows, columns=ENZYME_VALIDATION_COLUMNS)


def catpred_manifest_from_validation(validation: pd.DataFrame, only_ready: bool = False) -> pd.DataFrame:
    table = validation.fillna("")
    if only_ready and "catpred_ready" in table.columns:
        table = table.loc[table["catpred_ready"].apply(as_bool)]
    rows = []
    for _, row in table.iterrows():
        rows.append({column: row.get(column, "") for column in CATPRED_MANIFEST_COLUMNS})
    return pd.DataFrame(rows, columns=CATPRED_MANIFEST_COLUMNS)


def sequence_id_for(row: Mapping[str, object]) -> str:
    base = "_".join(
        part
        for part in [
            str(row.get("enzyme_validation_id", "")).strip(),
            str(row.get("protein_accession", "")).strip().replace(".", "_"),
        ]
        if part
    )
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", base) or "catpred_sequence"


def build_catpred_input_table(
    validation: pd.DataFrame,
    manifest: pd.DataFrame,
    substrates_config: Mapping[str, object],
    substrate_name: str = "sucrose",
) -> pd.DataFrame:
    substrate_config = (
        substrates_config.get("substrates", {}).get(substrate_name, {})
        if isinstance(substrates_config, Mapping)
        else {}
    )
    substrate_smiles = str(substrate_config.get("smiles", "") or "").strip()
    reaction_note = str(
        substrate_config.get(
            "reaction_note",
            "CatPred preparation row only; no prediction is implied until CatPred output is ingested.",
        )
    )
    validation_by_id = {
        str(row.get("enzyme_validation_id", "")): row.to_dict()
        for _, row in validation.fillna("").iterrows()
    }

    rows: list[dict[str, object]] = []
    for _, manifest_row in manifest.fillna("").iterrows():
        manifest_record = manifest_row.to_dict()
        validation_record = validation_by_id.get(str(manifest_record.get("enzyme_validation_id", "")), {})
        status = "ready"
        if not substrate_smiles:
            status = "needs_substrate_definition"
        if not as_bool(manifest_record.get("catpred_ready", "")):
            status = "not_catpred_ready"
        rows.append(
            {
                "protein_accession": manifest_record.get("protein_accession", ""),
                "enzyme_validation_id": manifest_record.get("enzyme_validation_id", ""),
                "sequence_id": sequence_id_for(manifest_record),
                "substrate_name": substrate_name,
                "substrate_smiles": substrate_smiles,
                "reaction_note": reaction_note,
                "enzyme_class": manifest_record.get("enzyme_class", validation_record.get("enzyme_class", "")),
                "candidate_id": manifest_record.get("candidate_id", ""),
                "organism_label": manifest_record.get("organism_label", ""),
                "catpred_input_status": status,
            }
        )
    return pd.DataFrame(rows, columns=CATPRED_INPUT_COLUMNS)


def _scale_high(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce")
    max_value = values.max(skipna=True)
    if pd.isna(max_value) or max_value <= 0:
        return pd.Series([0.0] * len(values), index=series.index)
    return (values / max_value).clip(lower=0.0, upper=1.0).fillna(0.0)


def _scale_low(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce")
    positives = values[values > 0]
    if positives.empty:
        return pd.Series([0.0] * len(values), index=series.index)
    min_value = positives.min(skipna=True)
    return (min_value / values).clip(lower=0.0, upper=1.0).fillna(0.0)


def _validation_join_columns(validation: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "enzyme_validation_id",
        "candidate_id",
        "protein_accession",
        "organism_label",
        "enzyme_class",
        "accession_link_confidence",
        "catpred_ready",
        "manual_verified",
    ]
    available = [column for column in columns if column in validation.columns]
    return validation.loc[:, available].copy()


def ingest_catpred_results(
    catpred_output: pd.DataFrame,
    validation: pd.DataFrame,
    column_map: Mapping[str, str] | None = None,
) -> pd.DataFrame:
    column_map = {
        "protein_accession": "protein_accession",
        "predicted_kcat": "predicted_kcat",
        "predicted_Km": "predicted_Km",
        "predicted_kcat_Km": "predicted_kcat_Km",
        "prediction_confidence": "prediction_confidence",
        **(dict(column_map or {})),
    }
    required_keys = ["protein_accession", "predicted_kcat", "predicted_Km", "predicted_kcat_Km"]
    missing = [column_map[key] for key in required_keys if column_map[key] not in catpred_output.columns]
    if missing:
        raise ValueError(f"CatPred output missing required columns: {', '.join(missing)}")

    output = catpred_output.copy()
    rename = {column_map[key]: key for key in required_keys if column_map[key] in output.columns}
    if column_map.get("prediction_confidence") in output.columns:
        rename[column_map["prediction_confidence"]] = "prediction_confidence"
    output = output.rename(columns=rename)
    output["catpred_kcat_score"] = _scale_high(output["predicted_kcat"])
    output["catpred_Km_score"] = _scale_low(output["predicted_Km"])
    output["catpred_efficiency_score"] = _scale_high(output["predicted_kcat_Km"])
    score_columns = ["catpred_kcat_score", "catpred_Km_score", "catpred_efficiency_score"]
    output["catpred_overall_score"] = output[score_columns].mean(axis=1).round(4)
    output["catpred_result_available"] = True

    validation_join = _validation_join_columns(validation.fillna(""))
    merged = validation_join.merge(output, on="protein_accession", how="inner")
    if merged.empty:
        output["enzyme_validation_id"] = ""
        output["candidate_id"] = ""
        output["organism_label"] = ""
        output["enzyme_class"] = ""
        output["accession_link_confidence"] = "unlinked"
        output["catpred_ready"] = False
        output["manual_verified"] = False
        merged = output
    merged["catpred_validation_level"] = merged["accession_link_confidence"].apply(
        lambda value: "catpred_validated" if value in READY_LINK_CONFIDENCE else "blocked_link_confidence"
    )
    return merged


def write_private_csv(dataframe: pd.DataFrame, path: str | Path) -> None:
    resolved = Path(path)
    resolved.parent.mkdir(parents=True, exist_ok=True)
    dataframe.to_csv(resolved, index=False)
