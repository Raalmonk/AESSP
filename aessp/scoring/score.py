from __future__ import annotations

from pathlib import Path
from typing import Mapping

import pandas as pd

from aessp.io import ensure_parent


COMPONENT_COLUMNS = [
    "literature_confidence_score",
    "product_structure_score",
    "target_mw_score",
    "branching_evidence_score",
    "enzyme_sequence_score",
    "availability_score",
    "safety_score",
    "processability_score",
]


def _has_value(value: object) -> bool:
    if value is None:
        return False
    if pd.isna(value):
        return False
    text = str(value).strip()
    return text.lower() not in {"", "nan", "none", "null", "unknown", "not_available_ncbi"}


def _as_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    return text in {"1", "true", "yes", "y", "available", "present"}


def _as_float(value: object, default: float = 0.0) -> float:
    try:
        if not _has_value(value):
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


def _clamp(value: float) -> float:
    return max(0.0, min(1.0, value))


def _availability_score(row: Mapping[str, object]) -> float:
    if _as_bool(row.get("available_in_china", "")):
        return 1.0
    status = str(row.get("availability_status", "")).lower()
    if "available" in status and "unavailable" not in status:
        return 0.8
    if "pending" in status or "manual" in status:
        return 0.3
    return 0.0


def _safety_score(row: Mapping[str, object]) -> float:
    notes = str(row.get("safety_notes", "")).lower()
    if not notes.strip():
        return 0.0
    if any(term in notes for term in ("pathogen", "toxin", "hazard", "risk")):
        return 0.0
    if any(term in notes for term in ("gras", "food", "starter", "low concern")):
        return 1.0
    return 0.5


def _component_scores(row: Mapping[str, object]) -> dict[str, float]:
    evidence_confidence = _clamp(_as_float(row.get("evidence_confidence", 0.0)))
    if evidence_confidence == 0.0 and _as_float(row.get("literature_evidence_count", 0.0)) > 0:
        evidence_confidence = _clamp(_as_float(row.get("literature_evidence_count", 0.0)) / 3.0)

    has_branching = any(
        _has_value(row.get(column, ""))
        for column in ("reported_branching", "reported_alpha_1_6", "reported_alpha_1_3")
    )
    has_mw = _has_value(row.get("reported_Mw", ""))
    has_sequence = _as_bool(row.get("protein_sequence_available", "")) or _has_value(
        row.get("protein_accession", "")
    )
    has_process = _has_value(row.get("reported_viscosity", "")) or _has_value(
        row.get("reported_yield", "")
    )
    return {
        "literature_confidence_score": evidence_confidence,
        "product_structure_score": 1.0 if has_branching else (0.5 if has_mw else 0.0),
        "target_mw_score": 1.0 if has_mw else 0.0,
        "branching_evidence_score": 1.0 if has_branching else 0.0,
        "enzyme_sequence_score": 1.0 if has_sequence else 0.0,
        "availability_score": _availability_score(row),
        "safety_score": _safety_score(row),
        "processability_score": 1.0 if has_process else 0.0,
    }


def score_candidates(dataframe: pd.DataFrame, weights_config: Mapping[str, object]) -> pd.DataFrame:
    weights = dict(weights_config.get("weights", {}) if isinstance(weights_config, Mapping) else {})
    uncertainty = dict(weights_config.get("uncertainty", {}) if isinstance(weights_config, Mapping) else {})
    missing_value_penalty = float(uncertainty.get("missing_value_penalty", 0.10))
    review_penalty_enabled = bool(uncertainty.get("auto_extracted_numeric_requires_review", True))

    rows: list[dict[str, object]] = []
    for _, row in dataframe.fillna("").iterrows():
        output = row.to_dict()
        components = _component_scores(output)
        missing_columns = [
            "reported_Mw",
            "reported_branching",
            "protein_accession",
            "availability_status",
            "safety_notes",
        ]
        missing_fraction = sum(not _has_value(output.get(column, "")) for column in missing_columns) / len(
            missing_columns
        )
        uncertainty_penalty = missing_value_penalty * missing_fraction
        if review_penalty_enabled and _as_bool(output.get("needs_manual_review", "")):
            uncertainty_penalty += missing_value_penalty

        weighted_total = sum(
            components[column] * float(weights.get(column, 0.0)) for column in COMPONENT_COLUMNS
        )
        output.update(components)
        output["uncertainty_penalty"] = round(uncertainty_penalty, 4)
        output["total_score"] = round(_clamp(weighted_total - uncertainty_penalty), 4)
        notes: list[str] = []
        if _as_bool(output.get("needs_manual_review", "")):
            notes.append("manual review required for automatically extracted evidence")
        if missing_fraction:
            notes.append("missing values reduced score")
        output["score_notes"] = "; ".join(notes)
        rows.append(output)

    scored = pd.DataFrame(rows)
    if "total_score" in scored:
        scored = scored.sort_values("total_score", ascending=False, kind="mergesort")
    return scored


def report_safe_dataframe(dataframe: pd.DataFrame) -> pd.DataFrame:
    forbidden_exact = {"sequence", "protein_sequence", "full_sequence", "fasta", "raw_fasta"}
    keep_columns: list[str] = []
    for column in dataframe.columns:
        lowered = column.lower()
        if lowered in forbidden_exact or "fasta" in lowered or "full_sequence" in lowered:
            continue
        keep_columns.append(column)
    return dataframe.loc[:, keep_columns]


def write_score_outputs(
    scored: pd.DataFrame,
    out_csv: str | Path,
    top20_path: str | Path,
    top8_path: str | Path,
) -> None:
    safe_scored = report_safe_dataframe(scored)
    ensure_parent(out_csv)
    safe_scored.to_csv(out_csv, index=False)

    ensure_parent(top20_path)
    safe_scored.head(20).to_csv(top20_path, index=False)

    ensure_parent(top8_path)
    safe_scored.head(8).to_csv(top8_path, index=False)
