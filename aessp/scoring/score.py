from __future__ import annotations

from pathlib import Path
from typing import Mapping
import re

import pandas as pd

from aessp.io import ensure_parent
from aessp.validation.catpred import READY_LINK_CONFIDENCE


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


def _has_numeric(value: object) -> bool:
    return bool(re.search(r"\d", str(value or "")))


def _manual_verified(row: Mapping[str, object]) -> bool:
    return _as_bool(row.get("manual_verified", ""))


def _pilot_ready(row: Mapping[str, object]) -> bool:
    return _as_bool(row.get("pilot_ready", "")) or _manual_verified(row)


def _availability_score(row: Mapping[str, object]) -> float:
    if _as_bool(row.get("available_in_china", "")):
        return 1.0
    status = str(row.get("availability_status", "")).lower()
    if "confirmed" in status and "china" in status and "unavailable" not in status:
        return 1.0
    if "outside china" in status and "available" in status and "unavailable" not in status:
        return 0.7
    if status in {"manual_review_needed", "manual review needed"}:
        return 0.1
    if "pending" in status or "manual" in status:
        return 0.1
    return 0.0


def _safety_score(row: Mapping[str, object]) -> float:
    notes = str(row.get("safety_notes", "")).lower()
    if not notes.strip():
        return 0.0
    if any(term in notes for term in ("pathogen", "toxin", "hazard", "risk")):
        return 0.0
    low_concern = any(term in notes for term in ("gras", "food", "starter", "industrial", "low concern"))
    if low_concern and _manual_verified(row):
        return 1.0
    if low_concern:
        return 0.5
    return 0.0


def _product_structure_score(row: Mapping[str, object]) -> float:
    if _manual_verified(row) and any(
        _has_value(row.get(column, ""))
        for column in ("reported_branching", "reported_alpha_1_6", "reported_alpha_1_3")
    ):
        return 1.0
    if any(_has_numeric(row.get(column, "")) for column in ("reported_branching", "reported_alpha_1_6", "reported_alpha_1_3")):
        return 0.7

    evidence = " ".join(
        str(row.get(column, "") or "")
        for column in ("branching_evidence_text", "nmr_evidence_text")
    ).lower()
    if not evidence.strip():
        return 0.0
    if "nmr" in evidence or "linkage" in evidence or "1,6" in evidence or "1,3" in evidence:
        return 0.4
    return 0.2


def _target_mw_score(row: Mapping[str, object]) -> float:
    if _manual_verified(row) and _has_value(row.get("reported_Mw", "")):
        return 1.0
    if _has_numeric(row.get("reported_Mw", "")):
        return 0.7
    if _has_value(row.get("mw_evidence_text", "")):
        return 0.2
    return 0.0


def _branching_evidence_score(row: Mapping[str, object]) -> float:
    if _manual_verified(row) and any(
        _has_value(row.get(column, ""))
        for column in ("reported_branching", "reported_alpha_1_6", "reported_alpha_1_3")
    ):
        return 1.0
    if any(_has_numeric(row.get(column, "")) for column in ("reported_branching", "reported_alpha_1_6", "reported_alpha_1_3")):
        return 0.7

    evidence = " ".join(
        str(row.get(column, "") or "")
        for column in ("branching_evidence_text", "nmr_evidence_text")
    ).lower()
    if not evidence.strip():
        return 0.0
    if "nmr" in evidence or "linkage" in evidence or "1,6" in evidence or "1,3" in evidence:
        return 0.3
    return 0.2


def _enzyme_sequence_score(row: Mapping[str, object]) -> float:
    if _manual_verified(row) and _has_value(row.get("protein_accession", "")):
        return 1.0
    if not (
        _as_bool(row.get("protein_sequence_available", ""))
        or _has_value(row.get("protein_accession", ""))
    ):
        return 0.0
    if _as_float(row.get("literature_evidence_count", 0.0)) > 0:
        return 0.8
    return 0.5


def _component_scores(row: Mapping[str, object]) -> dict[str, float]:
    evidence_confidence = _clamp(_as_float(row.get("evidence_confidence", 0.0)))
    if evidence_confidence == 0.0 and _as_float(row.get("literature_evidence_count", 0.0)) > 0:
        evidence_confidence = _clamp(_as_float(row.get("literature_evidence_count", 0.0)) / 3.0)

    has_process = any(
        _has_value(row.get(column, ""))
        for column in ("reported_viscosity", "reported_yield", "viscosity_evidence_text", "yield_evidence_text")
    )
    return {
        "literature_confidence_score": evidence_confidence,
        "product_structure_score": _product_structure_score(row),
        "target_mw_score": _target_mw_score(row),
        "branching_evidence_score": _branching_evidence_score(row),
        "enzyme_sequence_score": _enzyme_sequence_score(row),
        "availability_score": _availability_score(row),
        "safety_score": _safety_score(row),
        "processability_score": 0.5 if has_process and not _manual_verified(row) else (1.0 if has_process else 0.0),
    }


def _catpred_index(catpred_scores: pd.DataFrame | None) -> tuple[dict[str, dict[str, object]], set[str]]:
    if catpred_scores is None or catpred_scores.empty:
        return {}, set()

    valid: dict[str, dict[str, object]] = {}
    blocked_candidates: set[str] = set()
    for _, row in catpred_scores.fillna("").iterrows():
        record = row.to_dict()
        candidate_id = str(record.get("candidate_id", "")).strip()
        if not candidate_id:
            continue
        link_confidence = str(record.get("accession_link_confidence", "")).strip()
        if link_confidence not in READY_LINK_CONFIDENCE:
            blocked_candidates.add(candidate_id)
            continue
        if not _as_bool(record.get("catpred_result_available", "")):
            continue
        score = _as_float(record.get("catpred_overall_score", ""), default=-1.0)
        if score < 0:
            continue
        current = valid.get(candidate_id)
        if current is None or score > float(current.get("catpred_overall_score", -1.0)):
            valid[candidate_id] = record
    return valid, blocked_candidates


def score_candidates(
    dataframe: pd.DataFrame,
    weights_config: Mapping[str, object],
    catpred_scores: pd.DataFrame | None = None,
) -> pd.DataFrame:
    weights = dict(weights_config.get("weights", {}) if isinstance(weights_config, Mapping) else {})
    uncertainty = dict(weights_config.get("uncertainty", {}) if isinstance(weights_config, Mapping) else {})
    missing_value_penalty = float(uncertainty.get("missing_value_penalty", 0.10))
    review_penalty_enabled = bool(uncertainty.get("auto_extracted_numeric_requires_review", True))
    catpred_by_candidate, catpred_blocked_candidates = _catpred_index(catpred_scores)

    rows: list[dict[str, object]] = []
    for _, row in dataframe.fillna("").iterrows():
        output = row.to_dict()
        candidate_id = str(output.get("candidate_id", "")).strip()
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
        evidence_priority_score = round(_clamp(weighted_total - uncertainty_penalty), 4)
        catpred_record = catpred_by_candidate.get(candidate_id)
        output["catpred_kinetic_score"] = ""
        output["catpred_result_available"] = False
        output["catpred_validation_level"] = "not_yet_validated"
        if catpred_record is not None:
            catpred_score = _clamp(_as_float(catpred_record.get("catpred_overall_score", 0.0)))
            output["catpred_kinetic_score"] = round(catpred_score, 4)
            output["catpred_result_available"] = True
            output["catpred_validation_level"] = "catpred_validated"
            production_priority_score = round(_clamp((0.75 * evidence_priority_score) + (0.25 * catpred_score)), 4)
        else:
            if candidate_id in catpred_blocked_candidates:
                output["catpred_validation_level"] = "blocked_link_confidence"
            production_priority_score = evidence_priority_score
        output["evidence_priority_score"] = evidence_priority_score
        output["catpred_validated_score"] = output["catpred_kinetic_score"]
        output["production_priority_score"] = production_priority_score
        output["total_score"] = production_priority_score
        notes: list[str] = []
        if _as_bool(output.get("needs_manual_review", "")):
            notes.append("manual review required for automatically extracted evidence")
        if not _pilot_ready(output):
            notes.append("not pilot-ready without manual verification")
        if not output["catpred_result_available"]:
            notes.append("CatPred result not yet validated")
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
    if "pilot_ready" in safe_scored.columns or "manual_verified" in safe_scored.columns:
        ready_mask = safe_scored.apply(lambda row: _pilot_ready(row.to_dict()), axis=1)
    else:
        ready_mask = pd.Series(False, index=safe_scored.index)

    ensure_parent(out_csv)
    safe_scored.to_csv(out_csv, index=False)

    ensure_parent(top20_path)
    safe_scored.head(20).to_csv(top20_path, index=False)

    top8_path = Path(top8_path)
    not_ready_path = top8_path.with_name("top8_pilot_screening_NOT_READY.md")
    ready_candidates = safe_scored.loc[ready_mask].head(8)
    if len(ready_candidates) >= 8:
        ensure_parent(top8_path)
        ready_candidates.to_csv(top8_path, index=False)
        if not_ready_path.exists():
            not_ready_path.unlink()
        return

    if top8_path.exists():
        top8_path.unlink()
    ensure_parent(not_ready_path)
    not_ready_path.write_text(
        "\n".join(
            [
                "# Top 8 Pilot Screening Not Ready",
                "",
                f"Only {len(ready_candidates)} candidates are manual_verified=True or pilot_ready=True.",
                "A pilot-screening recommendation file is not generated from auto-extracted evidence alone.",
                "Use top20_manual_review.csv for manual curation first.",
                "",
            ]
        ),
        encoding="utf-8",
    )
