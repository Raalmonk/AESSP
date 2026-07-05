#!/usr/bin/env python
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from aessp.io import ensure_parent, read_csv_or_empty


def _as_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "y", "verified", "ready"}


def _has_value(value: object) -> bool:
    text = str(value or "").strip().lower()
    return text not in {"", "nan", "none", "null", "unknown"}


def _only_mention_level(row: pd.Series) -> bool:
    levels = [
        row.get("evidence_level", ""),
        row.get("structure_evidence_level", ""),
        row.get("enzyme_evidence_level", ""),
        row.get("availability_evidence_level", ""),
    ]
    return all(str(level or "mention_only") == "mention_only" for level in levels)


def build_audit_report(scores: pd.DataFrame) -> str:
    scores = scores.fillna("")
    total = len(scores)
    empty_genus_species = int(
        ((scores.get("genus", "") == "") | (scores.get("species", "") == "")).sum()
    ) if total else 0
    mention_only = int(scores.apply(_only_mention_level, axis=1).sum()) if total else 0
    accession_no_literature = 0
    if total:
        accession_no_literature = int(
            scores.apply(
                lambda row: _has_value(row.get("protein_accession", ""))
                and float(row.get("literature_evidence_count", 0) or 0) <= 0,
                axis=1,
            ).sum()
        )
    manual_verified = int(scores.get("manual_verified", pd.Series(dtype=object)).apply(_as_bool).sum())
    pilot_ready = int(
        scores.apply(
            lambda row: _as_bool(row.get("pilot_ready", "")) or _as_bool(row.get("manual_verified", "")),
            axis=1,
        ).sum()
    ) if total else 0

    warnings: list[str] = []
    if empty_genus_species:
        warnings.append(f"{empty_genus_species} candidates have empty genus or species.")
    if mention_only:
        warnings.append(f"{mention_only} candidates have mention-only evidence levels.")
    if accession_no_literature:
        warnings.append(f"{accession_no_literature} candidates have protein accession but no literature evidence.")
    if pilot_ready < 8:
        warnings.append("Fewer than 8 candidates are pilot-ready; do not use Top 8 as screening recommendations.")
    if not warnings:
        warnings.append("No conservative-pipeline warnings detected.")

    lines = [
        "# First Batch Audit",
        "",
        "## Counts",
        "",
        f"- candidates: {total}",
        f"- empty genus/species: {empty_genus_species}",
        f"- only mention-level evidence: {mention_only}",
        f"- protein accession without literature evidence: {accession_no_literature}",
        f"- manual_verified=True: {manual_verified}",
        f"- pilot_ready: {pilot_ready}",
        "",
        "## Top 20 By Score",
        "",
        "| rank | candidate_id | organism_label | total_score | evidence_level | structure | enzyme | availability |",
        "| --- | --- | --- | ---: | --- | --- | --- | --- |",
    ]
    top20 = scores.sort_values("total_score", ascending=False, kind="mergesort").head(20) if total else scores
    for rank, (_, row) in enumerate(top20.iterrows(), start=1):
        lines.append(
            "| {rank} | {candidate_id} | {organism_label} | {total_score} | {evidence_level} | "
            "{structure} | {enzyme} | {availability} |".format(
                rank=rank,
                candidate_id=row.get("candidate_id", ""),
                organism_label=row.get("organism_label", ""),
                total_score=row.get("total_score", ""),
                evidence_level=row.get("evidence_level", ""),
                structure=row.get("structure_evidence_level", ""),
                enzyme=row.get("enzyme_evidence_level", ""),
                availability=row.get("availability_evidence_level", ""),
            )
        )
    lines.extend(["", "## Warnings", ""])
    lines.extend(f"- {warning}" for warning in warnings)
    lines.append("")
    return "\n".join(lines)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Audit first-batch AESSP candidate scores.")
    parser.add_argument("--scores-csv", default="data/processed/dextran_candidate_scores.csv")
    parser.add_argument("--out-md", default="data/reports/first_batch_audit.md")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    scores = read_csv_or_empty(args.scores_csv)
    report = build_audit_report(scores)
    ensure_parent(args.out_md)
    Path(args.out_md).write_text(report, encoding="utf-8")
    print(report)


if __name__ == "__main__":
    main()
