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
from aessp.validation.catpred import as_bool


def _write_csv(dataframe: pd.DataFrame, path: str | Path) -> None:
    ensure_parent(path)
    dataframe.to_csv(path, index=False)


def _summary(validation: pd.DataFrame, catpred_scores: pd.DataFrame) -> str:
    validation = validation.fillna("")
    with_accession = int((validation.get("protein_accession", "").astype(str).str.strip() != "").sum()) if len(validation) else 0
    link_counts = validation.get("accession_link_confidence", pd.Series(dtype=object)).value_counts().to_dict()
    ready_count = int(validation.get("catpred_ready", pd.Series(dtype=object)).apply(as_bool).sum())
    blocked = validation.loc[~validation.get("catpred_ready", pd.Series(dtype=object)).apply(as_bool)] if len(validation) else validation
    result_count = int(catpred_scores.get("catpred_result_available", pd.Series(dtype=object)).apply(as_bool).sum()) if len(catpred_scores) else 0

    reason_counts: dict[str, int] = {}
    for reason_text in blocked.get("catpred_block_reason", []):
        for reason in str(reason_text or "unclassified").split(";"):
            cleaned = reason.strip() or "unclassified"
            reason_counts[cleaned] = reason_counts.get(cleaned, 0) + 1

    lines = [
        "# CatPred Validation Summary",
        "",
        "CatPred results are reported only when actual CatPred output has been ingested.",
        "",
        "## Counts",
        "",
        f"- candidate rows: {len(validation)}",
        f"- rows with protein accessions: {with_accession}",
        f"- exact_strain: {int(link_counts.get('exact_strain', 0))}",
        f"- exact_species: {int(link_counts.get('exact_species', 0))}",
        f"- genus_only: {int(link_counts.get('genus_only', 0))}",
        f"- CatPred-ready: {ready_count}",
        f"- blocked: {len(blocked)}",
        f"- CatPred results ingested: {result_count}",
        "",
        "## Blocked Reasons",
        "",
    ]
    if reason_counts:
        lines.extend(f"- {reason}: {count}" for reason, count in sorted(reason_counts.items()))
    else:
        lines.append("- none")

    lines.extend(["", "## Top CatPred-Validated Enzymes", ""])
    if catpred_scores.empty:
        lines.append("No CatPred result rows have been ingested.")
    else:
        top = catpred_scores.sort_values("catpred_overall_score", ascending=False, kind="mergesort").head(10)
        lines.append("| protein_accession | candidate_id | enzyme_class | catpred_overall_score |")
        lines.append("| --- | --- | --- | ---: |")
        for _, row in top.iterrows():
            lines.append(
                f"| {row.get('protein_accession', '')} | {row.get('candidate_id', '')} | "
                f"{row.get('enzyme_class', '')} | {row.get('catpred_overall_score', '')} |"
            )
    lines.append("")
    return "\n".join(lines)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Report CatPred validation readiness and ingested results.")
    parser.add_argument("--validation-csv", default="data/processed/enzyme_validation_table.csv")
    parser.add_argument("--catpred-scores-csv", default="data/processed/enzyme_catpred_scores.csv")
    parser.add_argument("--ready-out", default="data/reports/catpred_ready_enzymes.csv")
    parser.add_argument("--blocked-out", default="data/reports/catpred_blocked_enzymes.csv")
    parser.add_argument("--validated-out", default="data/reports/catpred_validated_candidates.csv")
    parser.add_argument("--summary-out", default="data/reports/catpred_validation_summary.md")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    validation = read_csv_or_empty(args.validation_csv)
    catpred_scores = read_csv_or_empty(args.catpred_scores_csv)
    ready_mask = validation.get("catpred_ready", pd.Series(dtype=object)).apply(as_bool) if len(validation) else pd.Series(dtype=bool)
    ready = validation.loc[ready_mask].copy() if len(validation) else validation
    blocked = validation.loc[~ready_mask].copy() if len(validation) else validation
    if len(catpred_scores):
        exact_mask = catpred_scores.get("accession_link_confidence", pd.Series(dtype=object)).isin(["exact_strain", "exact_species"])
        validated = catpred_scores.loc[exact_mask].copy()
    else:
        validated = catpred_scores

    _write_csv(ready, args.ready_out)
    _write_csv(blocked, args.blocked_out)
    _write_csv(validated, args.validated_out)
    ensure_parent(args.summary_out)
    Path(args.summary_out).write_text(_summary(validation, validated), encoding="utf-8")


if __name__ == "__main__":
    main()
