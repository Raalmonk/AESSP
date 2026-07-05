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
from aessp.validation.catpred import ingest_catpred_results


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Ingest real CatPred output and derive conservative scores.")
    parser.add_argument("--catpred-output", default="data/private/catpred/catpred_output.csv")
    parser.add_argument("--validation-csv", default="data/processed/enzyme_validation_table.csv")
    parser.add_argument("--out-csv", default="data/processed/enzyme_catpred_scores.csv")
    parser.add_argument("--accession-column", default="protein_accession")
    parser.add_argument("--kcat-column", default="predicted_kcat")
    parser.add_argument("--km-column", default="predicted_Km")
    parser.add_argument("--efficiency-column", default="predicted_kcat_Km")
    parser.add_argument("--confidence-column", default="prediction_confidence")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    catpred_output = pd.read_csv(args.catpred_output, keep_default_na=False)
    validation = read_csv_or_empty(args.validation_csv)
    column_map = {
        "protein_accession": args.accession_column,
        "predicted_kcat": args.kcat_column,
        "predicted_Km": args.km_column,
        "predicted_kcat_Km": args.efficiency_column,
        "prediction_confidence": args.confidence_column,
    }
    try:
        scores = ingest_catpred_results(catpred_output, validation, column_map=column_map)
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc
    ensure_parent(args.out_csv)
    scores.to_csv(args.out_csv, index=False)


if __name__ == "__main__":
    main()
