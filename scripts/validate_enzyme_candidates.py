#!/usr/bin/env python
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from aessp.io import ensure_parent, read_csv_or_empty
from aessp.validation.catpred import build_enzyme_validation_table


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate linked enzyme accessions before CatPred prep.")
    parser.add_argument("--candidates-csv", default="data/processed/dextran_candidate_table.csv")
    parser.add_argument("--protein-csv", default="data/processed/protein_records.csv")
    parser.add_argument("--manual-validation-csv", default="data/manual/manual_enzyme_validation.csv")
    parser.add_argument("--out-csv", default="data/processed/enzyme_validation_table.csv")
    parser.add_argument(
        "--allow-fasta-export",
        action="store_true",
        help="Allow CatPred readiness when sequence_hash is missing but explicit FASTA export will be run.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    candidates = read_csv_or_empty(args.candidates_csv)
    proteins = read_csv_or_empty(args.protein_csv)
    manual_path = Path(args.manual_validation_csv)
    manual = read_csv_or_empty(manual_path) if manual_path.exists() else None
    validation = build_enzyme_validation_table(
        candidates,
        proteins,
        manual_validation=manual,
        allow_fasta_export=args.allow_fasta_export,
    )
    ensure_parent(args.out_csv)
    validation.to_csv(args.out_csv, index=False)


if __name__ == "__main__":
    main()
