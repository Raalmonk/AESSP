#!/usr/bin/env python
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from aessp.io import read_csv_or_empty, read_yaml
from aessp.validation.catpred import build_catpred_input_table, write_private_csv


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare CatPred input table from validated enzyme manifest.")
    parser.add_argument("--validation-csv", default="data/processed/enzyme_validation_table.csv")
    parser.add_argument("--manifest-csv", default="data/private/catpred/catpred_input_manifest.csv")
    parser.add_argument("--substrates", default="configs/substrates.yaml")
    parser.add_argument("--substrate-name", default="sucrose")
    parser.add_argument("--out-csv", default="data/private/catpred/catpred_input_table.csv")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    validation = read_csv_or_empty(args.validation_csv)
    manifest = read_csv_or_empty(args.manifest_csv)
    substrates_config = read_yaml(args.substrates) if Path(args.substrates).exists() else {}
    table = build_catpred_input_table(validation, manifest, substrates_config, substrate_name=args.substrate_name)
    write_private_csv(table, args.out_csv)


if __name__ == "__main__":
    main()
