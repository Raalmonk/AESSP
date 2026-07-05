#!/usr/bin/env python
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from aessp.extraction.mentions import extract_mentions_from_text
from aessp.io import read_csv_or_empty, write_csv_records


MENTION_COLUMNS = [
    "record_id",
    "title",
    "species_mentioned",
    "strain_mentioned",
    "enzyme_terms",
    "product_terms",
    "mw_mentions",
    "branching_mentions",
    "nmr_mentions",
    "viscosity_mentions",
    "yield_mentions",
    "confidence",
    "needs_manual_review",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Extract rule-based candidate mentions from literature CSV.")
    parser.add_argument(
        "--literature-csv",
        default="data/processed/literature_records.csv",
        help="Input literature_records.csv path.",
    )
    parser.add_argument(
        "--out-csv",
        default="data/processed/candidate_mentions.csv",
        help="Output candidate_mentions.csv path.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    literature = read_csv_or_empty(args.literature_csv)
    records: list[dict[str, object]] = []
    for index, row in literature.iterrows():
        record_id = row.get("source_record_id", "") or row.get("pmid", "") or row.get("doi", "") or index
        records.append(
            extract_mentions_from_text(
                title=str(row.get("title", "")),
                abstract=str(row.get("abstract", "")),
                record_id=record_id,
            )
        )
    write_csv_records(records, args.out_csv, columns=MENTION_COLUMNS)


if __name__ == "__main__":
    main()
