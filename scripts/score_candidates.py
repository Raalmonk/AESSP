#!/usr/bin/env python
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from aessp.io import read_csv_or_empty, read_yaml
from aessp.scoring.score import score_candidates, write_score_outputs


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Score dextran candidates and write review reports.")
    parser.add_argument("--candidates-csv", default="data/processed/dextran_candidate_table.csv")
    parser.add_argument("--weights", default="configs/scoring_weights.yaml")
    parser.add_argument("--out-csv", default="data/processed/dextran_candidate_scores.csv")
    parser.add_argument("--top20", default="data/reports/top20_manual_review.csv")
    parser.add_argument("--top8", default="data/reports/top8_pilot_screening_recommendations.csv")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    candidates = read_csv_or_empty(args.candidates_csv)
    weights_config = read_yaml(args.weights)
    scored = score_candidates(candidates, weights_config)
    write_score_outputs(scored, args.out_csv, args.top20, args.top8)


if __name__ == "__main__":
    main()
