#!/usr/bin/env python
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from tqdm import tqdm

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from aessp.api.cache import FileCache
from aessp.api.ncbi import (
    NCBIClient,
    extract_esearch_ids,
    parse_fasta_hashes,
    parse_protein_summary_records,
)
from aessp.io import load_queries, utc_now_iso, write_csv_records, write_jsonl


PROTEIN_COLUMNS = [
    "query",
    "source_api",
    "accession",
    "protein_name",
    "gene_name",
    "organism",
    "taxonomy_id",
    "sequence_length",
    "sequence_hash",
    "ec_number",
    "reviewed_status",
    "source_url",
    "fetched_at",
    "notes",
]


def collect_protein_records(
    client: NCBIClient,
    query: str,
    max_results: int,
    fetched_at: str,
    fetch_fasta_hash: bool = False,
) -> list[dict[str, object]]:
    search_payload = client.esearch("protein", query, retmax=max_results)
    ids = extract_esearch_ids(search_payload)
    if not ids:
        return []
    summary_payload = client.esummary("protein", ids)
    records = parse_protein_summary_records(summary_payload, query=query, fetched_at=fetched_at)
    if fetch_fasta_hash:
        fasta_text = client.efetch_fasta("protein", ids)
        fasta_hashes = parse_fasta_hashes(fasta_text)
        for record in records:
            accession = str(record.get("accession", ""))
            hash_record = fasta_hashes.get(accession)
            if hash_record:
                record.update(hash_record)
                record["notes"] = "NCBI Protein metadata; sequence hashed but not written to reports"
    return records


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Collect dextran-related protein records from NCBI.")
    parser.add_argument("--queries", required=True, help="YAML query file.")
    parser.add_argument("--out-jsonl", required=True, help="Raw JSONL output path.")
    parser.add_argument("--out-csv", required=True, help="Processed CSV output path.")
    parser.add_argument("--email", default="", help="Contact email for NCBI policies.")
    parser.add_argument("--max-results-per-query", type=int, default=50)
    parser.add_argument("--cache-dir", default="data/cache")
    parser.add_argument(
        "--fetch-fasta-hash",
        action="store_true",
        help="Fetch FASTA only to compute sequence_hash; full sequences are not written.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    queries = load_queries(args.queries)
    cache = FileCache(args.cache_dir)
    client = NCBIClient(email=args.email, cache=cache)
    fetched_at = utc_now_iso()
    records: list[dict[str, object]] = []

    for entry in tqdm(queries, desc="protein queries"):
        query = entry["query"]
        try:
            records.extend(
                collect_protein_records(
                    client,
                    query,
                    args.max_results_per_query,
                    fetched_at,
                    fetch_fasta_hash=args.fetch_fasta_hash,
                )
            )
        except Exception as exc:
            print(f"NCBI Protein collection failed for query '{query}': {exc}", file=sys.stderr)

    normalized = [{column: record.get(column, "") for column in PROTEIN_COLUMNS} for record in records]
    write_jsonl(normalized, args.out_jsonl)
    write_csv_records(normalized, args.out_csv, columns=PROTEIN_COLUMNS)


if __name__ == "__main__":
    main()
