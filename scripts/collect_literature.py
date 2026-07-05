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
from aessp.api.crossref import CrossrefClient, parse_crossref_works
from aessp.api.ncbi import (
    NCBIClient,
    extract_esearch_ids,
    parse_pubmed_summary_records,
    parse_pubmed_xml_records,
)
from aessp.io import load_queries, utc_now_iso, write_csv_records, write_jsonl


LITERATURE_COLUMNS = [
    "query",
    "source_api",
    "title",
    "abstract",
    "authors",
    "year",
    "journal",
    "doi",
    "pmid",
    "url",
    "fetched_at",
    "source_record_id",
]


def collect_pubmed(
    client: NCBIClient,
    query: str,
    max_results: int,
    fetched_at: str,
) -> list[dict[str, object]]:
    search_payload = client.esearch("pubmed", query, retmax=max_results)
    ids = extract_esearch_ids(search_payload)
    if not ids:
        return []

    try:
        xml_text = client.efetch_pubmed_xml(ids)
        records = parse_pubmed_xml_records(xml_text, query=query, fetched_at=fetched_at)
        if records:
            return records
    except Exception as exc:  # pragma: no cover - fallback path depends on remote payload shape
        print(f"PubMed XML fetch failed for query '{query}': {exc}", file=sys.stderr)

    summary_payload = client.esummary("pubmed", ids)
    return parse_pubmed_summary_records(summary_payload, query=query, fetched_at=fetched_at)


def collect_crossref(
    client: CrossrefClient,
    query: str,
    max_results: int,
    fetched_at: str,
) -> list[dict[str, object]]:
    payload = client.query_works(query, rows=max_results)
    return parse_crossref_works(payload, query=query, fetched_at=fetched_at)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Collect dextran literature records from public APIs.")
    parser.add_argument("--queries", required=True, help="YAML query file.")
    parser.add_argument("--out-jsonl", required=True, help="Raw JSONL output path.")
    parser.add_argument("--out-csv", required=True, help="Processed CSV output path.")
    parser.add_argument("--email", default="", help="Contact email for API policies.")
    parser.add_argument("--max-results-per-query", type=int, default=20)
    parser.add_argument("--cache-dir", default="data/cache")
    parser.add_argument(
        "--sources",
        nargs="+",
        choices=["ncbi", "crossref"],
        default=["ncbi", "crossref"],
        help="APIs to query.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    queries = load_queries(args.queries)
    cache = FileCache(args.cache_dir)
    fetched_at = utc_now_iso()
    records: list[dict[str, object]] = []

    ncbi_client = NCBIClient(email=args.email, cache=cache) if "ncbi" in args.sources else None
    crossref_client = CrossrefClient(email=args.email, cache=cache) if "crossref" in args.sources else None

    for entry in tqdm(queries, desc="literature queries"):
        query = entry["query"]
        if ncbi_client is not None:
            try:
                records.extend(collect_pubmed(ncbi_client, query, args.max_results_per_query, fetched_at))
            except Exception as exc:
                print(f"NCBI PubMed collection failed for query '{query}': {exc}", file=sys.stderr)
        if crossref_client is not None:
            try:
                records.extend(collect_crossref(crossref_client, query, args.max_results_per_query, fetched_at))
            except Exception as exc:
                print(f"Crossref collection failed for query '{query}': {exc}", file=sys.stderr)

    normalized = [{column: record.get(column, "") for column in LITERATURE_COLUMNS} for record in records]
    write_jsonl(normalized, args.out_jsonl)
    write_csv_records(normalized, args.out_csv, columns=LITERATURE_COLUMNS)


if __name__ == "__main__":
    main()
