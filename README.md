# AESSP

AESSP is an API-first evidence pipeline for scouting dextran-producing strains and dextran-related enzymes. It collects public metadata through official APIs, extracts conservative rule-based mentions, builds a candidate table, and scores candidates for manual review and pilot-screening prioritization.

The project is intentionally not a biological truth engine. Automatically extracted yield, molecular-weight, branching, viscosity, or linkage evidence is treated as a lead for expert review, not as a verified fact.

## Quickstart

```bash
python -m pip install -e .
python -m pytest
```

Collect literature metadata:

```bash
python scripts/collect_literature.py \
  --queries configs/literature_queries.yaml \
  --out-jsonl data/raw/literature_records.jsonl \
  --out-csv data/processed/literature_records.csv \
  --email USER_EMAIL \
  --max-results-per-query 20
```

Collect protein metadata:

```bash
python scripts/collect_proteins.py \
  --queries configs/protein_queries.yaml \
  --out-jsonl data/raw/protein_records.jsonl \
  --out-csv data/processed/protein_records.csv \
  --email USER_EMAIL \
  --max-results-per-query 50
```

Build extracted evidence, candidates, and scores:

```bash
python scripts/extract_candidate_mentions.py
python scripts/build_candidate_table.py
python scripts/score_candidates.py
python scripts/audit_first_batch.py
```

Primary outputs are written under `data/processed/` and `data/reports/`.

`build_candidate_table.py` now writes two layers:

- `data/processed/dextran_candidate_table.csv` for strain candidates with literature product evidence.
- `data/processed/enzyme_candidate_table.csv` for protein/accession candidates that still need literature or manual linkage.

Automatic title/abstract extraction fills `mw_evidence_text`, `branching_evidence_text`,
`viscosity_evidence_text`, `yield_evidence_text`, and `nmr_evidence_text`.
It does not fill `reported_*` fields. Those are reserved for manual or otherwise
structured verification.

`score_candidates.py` always writes `dextran_candidate_scores.csv` and
`top20_manual_review.csv`. It writes real Top 8 pilot-screening recommendations
only when at least 8 candidates have `manual_verified=True` or `pilot_ready=True`;
otherwise it writes `top8_pilot_screening_NOT_READY.md`.

## API Compliance

- Uses NCBI E-utilities and Crossref REST APIs rather than scraping pages.
- NCBI requests include `tool`, optional `email`, and optional `NCBI_API_KEY`.
- NCBI requests are rate-limited to 3 requests per second by default.
- Crossref requests use a polite `User-Agent` and `mailto` parameter when `--email` is provided.
- Responses are cached under `data/cache/`; API keys and other secret-like parameters are excluded from cache keys.
- Login-protected, paywalled, or supplier-private pages are out of scope.

## Limitations

- The extraction layer is rule-based and conservative.
- Numeric evidence from titles or abstracts is always marked `needs_manual_review=True`.
- Mention-only evidence is intentionally scored conservatively and cannot produce pilot-ready recommendations by itself.
- Full protein sequences are not written to reports.
- NCBI Protein does not provide UniProt-style reviewed status, so that field is marked unavailable for NCBI records.
- Scores are prioritization aids for manual review and pilot screening, not final strain or process recommendations.
