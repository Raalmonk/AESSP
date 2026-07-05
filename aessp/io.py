from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import pandas as pd
import yaml


def ensure_parent(path: str | Path) -> Path:
    resolved = Path(path)
    resolved.parent.mkdir(parents=True, exist_ok=True)
    return resolved


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def read_yaml(path: str | Path) -> dict:
    with Path(path).open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def load_queries(path: str | Path) -> list[dict[str, str]]:
    data = read_yaml(path)
    queries = data.get("queries", [])
    normalized: list[dict[str, str]] = []
    for entry in queries:
        if isinstance(entry, str):
            normalized.append({"name": entry, "query": entry})
        elif isinstance(entry, Mapping) and entry.get("query"):
            query = str(entry["query"])
            normalized.append({"name": str(entry.get("name", query)), "query": query})
    return normalized


def normalize_cell(value: object) -> object:
    if isinstance(value, (list, tuple, set)):
        return "; ".join(str(item) for item in value if item not in (None, ""))
    if value is None:
        return ""
    return value


def normalize_record(record: Mapping[str, object], columns: Sequence[str] | None = None) -> dict:
    keys = columns if columns is not None else record.keys()
    return {key: normalize_cell(record.get(key, "")) for key in keys}


def write_jsonl(records: Iterable[Mapping[str, object]], path: str | Path) -> None:
    resolved = ensure_parent(path)
    with resolved.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(json.dumps(record, ensure_ascii=False, sort_keys=True) + "\n")


def write_csv_records(
    records: Iterable[Mapping[str, object]],
    path: str | Path,
    columns: Sequence[str] | None = None,
) -> pd.DataFrame:
    rows = [normalize_record(record, columns=columns) for record in records]
    dataframe = pd.DataFrame(rows, columns=columns)
    ensure_parent(path)
    dataframe.to_csv(path, index=False)
    return dataframe


def read_csv_or_empty(path: str | Path, columns: Sequence[str] | None = None) -> pd.DataFrame:
    resolved = Path(path)
    if not resolved.exists():
        return pd.DataFrame(columns=columns)
    return pd.read_csv(resolved, keep_default_na=False)


def append_unique(values: Iterable[object]) -> str:
    seen: set[str] = set()
    output: list[str] = []
    for value in values:
        if value is None:
            continue
        for part in str(value).split(";"):
            cleaned = part.strip()
            if cleaned and cleaned not in seen:
                seen.add(cleaned)
                output.append(cleaned)
    return "; ".join(output)
