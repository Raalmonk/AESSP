from __future__ import annotations

import html
import re

import requests

from aessp.api.cache import FileCache
from aessp.io import utc_now_iso


BASE_URL = "https://api.crossref.org"


def _clean_text(value: object) -> str:
    if value is None:
        return ""
    text = html.unescape(str(value))
    text = re.sub(r"<[^>]+>", " ", text)
    text = re.sub(r"\s+", " ", text).strip()
    return text


def _first_list_value(value: object) -> str:
    if isinstance(value, list) and value:
        return _clean_text(value[0])
    return _clean_text(value)


def _crossref_year(item: dict) -> str:
    for key in ("published-print", "published-online", "published", "issued", "created"):
        date_parts = item.get(key, {}).get("date-parts", [])
        if date_parts and date_parts[0]:
            year = date_parts[0][0]
            if year:
                return str(year)
    return ""


class CrossrefClient:
    def __init__(
        self,
        email: str = "",
        cache: FileCache | None = None,
        session: requests.Session | None = None,
        base_url: str = BASE_URL,
        timeout: float = 30.0,
        user_agent: str | None = None,
    ) -> None:
        self.email = email
        self.cache = cache or FileCache()
        self.session = session or requests.Session()
        self.base_url = base_url.rstrip("/")
        self.timeout = timeout
        self.user_agent = user_agent

    def _headers(self) -> dict[str, str]:
        if self.user_agent:
            return {"User-Agent": self.user_agent}
        if self.email:
            return {"User-Agent": f"AESSP/0.1 (mailto:{self.email})"}
        return {}

    def query_works(self, query: str, rows: int = 20, offset: int = 0) -> dict:
        params: dict[str, object] = {
            "offset": int(offset),
            "query": query,
            "rows": int(rows),
        }
        if self.email:
            params["mailto"] = self.email
        cached = self.cache.get_json("crossref", "works", params)
        if cached is not None:
            return cached  # type: ignore[return-value]

        response = self.session.get(
            f"{self.base_url}/works",
            params=params,
            headers=self._headers(),
            timeout=self.timeout,
        )
        response.raise_for_status()
        payload = response.json()
        self.cache.set_json("crossref", "works", params, payload)
        return payload


def parse_crossref_works(
    payload: dict,
    query: str = "",
    fetched_at: str | None = None,
) -> list[dict[str, object]]:
    fetched = fetched_at or utc_now_iso()
    items = payload.get("message", {}).get("items", []) or []
    records: list[dict[str, object]] = []
    for item in items:
        doi = _clean_text(item.get("DOI", ""))
        authors: list[str] = []
        for author in item.get("author", []) or []:
            name = " ".join(
                part
                for part in [
                    _clean_text(author.get("given", "")),
                    _clean_text(author.get("family", "")),
                ]
                if part
            )
            if not name:
                name = _clean_text(author.get("name", ""))
            if name:
                authors.append(name)
        url = _clean_text(item.get("URL", ""))
        if not url and doi:
            url = f"https://doi.org/{doi}"
        records.append(
            {
                "query": query,
                "source_api": "crossref",
                "title": _first_list_value(item.get("title", "")),
                "abstract": _clean_text(item.get("abstract", "")),
                "authors": "; ".join(authors),
                "year": _crossref_year(item),
                "journal": _first_list_value(item.get("container-title", "")),
                "doi": doi,
                "pmid": "",
                "url": url,
                "fetched_at": fetched,
                "source_record_id": doi or url,
            }
        )
    return records
