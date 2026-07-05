from __future__ import annotations

import hashlib
import html
import os
import re
import xml.etree.ElementTree as ET
from io import StringIO
from typing import Iterable

import requests
from Bio import SeqIO

from aessp.api.cache import FileCache
from aessp.api.rate_limit import RateLimiter
from aessp.io import utc_now_iso


SUPPORTED_DBS = {"pubmed", "protein", "nuccore"}
BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def _as_id_list(ids: str | Iterable[str]) -> list[str]:
    if isinstance(ids, str):
        return [part.strip() for part in ids.split(",") if part.strip()]
    return [str(part).strip() for part in ids if str(part).strip()]


def _validate_db(db: str) -> None:
    if db not in SUPPORTED_DBS:
        supported = ", ".join(sorted(SUPPORTED_DBS))
        raise ValueError(f"Unsupported NCBI db '{db}'. Supported values: {supported}")


def _clean_text(value: object) -> str:
    if value is None:
        return ""
    text = html.unescape(str(value))
    text = re.sub(r"\s+", " ", text).strip()
    return text


def _element_text(element: ET.Element | None) -> str:
    if element is None:
        return ""
    return _clean_text("".join(element.itertext()))


def _first_year(*values: object) -> str:
    for value in values:
        match = re.search(r"\b(18|19|20)\d{2}\b", str(value or ""))
        if match:
            return match.group(0)
    return ""


def _article_id(item: dict, id_type: str) -> str:
    for article_id in item.get("articleids", []) or []:
        if str(article_id.get("idtype", "")).lower() == id_type.lower():
            return str(article_id.get("value", "")).strip()
    return ""


class NCBIClient:
    def __init__(
        self,
        email: str = "",
        tool: str = "aessp",
        api_key: str | None = None,
        cache: FileCache | None = None,
        rate_limiter: RateLimiter | None = None,
        session: requests.Session | None = None,
        base_url: str = BASE_URL,
        timeout: float = 30.0,
    ) -> None:
        self.email = email
        self.tool = tool
        self.api_key = api_key if api_key is not None else os.getenv("NCBI_API_KEY", "")
        self.cache = cache or FileCache()
        self.rate_limiter = rate_limiter or RateLimiter(requests_per_second=3.0)
        self.session = session or requests.Session()
        self.base_url = base_url.rstrip("/")
        self.timeout = timeout

    def _with_common_params(self, params: dict[str, object]) -> dict[str, object]:
        request_params = dict(params)
        request_params["tool"] = self.tool
        if self.email:
            request_params["email"] = self.email
        if self.api_key:
            request_params["api_key"] = self.api_key
        return request_params

    def _get_json(self, endpoint: str, params: dict[str, object]) -> dict:
        request_params = self._with_common_params(params)
        cached = self.cache.get_json("ncbi", endpoint, request_params)
        if cached is not None:
            return cached  # type: ignore[return-value]

        self.rate_limiter.wait()
        response = self.session.get(
            f"{self.base_url}/{endpoint}",
            params=request_params,
            timeout=self.timeout,
        )
        response.raise_for_status()
        payload = response.json()
        self.cache.set_json("ncbi", endpoint, request_params, payload)
        return payload

    def _get_text(self, endpoint: str, params: dict[str, object]) -> str:
        request_params = self._with_common_params(params)
        cached = self.cache.get_text("ncbi", endpoint, request_params)
        if cached is not None:
            return cached

        self.rate_limiter.wait()
        response = self.session.get(
            f"{self.base_url}/{endpoint}",
            params=request_params,
            timeout=self.timeout,
        )
        response.raise_for_status()
        text = response.text
        self.cache.set_text("ncbi", endpoint, request_params, text)
        return text

    def esearch(
        self,
        db: str,
        term: str,
        retmax: int = 20,
        retstart: int = 0,
        sort: str | None = None,
    ) -> dict:
        _validate_db(db)
        params: dict[str, object] = {
            "db": db,
            "retmax": int(retmax),
            "retmode": "json",
            "retstart": int(retstart),
            "term": term,
        }
        if sort:
            params["sort"] = sort
        return self._get_json("esearch.fcgi", params)

    def esummary(self, db: str, ids: str | Iterable[str]) -> dict:
        _validate_db(db)
        id_list = _as_id_list(ids)
        if not id_list:
            return {"result": {"uids": []}}
        params = {
            "db": db,
            "id": ",".join(id_list),
            "retmode": "json",
        }
        return self._get_json("esummary.fcgi", params)

    def efetch(
        self,
        db: str,
        ids: str | Iterable[str],
        rettype: str | None = None,
        retmode: str = "text",
    ) -> str:
        _validate_db(db)
        id_list = _as_id_list(ids)
        if not id_list:
            return ""
        params: dict[str, object] = {
            "db": db,
            "id": ",".join(id_list),
            "retmode": retmode,
        }
        if rettype:
            params["rettype"] = rettype
        return self._get_text("efetch.fcgi", params)

    def efetch_pubmed_xml(self, ids: str | Iterable[str]) -> str:
        return self.efetch("pubmed", ids, retmode="xml")

    def efetch_fasta(self, db: str, ids: str | Iterable[str]) -> str:
        return self.efetch(db, ids, rettype="fasta", retmode="text")


def extract_esearch_ids(payload: dict) -> list[str]:
    return [str(item) for item in payload.get("esearchresult", {}).get("idlist", [])]


def parse_pubmed_summary_records(
    payload: dict,
    query: str = "",
    fetched_at: str | None = None,
) -> list[dict[str, object]]:
    fetched = fetched_at or utc_now_iso()
    result = payload.get("result", {})
    records: list[dict[str, object]] = []
    for uid in result.get("uids", []) or []:
        item = result.get(str(uid), {})
        authors = "; ".join(
            _clean_text(author.get("name", ""))
            for author in item.get("authors", []) or []
            if author.get("name")
        )
        doi = _article_id(item, "doi")
        pmid = _article_id(item, "pubmed") or str(uid)
        records.append(
            {
                "query": query,
                "source_api": "ncbi_pubmed",
                "title": _clean_text(item.get("title", "")),
                "abstract": "",
                "authors": authors,
                "year": _first_year(item.get("pubdate"), item.get("epubdate")),
                "journal": _clean_text(item.get("fulljournalname") or item.get("source", "")),
                "doi": doi,
                "pmid": pmid,
                "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else "",
                "fetched_at": fetched,
                "source_record_id": pmid,
            }
        )
    return records


def parse_pubmed_xml_records(
    xml_text: str,
    query: str = "",
    fetched_at: str | None = None,
) -> list[dict[str, object]]:
    if not xml_text.strip():
        return []
    fetched = fetched_at or utc_now_iso()
    root = ET.fromstring(xml_text)
    records: list[dict[str, object]] = []
    for article in root.findall(".//PubmedArticle"):
        medline = article.find("./MedlineCitation")
        article_node = medline.find("./Article") if medline is not None else None
        pmid = _element_text(medline.find("./PMID") if medline is not None else None)
        title = _element_text(article_node.find("./ArticleTitle") if article_node is not None else None)
        abstract_parts: list[str] = []
        if article_node is not None:
            for abstract_node in article_node.findall("./Abstract/AbstractText"):
                label = abstract_node.attrib.get("Label", "").strip()
                text = _element_text(abstract_node)
                if text and label:
                    abstract_parts.append(f"{label}: {text}")
                elif text:
                    abstract_parts.append(text)
        journal = _element_text(article_node.find("./Journal/Title") if article_node is not None else None)
        pub_date = article_node.find("./Journal/JournalIssue/PubDate") if article_node is not None else None
        year = ""
        if pub_date is not None:
            year = _first_year(
                _element_text(pub_date.find("./Year")),
                _element_text(pub_date.find("./MedlineDate")),
            )
        if not year and article_node is not None:
            year = _first_year(_element_text(article_node.find("./ArticleDate/Year")))

        author_names: list[str] = []
        if article_node is not None:
            for author in article_node.findall("./AuthorList/Author"):
                collective = _element_text(author.find("./CollectiveName"))
                if collective:
                    author_names.append(collective)
                    continue
                last = _element_text(author.find("./LastName"))
                fore = _element_text(author.find("./ForeName"))
                name = " ".join(part for part in [fore, last] if part)
                if name:
                    author_names.append(name)

        doi = ""
        pubmed_data = article.find("./PubmedData/ArticleIdList")
        if pubmed_data is not None:
            for article_id in pubmed_data.findall("./ArticleId"):
                if article_id.attrib.get("IdType", "").lower() == "doi":
                    doi = _element_text(article_id)
                    break

        records.append(
            {
                "query": query,
                "source_api": "ncbi_pubmed",
                "title": title,
                "abstract": " ".join(abstract_parts),
                "authors": "; ".join(author_names),
                "year": year,
                "journal": journal,
                "doi": doi,
                "pmid": pmid,
                "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else "",
                "fetched_at": fetched,
                "source_record_id": pmid,
            }
        )
    return records


def parse_protein_summary_records(
    payload: dict,
    query: str = "",
    fetched_at: str | None = None,
) -> list[dict[str, object]]:
    fetched = fetched_at or utc_now_iso()
    result = payload.get("result", {})
    records: list[dict[str, object]] = []
    for uid in result.get("uids", []) or []:
        item = result.get(str(uid), {})
        title = _clean_text(item.get("title", ""))
        accession = _clean_text(item.get("accessionversion") or item.get("caption") or uid)
        organism = _clean_text(item.get("organism", ""))
        if not organism:
            bracketed = re.findall(r"\[([^\[\]]+)\]", title)
            organism = bracketed[-1].strip() if bracketed else ""
        protein_name = re.sub(r"\s*\[[^\[\]]+\]\s*$", "", title).strip()
        protein_name = re.sub(r"^MULTISPECIES:\s*", "", protein_name, flags=re.IGNORECASE)
        ec_match = re.search(r"\bEC[:\s]*(\d+|-)\.(\d+|-)\.(\d+|-)\.(\d+|-)\b", title)
        ec_number = ".".join(ec_match.groups()) if ec_match else ""
        sequence_length = item.get("slen", "")
        records.append(
            {
                "query": query,
                "source_api": "ncbi_protein",
                "accession": accession,
                "protein_name": protein_name,
                "gene_name": _clean_text(item.get("genename", "")),
                "organism": organism,
                "taxonomy_id": _clean_text(item.get("taxid", "")),
                "sequence_length": sequence_length,
                "sequence_hash": "",
                "ec_number": ec_number,
                "reviewed_status": "not_available_ncbi",
                "source_url": f"https://www.ncbi.nlm.nih.gov/protein/{accession}" if accession else "",
                "fetched_at": fetched,
                "notes": f"NCBI Protein UID {uid}; full sequence not written to reports",
            }
        )
    return records


def sequence_sha256(sequence: str) -> str:
    compact = re.sub(r"\s+", "", sequence).upper()
    return hashlib.sha256(compact.encode("ascii")).hexdigest() if compact else ""


def parse_fasta_hashes(fasta_text: str) -> dict[str, dict[str, object]]:
    hashes: dict[str, dict[str, object]] = {}
    for record in SeqIO.parse(StringIO(fasta_text), "fasta"):
        sequence = str(record.seq)
        accession = record.id.split("|")[-1]
        hashes[accession] = {
            "sequence_hash": sequence_sha256(sequence),
            "sequence_length": len(sequence),
        }
        hashes[record.id] = hashes[accession]
    return hashes
