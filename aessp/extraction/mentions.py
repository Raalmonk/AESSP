from __future__ import annotations

import re
from typing import Iterable, Mapping


SPECIES_RE = re.compile(
    r"\b(Leuconostoc|Weissella|Lactobacillus|Limosilactobacillus|Levilactobacillus|"
    r"Lentilactobacillus|Streptococcus|Pediococcus|Lactococcus|Oenococcus|"
    r"Fructobacillus|Liquorilactobacillus)\s+([a-z][a-z-]+)\b"
)
STRAIN_PATTERNS = [
    re.compile(r"\bNRRL\s+[A-Z]-?\d+[A-Za-z0-9-]*\b", re.IGNORECASE),
    re.compile(
        r"\b(?:ATCC|DSMZ?|JCM|KCTC|LMG|CGMCC|CICC|NBRC|NCIMB|NCTC|TMW|KACC|BCC|MTCC)"
        r"\s*[A-Z]?\s*\d+[A-Za-z0-9-]*\b",
        re.IGNORECASE,
    ),
    re.compile(r"\bstrain\s+[A-Za-z0-9][A-Za-z0-9._/-]{1,24}\b", re.IGNORECASE),
]
NUMBER_RE = re.compile(
    r"\b\d+(?:[.,]\d+)?(?:\s*(?:x|×)\s*10\^?-?\d+)?\b|\b10\^?-?\d+\b",
    re.IGNORECASE,
)

ENZYME_TERMS = [
    "dextransucrase",
    "glucansucrase",
    "alternansucrase",
    "mutansucrase",
    "reuteransucrase",
    "levansucrase",
    "GH70",
]
PRODUCT_TERMS = [
    "dextran",
    "glucan",
    "exopolysaccharide",
    "EPS",
    "oligosaccharide",
    "polysaccharide",
]
MW_KEYWORDS = [
    "molecular weight",
    "molar mass",
    "weight-average",
    "weight average",
    "mw",
    "m_w",
    "kda",
    "mda",
    "da",
    "g/mol",
]
BRANCHING_KEYWORDS = [
    "branching",
    "branched",
    "branch",
    "linkage",
    "alpha-1,6",
    "alpha-1,3",
    "alpha-(1,6)",
    "alpha-(1,3)",
    "alpha 1,6",
    "alpha 1,3",
    "degree of branching",
]
NMR_KEYWORDS = ["nmr", "1h nmr", "13c nmr", "nuclear magnetic resonance"]
VISCOSITY_KEYWORDS = ["viscosity", "viscous", "rheology", "rheological"]
YIELD_KEYWORDS = ["yield", "g/l", "g l-1", "g l^-1", "productivity"]


def _clean(value: object) -> str:
    text = re.sub(r"\s+", " ", str(value or "")).strip()
    return text


def _unique(values: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    output: list[str] = []
    for value in values:
        cleaned = _clean(value)
        key = cleaned.lower()
        if cleaned and key not in seen:
            seen.add(key)
            output.append(cleaned)
    return output


def _join(values: Iterable[str]) -> str:
    return "; ".join(_unique(values))


def _sentences(text: str) -> list[str]:
    parts = re.split(r"(?<=[.!?])\s+", text)
    return [_clean(part) for part in parts if _clean(part)]


def _term_mentions(text: str, terms: Iterable[str]) -> list[str]:
    mentions: list[str] = []
    for term in terms:
        if re.search(rf"(?<![A-Za-z0-9]){re.escape(term)}(?![A-Za-z0-9])", text, re.IGNORECASE):
            mentions.append(term)
    return _unique(mentions)


def _keyword_sentences(
    text: str,
    keywords: Iterable[str],
    require_number: bool = False,
    max_length: int = 240,
) -> list[str]:
    matches: list[str] = []
    lowered_keywords = [keyword.lower() for keyword in keywords]
    for sentence in _sentences(text):
        lowered = sentence.lower()
        if not any(keyword in lowered for keyword in lowered_keywords):
            continue
        if require_number and not NUMBER_RE.search(sentence):
            continue
        if len(sentence) > max_length:
            sentence = sentence[: max_length - 3].rstrip() + "..."
        matches.append(sentence)
    return _unique(matches)


def _species_mentions(text: str) -> list[str]:
    return _unique(f"{match.group(1)} {match.group(2)}" for match in SPECIES_RE.finditer(text))


def _strain_mentions(text: str) -> list[str]:
    matches: list[str] = []
    for pattern in STRAIN_PATTERNS:
        matches.extend(match.group(0) for match in pattern.finditer(text))
    return _unique(matches)


def _confidence(fields: Mapping[str, str]) -> float:
    score = 0.0
    if fields.get("species_mentioned"):
        score += 0.20
    if fields.get("strain_mentioned"):
        score += 0.10
    if fields.get("enzyme_terms"):
        score += 0.20
    if fields.get("product_terms"):
        score += 0.15
    if fields.get("mw_mentions") or fields.get("branching_mentions") or fields.get("nmr_mentions"):
        score += 0.25
    if fields.get("viscosity_mentions") or fields.get("yield_mentions"):
        score += 0.10
    return round(min(score, 1.0), 3)


def extract_mentions_from_text(
    title: str = "",
    abstract: str = "",
    record_id: str | int = "",
) -> dict[str, object]:
    text = _clean(f"{title}. {abstract}")
    fields: dict[str, str] = {
        "species_mentioned": _join(_species_mentions(text)),
        "strain_mentioned": _join(_strain_mentions(text)),
        "enzyme_terms": _join(_term_mentions(text, ENZYME_TERMS)),
        "product_terms": _join(_term_mentions(text, PRODUCT_TERMS)),
        "mw_mentions": _join(_keyword_sentences(text, MW_KEYWORDS, require_number=True)),
        "branching_mentions": _join(_keyword_sentences(text, BRANCHING_KEYWORDS)),
        "nmr_mentions": _join(_keyword_sentences(text, NMR_KEYWORDS)),
        "viscosity_mentions": _join(_keyword_sentences(text, VISCOSITY_KEYWORDS)),
        "yield_mentions": _join(_keyword_sentences(text, YIELD_KEYWORDS, require_number=True)),
    }
    confidence = _confidence(fields)
    numeric_text = " ".join(
        fields[key]
        for key in (
            "mw_mentions",
            "branching_mentions",
            "nmr_mentions",
            "viscosity_mentions",
            "yield_mentions",
        )
    )
    needs_manual_review = bool(NUMBER_RE.search(numeric_text)) or confidence < 0.5
    return {
        "record_id": str(record_id),
        "title": _clean(title),
        **fields,
        "confidence": confidence,
        "needs_manual_review": needs_manual_review,
    }
