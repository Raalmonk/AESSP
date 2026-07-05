from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Mapping


SECRET_PARAM_RE = re.compile(
    r"(api[_-]?key|apikey|authorization|bearer|password|secret|token)",
    re.IGNORECASE,
)


def sanitize_params(params: Mapping[str, object] | None) -> dict[str, object]:
    cleaned: dict[str, object] = {}
    for key, value in (params or {}).items():
        if value is None or SECRET_PARAM_RE.search(str(key)):
            continue
        if isinstance(value, (list, tuple, set)):
            cleaned[str(key)] = [str(item) for item in value]
        else:
            cleaned[str(key)] = value
    return dict(sorted(cleaned.items()))


def cache_key(source: str, endpoint: str, params: Mapping[str, object] | None) -> str:
    payload = {
        "endpoint": endpoint,
        "params": sanitize_params(params),
        "source": source,
    }
    canonical = json.dumps(payload, ensure_ascii=True, sort_keys=True, default=str)
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def _safe_segment(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_") or "cache"


class FileCache:
    def __init__(self, root: str | Path = "data/cache", enabled: bool = True) -> None:
        self.root = Path(root)
        self.enabled = enabled

    def path_for(
        self,
        source: str,
        endpoint: str,
        params: Mapping[str, object] | None,
        suffix: str = ".json",
    ) -> Path:
        digest = cache_key(source, endpoint, params)
        return self.root / _safe_segment(source) / f"{digest}{suffix}"

    def get_json(
        self,
        source: str,
        endpoint: str,
        params: Mapping[str, object] | None,
    ) -> dict | list | None:
        if not self.enabled:
            return None
        path = self.path_for(source, endpoint, params, ".json")
        if not path.exists():
            return None
        with path.open("r", encoding="utf-8") as handle:
            return json.load(handle)

    def set_json(
        self,
        source: str,
        endpoint: str,
        params: Mapping[str, object] | None,
        payload: dict | list,
    ) -> Path | None:
        if not self.enabled:
            return None
        path = self.path_for(source, endpoint, params, ".json")
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8") as handle:
            json.dump(payload, handle, ensure_ascii=False, sort_keys=True)
        return path

    def get_text(
        self,
        source: str,
        endpoint: str,
        params: Mapping[str, object] | None,
    ) -> str | None:
        if not self.enabled:
            return None
        path = self.path_for(source, endpoint, params, ".txt")
        if not path.exists():
            return None
        return path.read_text(encoding="utf-8")

    def set_text(
        self,
        source: str,
        endpoint: str,
        params: Mapping[str, object] | None,
        payload: str,
    ) -> Path | None:
        if not self.enabled:
            return None
        path = self.path_for(source, endpoint, params, ".txt")
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(payload, encoding="utf-8")
        return path
