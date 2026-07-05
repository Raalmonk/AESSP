from __future__ import annotations

import responses

from aessp.api.cache import FileCache
from aessp.api.ncbi import BASE_URL, NCBIClient
from aessp.api.rate_limit import RateLimiter


def test_file_cache_ignores_api_key_in_cache_key(tmp_path):
    cache = FileCache(tmp_path)
    first = cache.path_for("ncbi", "esearch.fcgi", {"term": "dextran", "api_key": "SECRET"})
    second = cache.path_for("ncbi", "esearch.fcgi", {"term": "dextran", "api_key": "OTHER"})

    assert first == second
    assert "SECRET" not in str(first)

    cache.set_json("ncbi", "esearch.fcgi", {"term": "dextran", "api_key": "SECRET"}, {"ok": True})
    assert cache.get_json("ncbi", "esearch.fcgi", {"term": "dextran", "api_key": "OTHER"}) == {
        "ok": True
    }


def test_rate_limiter_is_deterministic_with_injected_clock():
    now = [0.0]
    slept: list[float] = []

    def clock() -> float:
        return now[0]

    def sleep(delay: float) -> None:
        slept.append(delay)
        now[0] += delay

    limiter = RateLimiter(requests_per_second=3.0, clock=clock, sleep_func=sleep)

    assert limiter.wait() == 0.0
    assert round(limiter.wait(), 3) == 0.333
    now[0] += 0.1
    assert round(limiter.wait(), 3) == 0.233
    assert len(slept) == 2


@responses.activate
def test_ncbi_client_uses_cache_for_mocked_request(tmp_path):
    cache = FileCache(tmp_path)
    client = NCBIClient(
        email="user@example.org",
        api_key="SECRET",
        cache=cache,
        rate_limiter=RateLimiter(enabled=False),
    )
    responses.add(
        responses.GET,
        f"{BASE_URL}/esearch.fcgi",
        json={"esearchresult": {"idlist": ["123"]}},
        status=200,
    )

    assert client.esearch("pubmed", "dextran")["esearchresult"]["idlist"] == ["123"]
    assert len(responses.calls) == 1
    assert "api_key=SECRET" in responses.calls[0].request.url

    responses.reset()
    assert client.esearch("pubmed", "dextran")["esearchresult"]["idlist"] == ["123"]
    assert len(responses.calls) == 0
