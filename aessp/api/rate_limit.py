from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Callable


Clock = Callable[[], float]
Sleeper = Callable[[float], None]


@dataclass
class RateLimiter:
    """Small deterministic rate limiter.

    The first request is allowed immediately. Later calls wait until at least
    ``1 / requests_per_second`` seconds have elapsed since the previous slot.
    Injecting ``clock`` and ``sleep_func`` makes tests deterministic.
    """

    requests_per_second: float = 3.0
    sleep_func: Sleeper = time.sleep
    clock: Clock = time.monotonic
    enabled: bool = True
    _next_allowed_at: float | None = field(default=None, init=False, repr=False)

    def __post_init__(self) -> None:
        if self.requests_per_second <= 0:
            raise ValueError("requests_per_second must be positive")

    @property
    def min_interval(self) -> float:
        return 1.0 / self.requests_per_second

    def wait(self) -> float:
        if not self.enabled:
            return 0.0

        now = self.clock()
        if self._next_allowed_at is None:
            self._next_allowed_at = now + self.min_interval
            return 0.0

        delay = self._next_allowed_at - now
        if delay > 0:
            self.sleep_func(delay)
            self._next_allowed_at += self.min_interval
            return delay

        self._next_allowed_at = now + self.min_interval
        return 0.0
