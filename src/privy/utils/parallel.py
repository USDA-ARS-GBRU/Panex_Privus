"""Thread/process pool helpers for Panex Privus.

The scan backends are currently serial, but this helper keeps deterministic
ordering for code paths that opt into simple thread-pool mapping.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable
from concurrent.futures import ThreadPoolExecutor
from typing import TypeVar

T = TypeVar("T")
R = TypeVar("R")


def map_with_threads(
    fn: Callable[[T], R],
    items: Iterable[T],
    threads: int = 1,
) -> list[R]:
    """Apply *fn* to each item, optionally using a thread pool.

    Results preserve the input order, matching built-in ``map`` semantics.
    """
    if threads == 1:
        return [fn(item) for item in items]
    with ThreadPoolExecutor(max_workers=threads) as executor:
        return list(executor.map(fn, items))
