"""Thread/process pool helpers for Panex Privus.

TODO (Phase 2): implement safe parallelisation helpers for contig-chunked
processing.  The design constraint is that parallelisation must not break
deterministic output ordering.
"""

from __future__ import annotations

from typing import Any, Callable, Iterable, TypeVar

T = TypeVar("T")
R = TypeVar("R")


def map_with_threads(
    fn: Callable[[T], R],
    items: Iterable[T],
    threads: int = 1,
) -> list[R]:
    """Apply *fn* to each item, optionally using a thread pool.

    Falls back to serial execution when *threads* == 1 (the default).

    TODO (Phase 2): implement concurrent.futures.ThreadPoolExecutor path.
    """
    if threads == 1:
        return [fn(item) for item in items]
    raise NotImplementedError(
        "Parallel execution (threads > 1) is not yet implemented.  "
        "Use --threads 1 (the default) for now."
    )
