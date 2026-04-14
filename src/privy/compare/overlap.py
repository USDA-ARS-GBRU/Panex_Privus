"""Interval overlap utilities for privy compare.

Re-exports :func:`~privy.core.intervals.reciprocal_overlap` and adds
source-aware tolerance wrappers.

TODO (Phase 5): implement :func:`find_overlapping_loci` index.
"""

from __future__ import annotations

from privy.core.intervals import reciprocal_overlap  # noqa: F401
from privy.core.locus import Locus


def find_overlapping_loci(
    query: Locus,
    candidates: list[Locus],
    min_reciprocal_overlap: float = 0.0,
) -> list[tuple[Locus, float]]:
    """Return all candidates overlapping *query* above *min_reciprocal_overlap*.

    Args:
        query: Query locus.
        candidates: Pool of candidate loci to search.
        min_reciprocal_overlap: Minimum reciprocal overlap fraction [0, 1].

    Returns:
        List of ``(locus, overlap_fraction)`` tuples, sorted by overlap desc.

    TODO (Phase 5): add an interval tree index for large candidate sets.
    """
    results: list[tuple[Locus, float]] = []
    for candidate in candidates:
        ro = reciprocal_overlap(query, candidate)
        if ro >= min_reciprocal_overlap:
            results.append((candidate, ro))
    return sorted(results, key=lambda t: t[1], reverse=True)
