"""Interval construction and merging for Panex Privus.

Single variants are often insufficient as biological objects.  This module
provides the logic to merge nearby passing loci into candidate genomic regions
for downstream scoring and reporting.

Coordinate convention: 0-based, half-open [start, end) — same as :class:`~privy.core.locus.Locus`.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator

from privy.core.locus import Locus, LocusType


def merge_loci_to_regions(
    loci: Iterable[Locus],
    merge_distance: int = 0,
    same_variant_class_only: bool = False,
    region_id_prefix: str = "REGION",
) -> list[Locus]:
    """Merge nearby passing loci into candidate genomic regions.

    Loci are sorted by ``(contig, start)`` internally.  Two loci are merged
    into the same region when they are on the same contig, their gap (in bp)
    is ≤ *merge_distance*, and optionally they share the same
    :class:`~privy.core.locus.LocusType`.

    Args:
        loci: Iterable of :class:`~privy.core.locus.Locus` objects.  Does not
            need to be pre-sorted.
        merge_distance: Maximum gap in bp between adjacent loci to merge.
            ``0`` means only overlapping or directly adjacent loci are merged.
        same_variant_class_only: If True, loci of different
            :class:`~privy.core.locus.LocusType` are never merged.
        region_id_prefix: Prefix for generated region IDs
            (e.g., ``"REGION"`` → ``"REGION000000"``).

    Returns:
        List of merged region :class:`~privy.core.locus.Locus` objects, sorted
        by ``(contig, start)``.  Single non-merged loci are also wrapped in a
        region-type Locus with a normalised ID.
    """
    sorted_loci: list[Locus] = sorted(loci, key=lambda loc: (loc.contig, loc.start))

    if not sorted_loci:
        return []

    regions: list[Locus] = []
    region_index = 0

    # Current merge group (accumulates loci that will become one region)
    current_group: list[Locus] = [sorted_loci[0]]

    for locus in sorted_loci[1:]:
        prev = current_group[-1]

        # ── Different contig: flush group ──────────────────────────────────
        if locus.contig != prev.contig:
            regions.append(_finalise_region(current_group, region_id_prefix, region_index))
            region_index += 1
            current_group = [locus]
            continue

        # ── Different variant class when restricted ────────────────────────
        if same_variant_class_only and locus.locus_type != current_group[0].locus_type:
            regions.append(_finalise_region(current_group, region_id_prefix, region_index))
            region_index += 1
            current_group = [locus]
            continue

        # ── Gap check ─────────────────────────────────────────────────────
        gap = locus.start - prev.end  # may be negative for overlapping loci
        if gap <= merge_distance:
            current_group.append(locus)
        else:
            regions.append(_finalise_region(current_group, region_id_prefix, region_index))
            region_index += 1
            current_group = [locus]

    # Flush the final group
    regions.append(_finalise_region(current_group, region_id_prefix, region_index))
    return regions


def _finalise_region(loci: list[Locus], prefix: str, index: int) -> Locus:
    """Collapse a list of loci into a single merged region Locus."""
    assert loci, "Cannot finalise an empty locus group."

    source_ids: list[str] = []
    for loc in loci:
        source_ids.extend(loc.source_ids)

    return Locus(
        locus_id=f"{prefix}{index:06d}",
        contig=loci[0].contig,
        start=loci[0].start,
        end=loci[-1].end,
        locus_type=LocusType.REGION,
        primary_source=loci[0].primary_source,
        source_ids=source_ids,
        metadata={"n_constituent_loci": str(len(loci))},
    )


def reciprocal_overlap(a: Locus, b: Locus) -> float:
    """Return the reciprocal overlap fraction between two loci.

    Reciprocal overlap = intersection / union of the two intervals.

    Args:
        a: First locus.
        b: Second locus.

    Returns:
        Overlap fraction in [0, 1].  Returns 0.0 if loci are on different
        contigs or do not overlap.
    """
    if a.contig != b.contig:
        return 0.0
    overlap_start = max(a.start, b.start)
    overlap_end = min(a.end, b.end)
    if overlap_end <= overlap_start:
        return 0.0
    intersection = overlap_end - overlap_start
    union = max(a.end, b.end) - min(a.start, b.start)
    return intersection / union if union > 0 else 0.0


def iter_contig_chunks(
    loci: list[Locus],
    chunk_size: int,
) -> Iterator[list[Locus]]:
    """Yield sub-lists of loci partitioned by contig, then by chunk.

    Used to feed the VCF scan backend one chunk at a time without holding
    all loci in memory.

    Args:
        loci: Loci sorted by ``(contig, start)``.
        chunk_size: Window size in bp.

    Yields:
        Lists of loci falling within the same contig chunk.
    """
    if not loci:
        return

    current_contig = loci[0].contig
    chunk_start = (loci[0].start // chunk_size) * chunk_size
    chunk: list[Locus] = []

    for locus in loci:
        if locus.contig != current_contig:
            if chunk:
                yield chunk
            current_contig = locus.contig
            chunk_start = (locus.start // chunk_size) * chunk_size
            chunk = [locus]
            continue

        chunk_index = (locus.start // chunk_size) * chunk_size
        if chunk_index != chunk_start:
            if chunk:
                yield chunk
            chunk_start = chunk_index
            chunk = []
        chunk.append(locus)

    if chunk:
        yield chunk


# ---------------------------------------------------------------------------
# Hierarchical base-assignment-by-class + run-length smoothing
# ---------------------------------------------------------------------------
#
# Adapted from geeViz's count_overlapsByGroup: assign every base in a window to
# at most one class by priority (no base counted twice), with the uncovered
# remainder collected as "missing".  This is the data behind stacked-area density
# / proportion tracks and presence/absence-by-class summaries.


def hierarchical_base_assignment(
    start: int,
    end: int,
    features: Iterable[tuple[int, int, str]],
    priority: list[str],
) -> dict[str, int]:
    """Assign each base in ``[start, end)`` to the highest-priority class covering it.

    Bases covered by several classes go to the earliest class in *priority*; bases
    covered by none accumulate in the ``"missing"`` bucket.  Returns base counts
    per class (every class in *priority* is present, even if 0) plus ``"missing"``;
    the counts sum to ``end - start``.

    Coordinate convention: 0-based half-open.  Exact (interval breakpoints, not
    per-base), so it scales to large windows.

    Args:
        features: ``(feature_start, feature_end, class)`` tuples (any coordinates;
            clipped to the window). Classes not in *priority* are ignored.
    """
    rank = {c: i for i, c in enumerate(priority)}
    clipped: list[tuple[int, int, str]] = []
    breakpoints = {start, end}
    for f_start, f_end, cls in features:
        if cls not in rank:
            continue
        lo = max(start, f_start)
        hi = min(end, f_end)
        if hi <= lo:
            continue
        clipped.append((lo, hi, cls))
        breakpoints.add(lo)
        breakpoints.add(hi)

    counts = {c: 0 for c in priority}
    counts["missing"] = 0
    points = sorted(p for p in breakpoints if start <= p <= end)
    for seg_start, seg_end in zip(points, points[1:], strict=False):
        if seg_end <= seg_start:
            continue
        best: str | None = None
        best_rank: int | None = None
        for lo, hi, cls in clipped:
            if lo <= seg_start and seg_end <= hi and (best_rank is None or rank[cls] < best_rank):
                best_rank = rank[cls]
                best = cls
        width = seg_end - seg_start
        if best is None:
            counts["missing"] += width
        else:
            counts[best] += width
    return counts


def class_proportions(counts: dict[str, int]) -> dict[str, float]:
    """Normalise per-class base counts to proportions summing to 1 (0 if empty)."""
    total = sum(counts.values())
    if total == 0:
        return {c: 0.0 for c in counts}
    return {c: n / total for c, n in counts.items()}


def run_length_smooth(labels: list[str], min_run: int) -> list[str]:
    """Relabel runs shorter than *min_run* to a neighbouring run's label.

    Removes salt-and-pepper noise from per-window classifications (geeViz's
    run-length filter): a short run adopts the previous surviving run's label, or
    the next run's label if it is the leading run.  Returns a new list; a no-op
    when *min_run* <= 1 or every run is short.
    """
    if min_run <= 1 or not labels:
        return list(labels)

    # Identify maximal runs as (label, start_index, length).
    runs: list[tuple[str, int, int]] = []
    i = 0
    n = len(labels)
    while i < n:
        j = i
        while j < n and labels[j] == labels[i]:
            j += 1
        runs.append((labels[i], i, j - i))
        i = j
    if all(length < min_run for _label, _start, length in runs):
        return list(labels)

    out = list(labels)
    for idx, (_label, start_i, length) in enumerate(runs):
        if length >= min_run:
            continue
        # Find the nearest long run: prefer previous, else next.
        replacement: str | None = None
        for j in range(idx - 1, -1, -1):
            if runs[j][2] >= min_run:
                replacement = runs[j][0]
                break
        if replacement is None:
            for j in range(idx + 1, len(runs)):
                if runs[j][2] >= min_run:
                    replacement = runs[j][0]
                    break
        if replacement is None:
            continue
        for k in range(start_i, start_i + length):
            out[k] = replacement
    return out
