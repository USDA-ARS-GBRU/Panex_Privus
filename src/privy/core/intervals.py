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
