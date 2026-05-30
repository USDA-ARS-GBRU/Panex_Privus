"""Per-path coordinate model for pangenome graphs.

A pangenome graph embeds many genomes as *paths* (P-lines) or *walks* (W-lines)
through a shared set of segments.  Each base therefore has two identities:

* an **unstable** graph identity ``(segment, offset)`` — changes if the graph is
  rebuilt or re-sorted;
* a **stable** reference identity ``(contig, coordinate)`` — survives graph edits.

:class:`PathCoordinateModel` precomputes, for every embedded path, the cumulative
path-local offset of each step, so it can answer in O(log n):

* where a segment occurs on a path (path-local span(s)), and
* what segment + within-segment offset a path-local position falls in,

and map either of those to a **stable** ``(contig, coordinate)``.  This is the
arithmetic substrate for projecting a region to *any* reference in the graph
(the cross-path projection itself lands in ``synteny/projection.py``).

Coordinate convention: 0-based, half-open ``[start, end)`` throughout — matching
:class:`~privy.core.locus.Locus` and the GFA W-line ``seq_start``/``seq_end``.

Assumptions / caveats
---------------------
* Path-local coordinate 0 corresponds to stable position ``base_offset``:
  ``seq_start`` for W-lines, and ``0`` for P-lines (whose name's final
  ``#``-component is taken as the stable contig).  This holds for pangenome
  P-lines that represent a whole sequence tiled by non-overlapping segments
  (PGGB / minigraph-cactus output).  rGFA ``SO``-based refinement is deferred.
* A segment may occur more than once on a path (tandem/CNV/inversion); all
  occurrences are indexed.
"""

from __future__ import annotations

import logging
from bisect import bisect_right
from dataclasses import dataclass, field

from privy.io.gfa import GfaGraph

log = logging.getLogger("privy.synteny.coordinates")


# ---------------------------------------------------------------------------
# Result objects
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class SegmentOccurrence:
    """One occurrence of a segment on a path, in path-local coordinates."""

    step_index: int
    start: int          # path-local, inclusive
    end: int            # path-local, exclusive
    orientation: str    # "+" or "-"


@dataclass(frozen=True)
class PathStepLocation:
    """Resolution of a path-local position to a segment + within-segment offset."""

    path_id: str
    path_position: int       # the queried path-local position
    segment: str
    orientation: str         # "+" or "-"
    step_index: int
    offset_in_segment: int   # 0-based offset within the segment's sequence


# ---------------------------------------------------------------------------
# Internal per-path record
# ---------------------------------------------------------------------------


@dataclass
class _PathIndex:
    steps: list[tuple[str, str]]          # (segment, orientation) in order
    starts: list[int]                     # path-local cumulative start per step
    total_length: int                     # full path length in bp
    contig: str                           # stable contig name
    base_offset: int                      # stable coordinate of path position 0
    seg_occurrences: dict[str, list[int]] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# PathCoordinateModel
# ---------------------------------------------------------------------------


class PathCoordinateModel:
    """Cumulative-offset index over every embedded path/walk of a GFA graph.

    Build with :meth:`from_graph`.  All public lookups are keyed by a unique
    *path id*: the P-line name for paths, and ``sample#hap#seqid`` (de-duplicated
    with an ``@k`` suffix on collision) for walks.
    """

    def __init__(self, indices: dict[str, _PathIndex]) -> None:
        self._paths = indices

    # -- construction ------------------------------------------------------

    @classmethod
    def from_graph(cls, graph: GfaGraph, *, strict: bool = True) -> PathCoordinateModel:
        """Build the model from a parsed :class:`~privy.io.gfa.GfaGraph`.

        Args:
            graph: A graph parsed by :func:`~privy.io.gfa.parse_gfa`.
            strict: If True (default), raise when a path references a segment
                absent from the graph.  If False, skip that path with a warning.
        """
        seg_len = {name: seg.length for name, seg in graph.segments.items()}
        indices: dict[str, _PathIndex] = {}

        # P-lines (classic paths).
        for name, path in graph.paths.items():
            steps = list(zip(path.segment_names, path.orientations, strict=True))
            built = cls._build_index(
                name, steps, seg_len, contig=_contig_from_path_name(name),
                base_offset=0, strict=strict,
            )
            if built is not None:
                indices[name] = built

        # W-lines (walks); construct unique ids.
        for walk in graph.walks:
            steps = [(step.segment, step.orient) for step in walk.steps]
            walk_id = _unique_id(f"{walk.sample}#{walk.hap_index}#{walk.seq_id}", indices)
            built = cls._build_index(
                walk_id, steps, seg_len, contig=walk.seq_id,
                base_offset=walk.seq_start, strict=strict,
            )
            if built is not None:
                indices[walk_id] = built

        return cls(indices)

    @staticmethod
    def _build_index(
        path_id: str,
        steps: list[tuple[str, str]],
        seg_len: dict[str, int],
        *,
        contig: str,
        base_offset: int,
        strict: bool,
    ) -> _PathIndex | None:
        starts: list[int] = []
        seg_occurrences: dict[str, list[int]] = {}
        cursor = 0
        for step_idx, (segment, _orient) in enumerate(steps):
            length = seg_len.get(segment)
            if length is None:
                msg = f"Path {path_id!r} references unknown segment {segment!r}."
                if strict:
                    raise KeyError(msg)
                log.warning("%s Skipping path.", msg)
                return None
            starts.append(cursor)
            seg_occurrences.setdefault(segment, []).append(step_idx)
            cursor += length
        return _PathIndex(
            steps=steps,
            starts=starts,
            total_length=cursor,
            contig=contig,
            base_offset=base_offset,
            seg_occurrences=seg_occurrences,
        )

    # -- container niceties ------------------------------------------------

    def __contains__(self, path_id: str) -> bool:
        return path_id in self._paths

    def __len__(self) -> int:
        return len(self._paths)

    def path_ids(self) -> list[str]:
        """Return all path ids in deterministic (insertion) order."""
        return list(self._paths)

    def _require(self, path_id: str) -> _PathIndex:
        index = self._paths.get(path_id)
        if index is None:
            raise KeyError(f"Unknown path id {path_id!r}.")
        return index

    # -- queries -----------------------------------------------------------

    def path_length(self, path_id: str) -> int:
        """Total length of *path_id* in bp."""
        return self._require(path_id).total_length

    def stable_contig(self, path_id: str) -> str:
        """Stable contig name that *path_id*'s coordinates project onto."""
        return self._require(path_id).contig

    def occurrences(self, path_id: str, segment: str) -> list[SegmentOccurrence]:
        """All path-local spans where *segment* occurs on *path_id*.

        Returns an empty list when the segment is not on the path.  Multiple
        results indicate tandem copies / CNV / a revisited node.
        """
        index = self._require(path_id)
        result: list[SegmentOccurrence] = []
        for step_idx in index.seg_occurrences.get(segment, ()):
            start = index.starts[step_idx]
            end = (
                index.starts[step_idx + 1]
                if step_idx + 1 < len(index.starts)
                else index.total_length
            )
            result.append(
                SegmentOccurrence(
                    step_index=step_idx,
                    start=start,
                    end=end,
                    orientation=index.steps[step_idx][1],
                )
            )
        return result

    def locate(self, path_id: str, position: int) -> PathStepLocation:
        """Resolve a path-local *position* to its segment + within-segment offset.

        Args:
            path_id: Path identifier.
            position: 0-based path-local coordinate, ``0 <= position < length``.

        Raises:
            IndexError: If *position* is outside ``[0, path_length)``.
        """
        index = self._require(path_id)
        if position < 0 or position >= index.total_length:
            raise IndexError(
                f"position {position} out of range [0, {index.total_length}) "
                f"for path {path_id!r}."
            )
        step_idx = bisect_right(index.starts, position) - 1
        segment, orient = index.steps[step_idx]
        return PathStepLocation(
            path_id=path_id,
            path_position=position,
            segment=segment,
            orientation=orient,
            step_index=step_idx,
            offset_in_segment=position - index.starts[step_idx],
        )

    def to_stable(self, path_id: str, position: int) -> tuple[str, int]:
        """Map a path-local *position* to a stable ``(contig, coordinate)``.

        Raises:
            IndexError: If *position* is outside ``[0, path_length)``.
        """
        index = self._require(path_id)
        if position < 0 or position >= index.total_length:
            raise IndexError(
                f"position {position} out of range [0, {index.total_length}) "
                f"for path {path_id!r}."
            )
        return index.contig, index.base_offset + position


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _contig_from_path_name(name: str) -> str:
    """Stable contig for a P-line: the final ``#``-delimited component (PanSN)."""
    return name.rsplit("#", 1)[-1] if "#" in name else name


def _unique_id(candidate: str, existing: dict[str, _PathIndex]) -> str:
    """Return *candidate*, or ``candidate@k`` if it collides with an existing id."""
    if candidate not in existing:
        return candidate
    k = 1
    while f"{candidate}@{k}" in existing:
        k += 1
    return f"{candidate}@{k}"
