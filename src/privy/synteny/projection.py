"""Coordinate projection across a pangenome graph — the flagship capability.

Given a position or region anchored to one genome (or defined directly in graph
node space), project it onto *any* other genome embedded in the graph.  This is
the mechanism behind "define a region once, see where it lands in every
reference" and behind lifting annotation tracks between assemblies.

The graph itself is the lift-over: two genomes share a segment exactly when they
contain the same sequence, so a base on a shared segment maps deterministically
between any two paths that traverse it — with orientation handled correctly so an
inverted traversal maps to the right base.

Two granularities are provided:

* :func:`project_coordinate` — a single base, orientation-aware, returning one
  candidate (or several, for a segment duplicated on the target → ``AMBIGUOUS``;
  none → ``ABSENT``).
* :func:`project_region` / :func:`project_node_set` — an interval (or a raw set
  of segments) projected onto many targets at once, returning a
  :class:`~privy.synteny.model.ProjectionMap`.

All coordinates are 0-based half-open.  Builds only on
:class:`~privy.synteny.coordinates.PathCoordinateModel`.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from enum import Enum

from privy.synteny.coordinates import PathCoordinateModel
from privy.synteny.model import GenomeInterval, ProjectionMap, split_pansn


class ProjectionStatus(str, Enum):
    """Outcome of projecting a single coordinate to a target path."""

    MAPPED = "mapped"        # exactly one target position
    AMBIGUOUS = "ambiguous"  # >1 target position (segment duplicated on target)
    ABSENT = "absent"        # segment not present on the target path


@dataclass(frozen=True)
class CoordinateProjection:
    """Result of projecting one source base onto a target path."""

    source_path: str
    source_position: int
    target_path: str
    segment: str
    status: ProjectionStatus
    # Stable (contig, coordinate) candidates on the target, in path order.
    targets: tuple[tuple[str, int], ...] = ()

    @property
    def best(self) -> tuple[str, int] | None:
        """First stable candidate, or None when ABSENT."""
        return self.targets[0] if self.targets else None


# ---------------------------------------------------------------------------
# Single-coordinate projection
# ---------------------------------------------------------------------------


def project_coordinate(
    model: PathCoordinateModel,
    source_path: str,
    position: int,
    target_path: str,
) -> CoordinateProjection:
    """Project path-local *position* on *source_path* onto *target_path*.

    Orientation-aware: if the shared segment is traversed in opposite directions
    on the two paths, the within-segment offset is flipped so the *same physical
    base* is returned.

    Raises:
        IndexError: If *position* is out of range on *source_path*.
        KeyError: If either path id is unknown.
    """
    loc = model.locate(source_path, position)
    seg_len = model.segment_length(loc.segment)
    if seg_len is None:  # pragma: no cover - guarded by model construction
        raise KeyError(f"unknown segment {loc.segment!r}")

    # Index of the physical base within the segment's own (forward) sequence.
    base_index = (
        loc.offset_in_segment
        if loc.orientation == "+"
        else seg_len - 1 - loc.offset_in_segment
    )

    targets: list[tuple[str, int]] = []
    for occ in model.occurrences(target_path, loc.segment):
        offset_on_target = base_index if occ.orientation == "+" else seg_len - 1 - base_index
        target_local = occ.start + offset_on_target
        targets.append(model.to_stable(target_path, target_local))

    if not targets:
        status = ProjectionStatus.ABSENT
    elif len(targets) == 1:
        status = ProjectionStatus.MAPPED
    else:
        status = ProjectionStatus.AMBIGUOUS

    return CoordinateProjection(
        source_path=source_path,
        source_position=position,
        target_path=target_path,
        segment=loc.segment,
        status=status,
        targets=tuple(targets),
    )


# ---------------------------------------------------------------------------
# Region / node-set projection
# ---------------------------------------------------------------------------


def _project_segments_to_path(
    model: PathCoordinateModel,
    segments: Iterable[str],
    target_path: str,
) -> GenomeInterval | None:
    """Bounding stable interval of *segments* on *target_path*, or None if absent.

    Takes the min start / max end over every occurrence of every requested
    segment, so reordered or inverted segments still yield the span they cover.
    """
    lo: int | None = None
    hi: int | None = None
    for segment in segments:
        for occ in model.occurrences(target_path, segment):
            lo = occ.start if lo is None else min(lo, occ.start)
            hi = occ.end if hi is None else max(hi, occ.end)
    if lo is None or hi is None:
        return None
    contig, stable_start = model.to_stable(target_path, lo)
    _, stable_last = model.to_stable(target_path, hi - 1)
    return GenomeInterval(
        genome=split_pansn(target_path)[0],
        contig=contig,
        start=stable_start,
        end=stable_last + 1,
    )


def project_node_set(
    model: PathCoordinateModel,
    segments: Iterable[str],
    *,
    targets: Iterable[str] | None = None,
    source_label: str = "node-set",
) -> ProjectionMap:
    """Project a set of graph segments onto many paths at once.

    This is the graph-native flagship: define a region *once* in node space, then
    fan out to every embedded reference simultaneously.

    Args:
        model: The coordinate model.
        segments: Segment ids defining the region in node space.
        targets: Paths to project onto (default: all paths in the model).
        source_label: Label recorded as the projection source.
    """
    segments = list(segments)
    target_ids = list(targets) if targets is not None else model.path_ids()
    projections: dict[str, GenomeInterval | None] = {
        tp: _project_segments_to_path(model, segments, tp) for tp in target_ids
    }
    return ProjectionMap(source=source_label, projections=projections)


def project_region(
    model: PathCoordinateModel,
    source_path: str,
    start: int,
    end: int,
    *,
    targets: Iterable[str] | None = None,
    include_source: bool = True,
) -> ProjectionMap:
    """Project a path-local interval ``[start, end)`` of *source_path* onto targets.

    The interval is resolved to the segments it covers on the source path, then
    those segments are located on each target (see :func:`project_node_set`).

    Args:
        targets: Paths to project onto (default: all paths except the source,
            or including it when *include_source* is True).
        include_source: If True (default) and *targets* is None, the source path
            is included in the output (it projects to itself).

    Raises:
        KeyError: If *source_path* is unknown.
    """
    covered = [step.segment for step in model.segments_in_range(source_path, start, end)]
    if targets is None:
        target_ids = [
            tp for tp in model.path_ids() if include_source or tp != source_path
        ]
    else:
        target_ids = list(targets)
    projections: dict[str, GenomeInterval | None] = {
        tp: _project_segments_to_path(model, covered, tp) for tp in target_ids
    }
    return ProjectionMap(
        source=f"{source_path}:{start}-{end}",
        projections=projections,
    )


# ---------------------------------------------------------------------------
# Annotation-track liftover
# ---------------------------------------------------------------------------


def lift_intervals(
    model: PathCoordinateModel,
    intervals: Iterable[GenomeInterval],
    source_path: str,
    target_path: str,
) -> list[GenomeInterval | None]:
    """Lift annotation intervals from *source_path* onto *target_path*.

    Each input interval is given in *source_path*'s stable contig coordinates
    (e.g. parsed from a GFF/BED on that genome).  Returns one projected
    :class:`~privy.synteny.model.GenomeInterval` per input — or ``None`` where the
    feature does not project (out of range on the source, or absent on the
    target).  This is the engine behind pairing projected regions with gene-model
    / annotation tracks.
    """
    results: list[GenomeInterval | None] = []
    for iv in intervals:
        if iv.end <= iv.start:
            results.append(None)
            continue
        try:
            local_start = model.to_path_local(source_path, iv.start)
            local_end = model.to_path_local(source_path, iv.end - 1) + 1
        except (IndexError, KeyError):
            results.append(None)
            continue
        pm = project_region(model, source_path, local_start, local_end, targets=[target_path])
        results.append(pm.projections[target_path])
    return results
