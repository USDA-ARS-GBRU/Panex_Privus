"""Graph-native synteny: collinear blocks + typed rearrangements from a graph.

Two genomes embedded in a pangenome graph are syntenic over a region exactly when
their paths traverse the same ordered run of segments.  This module derives those
blocks directly from path co-traversal (no realignment) and classifies each by
*conformation*, SyRI-style:

* **collinear** — shared segments in the same order and orientation;
* **inversion** — a maximal run traversed in the opposite order *and* orientation;
* **translocation** — a collinear block positioned out of order on the reference
  (overlaps an earlier block's reference span);
* **duplication** — a segment that occurs more than once (CNV) on the query (or
  reference) path.

The classification is a pragmatic graph-native adaptation of SyRI's
"longest-collinear-path, then type the residual by conformation"
(Goel et al. 2019, Genome Biology) operating over graph walks rather than a
pairwise whole-genome alignment.

All coordinates are 0-based half-open.  Builds on
:class:`~privy.synteny.coordinates.PathCoordinateModel`.
"""

from __future__ import annotations

from dataclasses import dataclass

from privy.synteny.coordinates import PathCoordinateModel
from privy.synteny.model import (
    Anchor,
    AnchorSource,
    BlockType,
    GenomeInterval,
    SyntenyBlock,
    split_pansn,
)


@dataclass(frozen=True)
class _AnchorRec:
    """Internal: one shared-segment anchor between a query and reference path."""

    segment: str
    q_start: int
    q_end: int
    r_start: int
    r_end: int
    rel: str   # relative orientation: "+" (same) or "-" (flipped)


def build_pairwise_blocks(
    model: PathCoordinateModel,
    query_path: str,
    ref_path: str,
    *,
    min_block_anchors: int = 1,
) -> list[SyntenyBlock]:
    """Detect typed synteny blocks of *query_path* against *ref_path*.

    Returns blocks sorted by query start.  Segments present only on the query
    (insertions) are skipped; segments occurring multiple times become
    DUPLICATION blocks; the remaining backbone is chained into COLLINEAR /
    INVERSION runs, then collinear runs that sit out of reference order are
    re-typed as TRANSLOCATION.

    Args:
        min_block_anchors: Drop collinear/inversion blocks supported by fewer
            than this many anchors (duplications are always reported).
    """
    backbone: list[_AnchorRec] = []
    duplications: list[_AnchorRec] = []
    seen: set[str] = set()

    for step in model.iter_steps(query_path):
        ref_occs = model.occurrences(ref_path, step.segment)
        if not ref_occs:
            continue  # query-private segment (insertion)
        rocc = ref_occs[0]
        rel = "+" if step.orientation == rocc.orientation else "-"
        rec = _AnchorRec(
            segment=step.segment,
            q_start=step.start, q_end=step.end,
            r_start=rocc.start, r_end=rocc.end,
            rel=rel,
        )
        is_dup = len(ref_occs) > 1 or len(model.occurrences(query_path, step.segment)) > 1
        if is_dup and step.segment in seen:
            duplications.append(rec)
        else:
            seen.add(step.segment)
            backbone.append(rec)

    runs = _chain(backbone)
    blocks = [
        _finalize_block(model, query_path, ref_path, run, idx)
        for idx, run in enumerate(runs)
    ]
    blocks = _reclassify_translocations(blocks)
    if min_block_anchors > 1:
        blocks = [b for b in blocks if b.n_anchors >= min_block_anchors]

    next_id = len(blocks)
    for offset, rec in enumerate(duplications):
        blocks.append(
            _finalize_block(model, query_path, ref_path, [rec], next_id + offset,
                            forced_type=BlockType.DUPLICATION)
        )

    blocks.sort(key=lambda b: (b.query.start, b.query.end))
    return blocks


def _chain(anchors: list[_AnchorRec]) -> list[list[_AnchorRec]]:
    """Group anchors (in query order) into maximal consistent collinear/inverted runs."""
    runs: list[list[_AnchorRec]] = []
    current: list[_AnchorRec] = []
    for a in anchors:
        if not current:
            current = [a]
            continue
        prev = current[-1]
        same_rel = a.rel == current[0].rel
        increasing = a.r_start > prev.r_start
        decreasing = a.r_start < prev.r_start
        continues = same_rel and (increasing if current[0].rel == "+" else decreasing)
        if continues:
            current.append(a)
        else:
            runs.append(current)
            current = [a]
    if current:
        runs.append(current)
    return runs


def _reclassify_translocations(blocks: list[SyntenyBlock]) -> list[SyntenyBlock]:
    """Re-type collinear blocks that sit out of reference order as TRANSLOCATION.

    Walks blocks in query order; a COLLINEAR block whose reference start falls
    before the running maximum reference end of the accepted backbone is moved
    (a translocation).  Inversions/duplications are left untouched.
    """
    running_max_ref_end = -1
    out: list[SyntenyBlock] = []
    for block in blocks:
        if block.block_type is BlockType.COLLINEAR:
            if block.target.start < running_max_ref_end:
                out.append(_retype(block, BlockType.TRANSLOCATION))
                continue
            running_max_ref_end = max(running_max_ref_end, block.target.end)
        out.append(block)
    return out


def _retype(block: SyntenyBlock, block_type: BlockType) -> SyntenyBlock:
    return SyntenyBlock(
        block_id=block.block_id,
        query=block.query,
        target=block.target,
        strand=block.strand,
        block_type=block_type,
        anchors=block.anchors,
        score=block.score,
        e_value=block.e_value,
    )


def _finalize_block(
    model: PathCoordinateModel,
    query_path: str,
    ref_path: str,
    run: list[_AnchorRec],
    idx: int,
    *,
    forced_type: BlockType | None = None,
) -> SyntenyBlock:
    """Collapse a run of anchors into a typed :class:`SyntenyBlock`."""
    q_lo = min(a.q_start for a in run)
    q_hi = max(a.q_end for a in run)
    r_lo = min(a.r_start for a in run)
    r_hi = max(a.r_end for a in run)
    rel = run[0].rel

    q_contig, q_start = model.to_stable(query_path, q_lo)
    _, q_last = model.to_stable(query_path, q_hi - 1)
    r_contig, r_start = model.to_stable(ref_path, r_lo)
    _, r_last = model.to_stable(ref_path, r_hi - 1)

    query_iv = GenomeInterval(split_pansn(query_path)[0], q_contig, q_start, q_last + 1)
    target_iv = GenomeInterval(split_pansn(ref_path)[0], r_contig, r_start, r_last + 1)

    if forced_type is not None:
        block_type = forced_type
    elif rel == "-":
        block_type = BlockType.INVERSION
    else:
        block_type = BlockType.COLLINEAR

    anchors = tuple(
        _anchor_for(model, query_path, ref_path, rec) for rec in run
    )
    return SyntenyBlock(
        block_id=f"{query_path}~{ref_path}:B{idx}",
        query=query_iv,
        target=target_iv,
        strand=rel,
        block_type=block_type,
        anchors=anchors,
        score=float(len(run)),
    )


def _anchor_for(
    model: PathCoordinateModel,
    query_path: str,
    ref_path: str,
    rec: _AnchorRec,
) -> Anchor:
    q_contig, q_start = model.to_stable(query_path, rec.q_start)
    _, q_last = model.to_stable(query_path, rec.q_end - 1)
    r_contig, r_start = model.to_stable(ref_path, rec.r_start)
    _, r_last = model.to_stable(ref_path, rec.r_end - 1)
    return Anchor(
        query=GenomeInterval(split_pansn(query_path)[0], q_contig, q_start, q_last + 1),
        target=GenomeInterval(split_pansn(ref_path)[0], r_contig, r_start, r_last + 1),
        strand=rec.rel,
        source=AnchorSource.GRAPH,
        name=rec.segment,
    )
