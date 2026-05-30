"""Collinearity chaining of alignment anchors (DAGchainer / MCScanX style).

Where :mod:`privy.synteny.graph_blocks` derives synteny from exact shared graph
segments, this module chains *alignment* anchors — approximate, gappy mappings
read from a PAF (``odgi untangle`` / ``minimap2`` / ``wfmash``) or a gene-pair
table — into collinear blocks.  It is the bridge for non-graph inputs and for
repeat-heavy loci where co-traversal is ambiguous.

The chainer is a pure-Python dynamic program over anchor points, in the lineage of
DAGchainer (Haas et al. 2004) and MCScanX (Wang et al. 2012): anchors that
preserve order on both axes are linked, rewarding matches and penalising the
indel implied by unequal query/target gaps; maximal high-scoring chains with
enough anchors become blocks (forward → collinear, reverse → inversion).

Default parameters echo MCScanX (``match_score=50``, ``gap`` penalised) but use
base-pair-aware gap costs (DAGchainer ``-D`` style) since anchors carry bp
coordinates.  Permutation E-values are not yet computed (score + anchor count are
reported); that refinement is deferred.

All coordinates are 0-based half-open.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path

from privy.io.paf import parse_paf
from privy.synteny.model import (
    Anchor,
    BlockType,
    GenomeInterval,
    SyntenyBlock,
)


@dataclass(frozen=True)
class ChainParams:
    """Tunables for collinearity chaining (MCScanX/DAGchainer-style)."""

    match_score: float = 50.0       # reward per chained anchor (MCScanX MATCH_SCORE)
    gap_open: float = 0.0           # fixed cost per inter-anchor link
    gap_extend: float = 1e-4        # cost per bp of indel (|Δquery − Δtarget|)
    max_gap: int = 200_000          # max bp between consecutive anchors per axis
    min_anchors: int = 3            # min anchors for a reported block (MCScanX MATCH_SIZE=5)
    min_score: float = 0.0          # min chain score to report


# ---------------------------------------------------------------------------
# Repeat suppression (optional pre-filter)
# ---------------------------------------------------------------------------


def suppress_repetitive_anchors(
    anchors: Sequence[Anchor],
    *,
    bin_size: int = 10_000,
    max_per_bin: int = 20,
) -> list[Anchor]:
    """Drop anchors landing in over-dense target bins (DEEPSPACE-style repeat masking).

    Counts anchors per ``(target_contig, target_start // bin_size)`` bin and removes
    every anchor in bins exceeding *max_per_bin* — a cheap guard against
    repeat-induced spurious anchors, important for plant genomes.
    """
    counts: dict[tuple[str, int], int] = {}
    for a in anchors:
        key = (a.target.contig, a.target.start // bin_size)
        counts[key] = counts.get(key, 0) + 1
    return [
        a for a in anchors
        if counts[(a.target.contig, a.target.start // bin_size)] <= max_per_bin
    ]


# ---------------------------------------------------------------------------
# Chaining
# ---------------------------------------------------------------------------


def chain_anchors(
    anchors: Iterable[Anchor],
    params: ChainParams | None = None,
) -> list[SyntenyBlock]:
    """Chain anchors into collinear/inverted blocks.

    Anchors are grouped by ``(query_genome, query_contig, target_genome,
    target_contig)`` and strand; each group is chained independently.  Returns
    blocks sorted by query start.
    """
    params = params or ChainParams()
    groups: dict[tuple[str, str, str, str, str], list[Anchor]] = {}
    for a in anchors:
        key = (a.query.genome, a.query.contig, a.target.genome, a.target.contig, a.strand)
        groups.setdefault(key, []).append(a)

    blocks: list[SyntenyBlock] = []
    block_idx = 0
    for key, group in groups.items():
        strand = key[4]
        for chain in _chain_one_group(group, strand, params):
            blocks.append(_block_from_chain(chain, strand, block_idx))
            block_idx += 1

    blocks.sort(key=lambda b: (b.query.contig, b.query.start, b.query.end))
    return blocks


def _chain_one_group(
    group: list[Anchor],
    strand: str,
    params: ChainParams,
) -> list[list[Anchor]]:
    """Run the chaining DP on one (contig-pair, strand) group; return anchor chains."""
    anchors = sorted(group, key=lambda a: (a.query.start, a.target.start))
    n = len(anchors)
    if n == 0:
        return []

    score = [params.match_score] * n
    pred = [-1] * n
    for i in range(n):
        qi, ti = anchors[i].query.start, anchors[i].target.start
        for j in range(i):
            dq = qi - anchors[j].query.start
            if dq <= 0:
                continue
            if strand == "+":
                dt = ti - anchors[j].target.start
            else:
                dt = anchors[j].target.start - ti
            if dt <= 0 or dq > params.max_gap or dt > params.max_gap:
                continue
            indel = abs(dq - dt)
            cand = score[j] + params.match_score - (params.gap_open + params.gap_extend * indel)
            if cand > score[i]:
                score[i] = cand
                pred[i] = j

    return _extract_chains(anchors, score, pred, params)


def _extract_chains(
    anchors: list[Anchor],
    score: list[float],
    pred: list[int],
    params: ChainParams,
) -> list[list[Anchor]]:
    """Greedily peel maximal high-scoring chains, highest score first."""
    order = sorted(range(len(anchors)), key=lambda i: score[i], reverse=True)
    used = [False] * len(anchors)
    chains: list[list[Anchor]] = []
    for endpoint in order:
        if used[endpoint]:
            continue
        idxs: list[int] = []
        i = endpoint
        chain_score = score[endpoint]
        while i != -1 and not used[i]:
            idxs.append(i)
            used[i] = True
            i = pred[i]
        idxs.reverse()
        if len(idxs) >= params.min_anchors and chain_score >= params.min_score:
            chains.append([anchors[k] for k in idxs])
    return chains


def _block_from_chain(chain: list[Anchor], strand: str, idx: int) -> SyntenyBlock:
    q_genome = chain[0].query.genome
    q_contig = chain[0].query.contig
    t_genome = chain[0].target.genome
    t_contig = chain[0].target.contig
    q_start = min(a.query.start for a in chain)
    q_end = max(a.query.end for a in chain)
    t_start = min(a.target.start for a in chain)
    t_end = max(a.target.end for a in chain)
    block_type = BlockType.COLLINEAR if strand == "+" else BlockType.INVERSION
    return SyntenyBlock(
        block_id=f"chain:B{idx}",
        query=GenomeInterval(q_genome, q_contig, q_start, q_end),
        target=GenomeInterval(t_genome, t_contig, t_start, t_end),
        strand=strand,
        block_type=block_type,
        anchors=tuple(chain),
        score=float(len(chain)),
    )


# ---------------------------------------------------------------------------
# PAF convenience
# ---------------------------------------------------------------------------


def chain_paf(
    paf_path: Path,
    params: ChainParams | None = None,
    *,
    pansn_delimiter: str = "#",
    suppress_repeats: bool = False,
    bin_size: int = 10_000,
    max_per_bin: int = 20,
) -> list[SyntenyBlock]:
    """Read a PAF, convert rows to anchors, optionally mask repeats, and chain.

    Args:
        suppress_repeats: Apply :func:`suppress_repetitive_anchors` before chaining.
    """
    anchors = [
        Anchor.from_paf(rec, pansn_delimiter=pansn_delimiter)
        for rec in parse_paf(Path(paf_path))
    ]
    if suppress_repeats:
        anchors = suppress_repetitive_anchors(
            anchors, bin_size=bin_size, max_per_bin=max_per_bin
        )
    return chain_anchors(anchors, params)
