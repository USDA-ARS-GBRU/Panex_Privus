"""Unit tests for src/privy/synteny/chain.py."""

from __future__ import annotations

from pathlib import Path

from privy.io.paf import PafRecord, write_paf
from privy.synteny.chain import (
    ChainParams,
    chain_anchors,
    chain_paf,
    suppress_repetitive_anchors,
)
from privy.synteny.model import Anchor, BlockType


def _anchor(q0, q1, t0, t1, strand="+", q="qA#0#chr1", t="tB#0#chr1"):
    rec = PafRecord(q, 100_000, q0, q1, strand, t, 100_000, t0, t1, q1 - q0, q1 - q0, 60)
    return Anchor.from_paf(rec)


# ---------------------------------------------------------------------------
# Forward / collinear chaining
# ---------------------------------------------------------------------------


class TestForwardChain:
    def test_single_collinear_chain(self):
        anchors = [_anchor(i * 100, i * 100 + 80, i * 100, i * 100 + 80) for i in range(5)]
        blocks = chain_anchors(anchors, ChainParams(min_anchors=3))
        assert len(blocks) == 1
        b = blocks[0]
        assert b.block_type is BlockType.COLLINEAR
        assert b.strand == "+"
        assert b.n_anchors == 5
        assert (b.query.start, b.query.end) == (0, 480)

    def test_min_anchors_filters_short_chains(self):
        anchors = [_anchor(i * 100, i * 100 + 80, i * 100, i * 100 + 80) for i in range(2)]
        assert chain_anchors(anchors, ChainParams(min_anchors=3)) == []

    def test_large_gap_breaks_into_two_chains(self):
        near = [_anchor(i * 100, i * 100 + 80, i * 100, i * 100 + 80) for i in range(4)]
        far = [
            _anchor(500_000 + i * 100, 500_000 + i * 100 + 80,
                    500_000 + i * 100, 500_000 + i * 100 + 80)
            for i in range(4)
        ]
        blocks = chain_anchors(near + far, ChainParams(min_anchors=3, max_gap=200_000))
        assert len(blocks) == 2   # gap of ~500kb > max_gap -> two separate blocks


# ---------------------------------------------------------------------------
# Reverse / inversion chaining
# ---------------------------------------------------------------------------


class TestReverseChain:
    def test_reverse_anchors_form_inversion(self):
        # query increasing, target decreasing, strand "-"
        anchors = [
            _anchor(i * 100, i * 100 + 80, (4 - i) * 100, (4 - i) * 100 + 80, strand="-")
            for i in range(5)
        ]
        blocks = chain_anchors(anchors, ChainParams(min_anchors=3))
        assert len(blocks) == 1
        assert blocks[0].block_type is BlockType.INVERSION
        assert blocks[0].strand == "-"
        assert blocks[0].n_anchors == 5


# ---------------------------------------------------------------------------
# Grouping by contig pair
# ---------------------------------------------------------------------------


class TestGrouping:
    def test_separate_contig_pairs_chain_independently(self):
        a = [_anchor(i * 100, i * 100 + 80, i * 100, i * 100 + 80, t="tB#0#chr1") for i in range(4)]
        b = [_anchor(i * 100, i * 100 + 80, i * 100, i * 100 + 80, t="tB#0#chr2") for i in range(4)]
        blocks = chain_anchors(a + b, ChainParams(min_anchors=3))
        contigs = sorted(blk.target.contig for blk in blocks)
        assert contigs == ["chr1", "chr2"]


# ---------------------------------------------------------------------------
# Repeat suppression
# ---------------------------------------------------------------------------


class TestRepeatSuppression:
    def test_dense_target_bin_dropped(self):
        # 30 anchors piled into one 10kb target bin + 3 spread out
        dense = [_anchor(i, i + 5, 0, 5) for i in range(30)]
        spread = [_anchor(1000 + i * 100, 1000 + i * 100 + 80, 50_000 + i * 100,
                          50_000 + i * 100 + 80) for i in range(3)]
        kept = suppress_repetitive_anchors(dense + spread, bin_size=10_000, max_per_bin=20)
        # the dense bin (>20) is dropped; the spread anchors survive
        assert len(kept) == 3
        assert all(a.target.start >= 50_000 for a in kept)


# ---------------------------------------------------------------------------
# PAF convenience
# ---------------------------------------------------------------------------


class TestChainPaf:
    def test_chain_from_paf_file(self, tmp_path: Path):
        recs = [
            PafRecord("qA#0#chr1", 100_000, i * 100, i * 100 + 80, "+",
                      "tB#0#chr1", 100_000, i * 100, i * 100 + 80, 80, 80, 60)
            for i in range(5)
        ]
        paf = tmp_path / "aln.paf"
        write_paf(recs, paf)
        blocks = chain_paf(paf, ChainParams(min_anchors=3))
        assert len(blocks) == 1
        assert blocks[0].n_anchors == 5
        assert blocks[0].query.genome == "qA"
        assert blocks[0].target.genome == "tB"
