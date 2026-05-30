"""Unit tests for src/privy/synteny/graph_blocks.py."""

from __future__ import annotations

from pathlib import Path

from privy.io.gfa import parse_gfa
from privy.synteny.coordinates import PathCoordinateModel
from privy.synteny.graph_blocks import build_pairwise_blocks
from privy.synteny.model import BlockType
from privy.synthetic import (
    collinear_pangenome,
    duplication_pangenome,
    inversion_pangenome,
    translocation_pangenome,
)


def _model(pg, tmp_path: Path) -> PathCoordinateModel:
    return PathCoordinateModel.from_graph(parse_gfa(pg.write(tmp_path / "g.gfa")))


def _types(blocks) -> list[BlockType]:
    return [b.block_type for b in blocks]


class TestCollinear:
    def test_single_collinear_block(self, tmp_path):
        m = _model(collinear_pangenome(n_genomes=2, n_segments=6, seg_len=10), tmp_path)
        blocks = build_pairwise_blocks(m, "sample1#0#chr1", "sample0#0#chr1")
        assert len(blocks) == 1
        b = blocks[0]
        assert b.block_type is BlockType.COLLINEAR
        assert b.strand == "+"
        assert b.n_anchors == 6
        assert b.query.contig == "chr1"
        assert (b.target.start, b.target.end) == (0, 60)


class TestInversion:
    def test_inversion_block_detected(self, tmp_path):
        m = _model(inversion_pangenome(seg_len=10), tmp_path)
        blocks = build_pairwise_blocks(m, "sample3#0#chr1", "sample0#0#chr1")
        types = _types(blocks)
        assert types.count(BlockType.INVERSION) == 1
        assert types.count(BlockType.COLLINEAR) == 2   # flanks
        inv = next(b for b in blocks if b.block_type is BlockType.INVERSION)
        assert inv.strand == "-"
        # the inverted run covers s3,s4 -> reference span chr1:20-40
        assert (inv.target.start, inv.target.end) == (20, 40)
        assert {a.name for a in inv.anchors} == {"s3", "s4"}

    def test_self_comparison_is_collinear(self, tmp_path):
        m = _model(inversion_pangenome(), tmp_path)
        blocks = build_pairwise_blocks(m, "sample0#0#chr1", "sample0#0#chr1")
        assert _types(blocks) == [BlockType.COLLINEAR]


class TestTranslocation:
    def test_moved_block_typed_translocation(self, tmp_path):
        m = _model(translocation_pangenome(seg_len=10), tmp_path)
        # sample2 reorders to s1,s4,s5,s2,s3
        blocks = build_pairwise_blocks(m, "sample2#0#chr1", "sample0#0#chr1")
        types = _types(blocks)
        assert BlockType.TRANSLOCATION in types
        # the translocated block is the moved s2,s3
        trans = next(b for b in blocks if b.block_type is BlockType.TRANSLOCATION)
        assert {a.name for a in trans.anchors} == {"s2", "s3"}

    def test_collinear_genome_has_no_translocation(self, tmp_path):
        m = _model(translocation_pangenome(), tmp_path)
        blocks = build_pairwise_blocks(m, "sample1#0#chr1", "sample0#0#chr1")
        assert _types(blocks) == [BlockType.COLLINEAR]


class TestDuplication:
    def test_cnv_segment_becomes_duplication_block(self, tmp_path):
        m = _model(duplication_pangenome(seg_len=10), tmp_path)
        # sample2 has s2 duplicated
        blocks = build_pairwise_blocks(m, "sample2#0#chr1", "sample0#0#chr1")
        types = _types(blocks)
        assert types.count(BlockType.DUPLICATION) == 1
        dup = next(b for b in blocks if b.block_type is BlockType.DUPLICATION)
        assert {a.name for a in dup.anchors} == {"s2"}
        # the backbone (deduped) is a single collinear block
        assert types.count(BlockType.COLLINEAR) == 1


class TestInsertionsSkipped:
    def test_query_private_segment_skipped(self, tmp_path):
        from privy.synthetic import SyntheticPangenome

        pg = SyntheticPangenome()
        for s in ("s1", "s2", "s3", "x1"):
            pg.add_segment(s, 10)
        pg.add_genome("sample0#0#chr1", [("s1", "+"), ("s2", "+"), ("s3", "+")])
        # query has an extra insertion x1 not on the reference
        pg.add_genome("sample1#0#chr1", [("s1", "+"), ("x1", "+"), ("s2", "+"), ("s3", "+")])
        m = _model(pg, tmp_path)
        blocks = build_pairwise_blocks(m, "sample1#0#chr1", "sample0#0#chr1")
        # x1 is skipped; s1,s2,s3 remain collinear
        assert _types(blocks) == [BlockType.COLLINEAR]
        names = {a.name for b in blocks for a in b.anchors}
        assert "x1" not in names

    def test_min_block_anchors_filter(self, tmp_path):
        m = _model(collinear_pangenome(n_genomes=2, n_segments=6), tmp_path)
        blocks = build_pairwise_blocks(
            m, "sample1#0#chr1", "sample0#0#chr1", min_block_anchors=10
        )
        assert blocks == []   # the only block has 6 anchors < 10
