"""Unit tests for src/privy/synteny/model.py."""

from __future__ import annotations

import pytest

from privy.io.paf import PafRecord
from privy.synteny.model import (
    Anchor,
    AnchorSource,
    BlockType,
    GenomeInterval,
    ProjectionMap,
    ReferenceRange,
    SyntenyBlock,
    SyntenyRegion,
    split_pansn,
)

# ---------------------------------------------------------------------------
# GenomeInterval
# ---------------------------------------------------------------------------


class TestGenomeInterval:
    def test_length_and_fields(self):
        iv = GenomeInterval("sampleA", "chr1", 100, 250)
        assert iv.length == 150
        assert iv.genome == "sampleA"
        assert iv.contig == "chr1"

    def test_validates_bounds(self):
        with pytest.raises(ValueError):
            GenomeInterval("g", "c", -1, 10)
        with pytest.raises(ValueError):
            GenomeInterval("g", "c", 10, 5)

    def test_overlaps(self):
        a = GenomeInterval("g", "chr1", 0, 100)
        assert a.overlaps(GenomeInterval("g", "chr1", 50, 150))
        assert not a.overlaps(GenomeInterval("g", "chr1", 100, 200))   # half-open touch
        assert not a.overlaps(GenomeInterval("g", "chr2", 50, 150))    # other contig

    def test_frozen_hashable(self):
        iv = GenomeInterval("g", "c", 0, 10)
        assert iv in {iv}


# ---------------------------------------------------------------------------
# Anchor + from_paf
# ---------------------------------------------------------------------------


class TestAnchor:
    def test_defaults_and_validation(self):
        a = Anchor(GenomeInterval("g1", "c", 0, 10), GenomeInterval("g2", "c", 5, 15))
        assert a.strand == "+"
        assert a.source is AnchorSource.GRAPH
        assert a.is_reverse is False
        with pytest.raises(ValueError):
            Anchor(GenomeInterval("g1", "c", 0, 10), GenomeInterval("g2", "c", 0, 10), strand="*")

    def test_from_paf_pansn_names(self):
        rec = PafRecord(
            "sampleA#0#chr1", 1000, 100, 600, "+",
            "sampleB#1#chr1", 5000, 1100, 1600, 480, 500, 60,
        )
        a = Anchor.from_paf(rec)
        assert a.query == GenomeInterval("sampleA", "chr1", 100, 600)
        assert a.target == GenomeInterval("sampleB", "chr1", 1100, 1600)
        assert a.strand == "+"
        assert a.source is AnchorSource.PAF
        assert a.identity == pytest.approx(0.96)
        assert a.score == pytest.approx(60.0)

    def test_from_paf_plain_names_reverse(self):
        rec = PafRecord("q", 800, 0, 300, "-", "t", 4000, 200, 500, 290, 300, 30)
        a = Anchor.from_paf(rec)
        assert a.query.genome == "q" and a.query.contig == "q"
        assert a.is_reverse is True


# ---------------------------------------------------------------------------
# SyntenyBlock / SyntenyRegion
# ---------------------------------------------------------------------------


class TestSyntenyBlock:
    def _block(self, bid="B1", strand="+", btype=BlockType.COLLINEAR):
        anchors = (
            Anchor(GenomeInterval("g1", "c1", 0, 10), GenomeInterval("g2", "c1", 0, 10)),
            Anchor(GenomeInterval("g1", "c1", 20, 30), GenomeInterval("g2", "c1", 20, 30)),
        )
        return SyntenyBlock(
            block_id=bid,
            query=GenomeInterval("g1", "c1", 0, 30),
            target=GenomeInterval("g2", "c1", 0, 30),
            strand=strand,
            block_type=btype,
            anchors=anchors,
        )

    def test_n_anchors_and_type(self):
        b = self._block(btype=BlockType.INVERSION, strand="-")
        assert b.n_anchors == 2
        assert b.block_type is BlockType.INVERSION
        assert b.block_type.value == "inversion"

    def test_region_aggregates(self):
        region = SyntenyRegion(
            region_id="R1",
            reference=GenomeInterval("g2", "c1", 0, 30),
            blocks=(self._block("B1"), self._block("B2")),
        )
        assert region.n_blocks == 2
        assert region.genomes == ("g1", "g2")


# ---------------------------------------------------------------------------
# ReferenceRange
# ---------------------------------------------------------------------------


class TestReferenceRange:
    def test_distinct_haplotypes_and_length(self):
        rr = ReferenceRange(
            range_id="chr1_0000000100",
            contig="chr1",
            start=100,
            end=200,
            haplotypes={"g1": "md5a", "g2": "md5a", "g3": "md5b"},
        )
        assert rr.length == 100
        assert rr.n_distinct_haplotypes == 2   # md5a, md5b

    def test_validates_bounds(self):
        with pytest.raises(ValueError):
            ReferenceRange("r", "c", 200, 100)


# ---------------------------------------------------------------------------
# ProjectionMap
# ---------------------------------------------------------------------------


class TestProjectionMap:
    def test_present_and_absent(self):
        pm = ProjectionMap(
            source="nodeset:42",
            projections={
                "g1": GenomeInterval("g1", "chr1", 100, 200),
                "g2": None,
                "g3": GenomeInterval("g3", "chr1", 500, 600),
            },
        )
        assert pm.present_in() == ("g1", "g3")
        assert pm.absent_in() == ("g2",)


# ---------------------------------------------------------------------------
# split_pansn
# ---------------------------------------------------------------------------


class TestSplitPansn:
    def test_three_field(self):
        assert split_pansn("HG002#1#chr1") == ("HG002", 1, "chr1")

    def test_three_field_nonnumeric_hap(self):
        assert split_pansn("S#mat#chr1") == ("S", None, "chr1")

    def test_two_field(self):
        assert split_pansn("sample#chr1") == ("sample", None, "chr1")

    def test_plain(self):
        assert split_pansn("chr1") == ("chr1", None, "chr1")

    def test_custom_delimiter(self):
        assert split_pansn("s|0|c", delimiter="|") == ("s", 0, "c")
