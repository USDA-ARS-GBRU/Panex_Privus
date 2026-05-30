"""Unit tests for src/privy/synteny/build.py."""

from __future__ import annotations

from pathlib import Path

from privy.io.gfa import parse_gfa
from privy.synteny.build import build_synteny, tag_region_privacy
from privy.synteny.model import BlockType
from privy.synthetic import (
    collinear_pangenome,
    inversion_pangenome,
    presence_absence_pangenome,
)


def _model(pg, tmp_path: Path):
    return parse_gfa(pg.write(tmp_path / "g.gfa"))


def _coord_model(pg, tmp_path: Path):
    from privy.synteny.coordinates import PathCoordinateModel

    return PathCoordinateModel.from_graph(_model(pg, tmp_path))


# ---------------------------------------------------------------------------
# build_synteny / region grouping
# ---------------------------------------------------------------------------


class TestBuildSynteny:
    def test_blocks_for_all_queries(self, tmp_path):
        m = _coord_model(collinear_pangenome(n_genomes=4, n_segments=5), tmp_path)
        result = build_synteny(m, "sample0#0#chr1")
        # 3 queries, each a single collinear block
        assert result.reference == "sample0#0#chr1"
        assert len(result.blocks) == 3
        assert all(b.block_type is BlockType.COLLINEAR for b in result.blocks)

    def test_explicit_query_subset(self, tmp_path):
        m = _coord_model(collinear_pangenome(n_genomes=4), tmp_path)
        result = build_synteny(m, "sample0#0#chr1", ["sample1#0#chr1"])
        assert len(result.blocks) == 1
        assert result.blocks[0].query.genome == "sample1"

    def test_collinear_blocks_group_into_one_region(self, tmp_path):
        m = _coord_model(collinear_pangenome(n_genomes=4, n_segments=5), tmp_path)
        result = build_synteny(m, "sample0#0#chr1")
        # all blocks span the same reference range -> one merged region
        assert len(result.regions) == 1
        assert result.regions[0].reference.contig == "chr1"
        assert result.regions[0].n_blocks == 3


class TestGroupRegions:
    def test_disjoint_blocks_make_separate_regions(self, tmp_path):
        # build blocks manually via inversion graph then regroup by hand is overkill;
        # use two non-overlapping synthetic queries through different segments.
        from privy.synthetic import SyntheticPangenome

        pg = SyntheticPangenome()
        for s in ("s1", "s2", "s3", "s4"):
            pg.add_segment(s, 10)
        pg.add_genome("sample0#0#chr1", [(s, "+") for s in ("s1", "s2", "s3", "s4")])
        # query covering only s1 and only s4 would still chain; instead give a query
        # that shares just s1 (region A) — separate from a query sharing just s4.
        pg.add_genome("sampleA#0#chr1", [("s1", "+")])
        pg.add_genome("sampleB#0#chr1", [("s4", "+")])
        m = _coord_model(pg, tmp_path)
        result = build_synteny(m, "sample0#0#chr1")
        # one block at chr1:0-10 (s1), one at chr1:30-40 (s4) -> two regions
        assert len(result.regions) == 2
        spans = sorted((r.reference.start, r.reference.end) for r in result.regions)
        assert spans == [(0, 10), (30, 40)]


# ---------------------------------------------------------------------------
# private-region tagging (the differentiator)
# ---------------------------------------------------------------------------


class TestTagRegionPrivacy:
    def test_detects_target_private_region(self, tmp_path):
        # targets have s2,s3; off-targets deleted them
        pg = presence_absence_pangenome(seg_len=10)
        m = _coord_model(pg, tmp_path)
        result = build_synteny(m, "sample0#0#chr1")
        privacy = tag_region_privacy(
            m,
            result.regions,
            targets=["sample1#0#chr1"],            # query target (sample0 is reference)
            off_targets=["sample2#0#chr1", "sample3#0#chr1"],
        )
        (region_id,) = privacy
        verdict = privacy[region_id]
        assert verdict.target_private is True
        assert verdict.target_present == ("sample1#0#chr1",)
        assert verdict.offtarget_present == ()   # deletion -> not present at full fraction

    def test_collinear_shared_region_is_not_private(self, tmp_path):
        m = _coord_model(collinear_pangenome(n_genomes=4, n_segments=5), tmp_path)
        result = build_synteny(m, "sample0#0#chr1")
        privacy = tag_region_privacy(
            m,
            result.regions,
            targets=["sample1#0#chr1"],
            off_targets=["sample2#0#chr1", "sample3#0#chr1"],
        )
        verdict = next(iter(privacy.values()))
        # everyone shares the region -> off-targets present -> not private
        assert verdict.target_private is False
        assert set(verdict.offtarget_present) == {"sample2#0#chr1", "sample3#0#chr1"}

    def test_inversion_region_present_in_both_cohorts(self, tmp_path):
        m = _coord_model(inversion_pangenome(seg_len=10), tmp_path)
        result = build_synteny(m, "sample0#0#chr1")
        privacy = tag_region_privacy(
            m, result.regions,
            targets=["sample3#0#chr1"],                 # the inverted genome
            off_targets=["sample1#0#chr1", "sample2#0#chr1"],
        )
        # all genomes still traverse the segments (just inverted) -> present in both
        assert all(not v.target_private for v in privacy.values())
