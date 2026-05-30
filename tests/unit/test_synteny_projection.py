"""Unit tests for src/privy/synteny/projection.py."""

from __future__ import annotations

from pathlib import Path

import pytest

from privy.io.gfa import parse_gfa
from privy.synteny.coordinates import PathCoordinateModel
from privy.synteny.model import GenomeInterval
from privy.synteny.projection import (
    CoordinateProjection,
    ProjectionStatus,
    lift_intervals,
    project_coordinate,
    project_node_set,
    project_region,
)
from privy.synthetic import (
    collinear_pangenome,
    duplication_pangenome,
    inversion_pangenome,
)


def _model(pg, tmp_path: Path) -> PathCoordinateModel:
    return PathCoordinateModel.from_graph(parse_gfa(pg.write(tmp_path / "g.gfa")))


# ---------------------------------------------------------------------------
# Single-coordinate projection
# ---------------------------------------------------------------------------


class TestProjectCoordinate:
    def test_collinear_identity(self, tmp_path):
        m = _model(collinear_pangenome(n_genomes=3, n_segments=6, seg_len=10), tmp_path)
        proj = project_coordinate(m, "sample0#0#chr1", 25, "sample1#0#chr1")
        assert proj.status is ProjectionStatus.MAPPED
        assert proj.best == ("chr1", 25)   # collinear -> identical coordinate

    def test_projects_into_inversion_keeps_same_base(self, tmp_path):
        # sample3 has s3,s4 inverted. seg_len=10, segments s1..s6.
        # s3 spans path-local [20,30) on a forward genome; pick position 22 (offset 2 into s3).
        m = _model(inversion_pangenome(seg_len=10), tmp_path)
        # forward reference position 22 -> s3 offset 2
        ref = project_coordinate(m, "sample0#0#chr1", 22, "sample0#0#chr1")
        assert ref.segment == "s3"
        assert ref.best == ("chr1", 22)
        # project the SAME physical base onto the inverted genome (sample3)
        proj = project_coordinate(m, "sample0#0#chr1", 22, "sample3#0#chr1")
        assert proj.status is ProjectionStatus.MAPPED
        assert proj.segment == "s3"
        # On sample3, order is s1,s2,s4,s3,s5,s6 -> s3 occupies path-local [30,40),
        # reversed. Base index 2 in s3 maps to offset (10-1-2)=7 -> 30+7 = 37.
        assert proj.best == ("chr1", 37)

    def test_inversion_roundtrip_returns_to_origin(self, tmp_path):
        m = _model(inversion_pangenome(seg_len=10), tmp_path)
        fwd = project_coordinate(m, "sample0#0#chr1", 22, "sample3#0#chr1")
        (contig, coord) = fwd.best
        back = project_coordinate(m, "sample3#0#chr1", coord, "sample0#0#chr1")
        assert back.best == ("chr1", 22)   # round-trip is exact

    def test_absent_segment(self, tmp_path):
        # build two genomes that share no segments
        from privy.synthetic import SyntheticPangenome

        pg = SyntheticPangenome()
        pg.add_segment("a1", 10).add_segment("b1", 10)
        pg.add_genome("g0#0#chr1", [("a1", "+")], cohort="target")
        pg.add_genome("g1#0#chr1", [("b1", "+")], cohort="offtarget")
        m = _model(pg, tmp_path)
        proj = project_coordinate(m, "g0#0#chr1", 3, "g1#0#chr1")
        assert proj.status is ProjectionStatus.ABSENT
        assert proj.best is None

    def test_duplication_is_ambiguous(self, tmp_path):
        # sample2 has s2 duplicated; project a base of s2 from the reference onto it
        m = _model(duplication_pangenome(seg_len=10), tmp_path)
        proj = project_coordinate(m, "sample0#0#chr1", 13, "sample2#0#chr1")  # s2 offset 3
        assert proj.segment == "s2"
        assert proj.status is ProjectionStatus.AMBIGUOUS
        assert len(proj.targets) == 2
        # the two copies sit 10 bp apart (tandem)
        coords = sorted(c for _, c in proj.targets)
        assert coords[1] - coords[0] == 10

    def test_out_of_range_raises(self, tmp_path):
        m = _model(collinear_pangenome(), tmp_path)
        with pytest.raises(IndexError):
            project_coordinate(m, "sample0#0#chr1", 10_000, "sample1#0#chr1")


# ---------------------------------------------------------------------------
# Region / node-set projection
# ---------------------------------------------------------------------------


class TestProjectRegion:
    def test_region_projects_to_all(self, tmp_path):
        m = _model(collinear_pangenome(n_genomes=3, n_segments=6, seg_len=10), tmp_path)
        pm = project_region(m, "sample0#0#chr1", 15, 35)   # covers s2,s3,s4 (10-40)
        # collinear: every genome gets the spanning interval [10,40)
        for target in ("sample0#0#chr1", "sample1#0#chr1", "sample2#0#chr1"):
            assert pm.projections[target] == GenomeInterval(
                genome=target.split("#")[0], contig="chr1", start=10, end=40
            )
        assert pm.present_in() == ("sample0#0#chr1", "sample1#0#chr1", "sample2#0#chr1")

    def test_region_excludes_source_when_requested(self, tmp_path):
        m = _model(collinear_pangenome(n_genomes=3), tmp_path)
        pm = project_region(m, "sample0#0#chr1", 0, 20, include_source=False)
        assert "sample0#0#chr1" not in pm.projections

    def test_region_absent_target(self, tmp_path):
        from privy.synthetic import SyntheticPangenome

        pg = SyntheticPangenome()
        pg.add_segment("a1", 10).add_segment("a2", 10).add_segment("b1", 10)
        pg.add_genome("g0#0#chr1", [("a1", "+"), ("a2", "+")])
        pg.add_genome("g1#0#chr1", [("b1", "+")])
        m = _model(pg, tmp_path)
        pm = project_region(m, "g0#0#chr1", 0, 20, targets=["g1#0#chr1"])
        assert pm.projections["g1#0#chr1"] is None
        assert pm.absent_in() == ("g1#0#chr1",)


class TestProjectNodeSet:
    def test_node_set_to_all_paths(self, tmp_path):
        m = _model(collinear_pangenome(n_genomes=3, n_segments=6, seg_len=10), tmp_path)
        pm = project_node_set(m, ["s3", "s4"], source_label="region-of-interest")
        assert pm.source == "region-of-interest"
        # s3,s4 -> path-local [20,40) -> stable chr1:20-40 on every collinear genome
        for target in m.path_ids():
            assert pm.projections[target] == GenomeInterval(
                genome=target.split("#")[0], contig="chr1", start=20, end=40
            )

    def test_node_set_subset_targets(self, tmp_path):
        m = _model(collinear_pangenome(n_genomes=4), tmp_path)
        pm = project_node_set(m, ["s1"], targets=["sample2#0#chr1"])
        assert list(pm.projections) == ["sample2#0#chr1"]

    def test_returns_coordinateprojection_type(self, tmp_path):
        m = _model(collinear_pangenome(), tmp_path)
        proj = project_coordinate(m, "sample0#0#chr1", 0, "sample1#0#chr1")
        assert isinstance(proj, CoordinateProjection)


# ---------------------------------------------------------------------------
# Annotation-track liftover
# ---------------------------------------------------------------------------


class TestLiftIntervals:
    def test_collinear_lift_identity(self, tmp_path):
        m = _model(collinear_pangenome(n_genomes=3, n_segments=6, seg_len=10), tmp_path)
        feats = [
            GenomeInterval("sample0", "chr1", 5, 25),    # a "gene" on s1-s3
            GenomeInterval("sample0", "chr1", 40, 50),   # on s5
        ]
        lifted = lift_intervals(m, feats, "sample0#0#chr1", "sample1#0#chr1")
        assert lifted[0] == GenomeInterval("sample1", "chr1", 0, 30)   # spans s1,s2,s3
        assert lifted[1] == GenomeInterval("sample1", "chr1", 40, 50)  # s5

    def test_lift_through_inversion(self, tmp_path):
        m = _model(inversion_pangenome(seg_len=10), tmp_path)
        # a feature on s3 (chr1:20-30 on the reference) lifts onto the inverted genome
        feats = [GenomeInterval("sample0", "chr1", 20, 30)]
        lifted = lift_intervals(m, feats, "sample0#0#chr1", "sample3#0#chr1")
        # s3 sits at path-local [30,40) on sample3 -> stable chr1:30-40
        assert lifted[0] == GenomeInterval("sample3", "chr1", 30, 40)

    def test_lift_absent_returns_none(self, tmp_path):
        from privy.synthetic import SyntheticPangenome

        pg = SyntheticPangenome()
        pg.add_segment("a1", 10).add_segment("a2", 10).add_segment("b1", 10)
        pg.add_genome("g0#0#chr1", [("a1", "+"), ("a2", "+")])
        pg.add_genome("g1#0#chr1", [("b1", "+")])
        m = _model(pg, tmp_path)
        lifted = lift_intervals(m, [GenomeInterval("g0", "chr1", 0, 20)], "g0#0#chr1", "g1#0#chr1")
        assert lifted == [None]

    def test_lift_out_of_range_returns_none(self, tmp_path):
        m = _model(collinear_pangenome(), tmp_path)
        feats = [GenomeInterval("sample0", "chr1", 5_000, 5_010)]
        lifted = lift_intervals(m, feats, "sample0#0#chr1", "sample1#0#chr1")
        assert lifted == [None]


class TestToPathLocal:
    def test_inverse_of_to_stable_pline(self, tmp_path):
        m = _model(collinear_pangenome(), tmp_path)
        assert m.to_path_local("sample0#0#chr1", 25) == 25

    def test_out_of_range(self, tmp_path):
        m = _model(collinear_pangenome(), tmp_path)
        with pytest.raises(IndexError):
            m.to_path_local("sample0#0#chr1", 10_000)
