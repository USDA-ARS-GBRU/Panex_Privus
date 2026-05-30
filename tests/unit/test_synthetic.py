"""Unit tests for src/privy/synthetic.py.

These also serve as integration smoke tests: every synthetic graph is parsed by
the real GFA parser and fed to PathCoordinateModel, confirming the P0 stack works
on realistic-shaped inputs and recovers the planted structure.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from privy.io.gfa import parse_gfa
from privy.synteny.coordinates import PathCoordinateModel
from privy.synthetic import (
    SyntheticPangenome,
    allopolyploid_pangenome,
    collinear_pangenome,
    duplication_pangenome,
    inversion_pangenome,
)


def _parse(pg: SyntheticPangenome, tmp_path: Path, *, use_walks: bool = False):
    return parse_gfa(pg.write(tmp_path / "g.gfa", use_walks=use_walks))


# ---------------------------------------------------------------------------
# Determinism & builder basics
# ---------------------------------------------------------------------------


class TestBuilder:
    def test_deterministic_output(self):
        assert collinear_pangenome().to_gfa() == collinear_pangenome().to_gfa()

    def test_rejects_unknown_segment(self):
        pg = SyntheticPangenome().add_segment("s1", 10)
        with pytest.raises(ValueError, match="unknown segment"):
            pg.add_genome("g#0#c", [("s2", "+")])

    def test_rejects_duplicate_segment(self):
        pg = SyntheticPangenome().add_segment("s1", 10)
        with pytest.raises(ValueError, match="duplicate"):
            pg.add_segment("s1", 5)

    def test_cohort_accessors(self):
        pg = collinear_pangenome(n_genomes=4, n_target=2)
        assert pg.cohort("target") == ["sample0#0#chr1", "sample1#0#chr1"]
        assert pg.cohort("offtarget") == ["sample2#0#chr1", "sample3#0#chr1"]


# ---------------------------------------------------------------------------
# Collinear: parses, PanSN, rGFA tags, coordinate model
# ---------------------------------------------------------------------------


class TestCollinear:
    def test_parses_and_builds_model(self, tmp_path):
        pg = collinear_pangenome(n_genomes=4, n_segments=6, seg_len=10)
        graph = _parse(pg, tmp_path)
        assert len(graph.paths) == 4
        model = PathCoordinateModel.from_graph(graph)
        assert model.path_length("sample0#0#chr1") == 60   # 6 * 10
        # collinear: every genome maps the same path position to the same contig pos
        assert model.to_stable("sample0#0#chr1", 25) == ("chr1", 25)
        assert model.to_stable("sample3#0#chr1", 25) == ("chr1", 25)

    def test_reference_has_rgfa_tags(self, tmp_path):
        pg = collinear_pangenome()
        graph = _parse(pg, tmp_path)
        s1 = graph.segments["s1"]
        assert s1.ref_contig == "chr1"
        assert s1.ref_start == 0
        assert graph.segments["s2"].ref_start == 10

    def test_walk_mode_parses(self, tmp_path):
        pg = collinear_pangenome()
        graph = _parse(pg, tmp_path, use_walks=True)
        assert len(graph.walks) == 4
        model = PathCoordinateModel.from_graph(graph)
        assert "sample0#0#chr1" in model


# ---------------------------------------------------------------------------
# Planted structural events recovered via the coordinate model
# ---------------------------------------------------------------------------


class TestInversion:
    def test_inverted_run_has_reverse_orientation(self, tmp_path):
        pg = inversion_pangenome()
        model = PathCoordinateModel.from_graph(_parse(pg, tmp_path))
        # target genome (sample3) carries inverted s3/s4
        (occ_s3,) = model.occurrences("sample3#0#chr1", "s3")
        (occ_s4,) = model.occurrences("sample3#0#chr1", "s4")
        assert occ_s3.orientation == "-"
        assert occ_s4.orientation == "-"
        # s4 now precedes s3 on the inverted path
        assert occ_s4.start < occ_s3.start
        # an off-target genome keeps forward orientation
        (ref_s3,) = model.occurrences("sample0#0#chr1", "s3")
        assert ref_s3.orientation == "+"


class TestDuplication:
    def test_tandem_duplication_two_occurrences(self, tmp_path):
        pg = duplication_pangenome()
        model = PathCoordinateModel.from_graph(_parse(pg, tmp_path))
        dup = model.occurrences("sample2#0#chr1", "s2")
        assert len(dup) == 2   # CNV: s2 appears twice on the target genome
        # the off-target reference genome has exactly one copy
        assert len(model.occurrences("sample0#0#chr1", "s2")) == 1


# ---------------------------------------------------------------------------
# Allopolyploid: two subgenomes per sample
# ---------------------------------------------------------------------------


class TestAllopolyploid:
    def test_two_subgenome_paths_per_sample(self, tmp_path):
        pg = allopolyploid_pangenome()
        graph = _parse(pg, tmp_path)
        model = PathCoordinateModel.from_graph(graph)
        ids = set(model.path_ids())
        assert "sample0#0#chrA" in ids
        assert "sample0#1#chrD" in ids
        # subgenomes use disjoint segments
        assert model.occurrences("sample0#0#chrA", "d1") == []
        assert model.stable_contig("sample0#1#chrD") == "chrD"

    def test_target_cohort_labelled(self):
        pg = allopolyploid_pangenome()
        assert "sample0#0#chrA" in pg.cohort("target")
        assert "sample1#0#chrA" in pg.cohort("offtarget")
