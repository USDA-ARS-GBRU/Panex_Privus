"""Unit tests for src/privy/popgen/relationship.py (dosage matrix + VanRaden GRM)."""

from __future__ import annotations

from pathlib import Path

import pytest

from privy.io.gfa import parse_gfa
from privy.microhap import detect_microhaplotypes
from privy.polyploid import group_paths_by_sample
from privy.popgen.relationship import build_dosage_matrix, vanraden_grm
from privy.synteny.coordinates import PathCoordinateModel
from privy.synthetic import autopolyploid_pangenome


def _loci_and_groups(tmp_path: Path, ploidy: int = 2):
    graph = parse_gfa(autopolyploid_pangenome(ploidy=ploidy, seg_len=10).write(tmp_path / "g.gfa"))
    model = PathCoordinateModel.from_graph(graph)
    loci = detect_microhaplotypes(graph, model, "sample0#0#chr1")
    return loci, group_paths_by_sample(model.path_ids())


class TestDosageMatrix:
    def test_shape_and_dosage_gradient(self, tmp_path):
        loci, groups = _loci_and_groups(tmp_path, ploidy=2)
        dm = build_dosage_matrix(loci, groups)
        assert dm.ploidy == 2
        assert dm.samples == ("sample0", "sample1", "sample2")
        assert len(dm.locus_ids) == 1
        # one column; dosages 0,1,2 across samples
        assert [row[0] for row in dm.matrix] == [0, 1, 2]


class TestVanRadenGrm:
    def test_grm_is_symmetric_and_expected(self, tmp_path):
        loci, groups = _loci_and_groups(tmp_path, ploidy=2)
        dm = build_dosage_matrix(loci, groups)
        samples, grm = vanraden_grm(dm)
        assert samples == ("sample0", "sample1", "sample2")
        # symmetry
        for i in range(3):
            for j in range(3):
                assert grm[i][j] == pytest.approx(grm[j][i])
        # one locus, dosages [0,1,2], p=0.5, P=1, Z=[-1,0,1], denom=2*0.5*0.5=0.5
        # G = ZZ'/denom -> diagonal [2, 0, 2]
        assert grm[0][0] == pytest.approx(2.0)
        assert grm[1][1] == pytest.approx(0.0)
        assert grm[2][2] == pytest.approx(2.0)
        assert grm[0][2] == pytest.approx(-2.0)

    def test_empty_raises(self, tmp_path):
        from privy.popgen.relationship import DosageMatrix

        with pytest.raises(ValueError):
            vanraden_grm(DosageMatrix(samples=(), locus_ids=(), ploidy=2, matrix=[]))

    def test_monomorphic_returns_zeros(self, tmp_path):
        from privy.popgen.relationship import DosageMatrix

        # all samples identical dosage -> no variance -> zeros
        dm = DosageMatrix(
            samples=("a", "b"), locus_ids=("L1",), ploidy=2, matrix=[[2], [2]]
        )
        _samples, grm = vanraden_grm(dm)
        assert grm == [[0.0, 0.0], [0.0, 0.0]]
