"""Unit tests for src/privy/popgen/structure.py (PCA / DAPC)."""

from __future__ import annotations

import pytest

from privy.popgen.relationship import DosageMatrix
from privy.popgen.structure import dapc, labels_from_cohorts, pca
from privy.utils.optional import MissingDependencyError, is_available


def _dm() -> DosageMatrix:
    # 4 samples, 3 loci; two clear groups: {a,b} low dosage, {c,d} high dosage
    return DosageMatrix(
        samples=("a", "b", "c", "d"),
        locus_ids=("L1", "L2", "L3"),
        ploidy=2,
        matrix=[
            [0, 0, 1],
            [0, 1, 0],
            [2, 2, 1],
            [2, 1, 2],
        ],
    )


class TestPca:
    def test_shape_and_variance(self):
        result = pca(_dm(), n_components=2)
        assert result.samples == ("a", "b", "c", "d")
        assert len(result.coords) == 4
        assert all(len(row) == result.n_components for row in result.coords)
        assert result.n_components >= 1
        # explained variance ratios are in [0,1] and sorted descending
        evr = result.explained_variance_ratio
        assert all(0.0 <= v <= 1.0 for v in evr)
        assert evr == sorted(evr, reverse=True)

    def test_groups_separate_on_pc1(self):
        result = pca(_dm(), n_components=1)
        pc1 = [row[0] for row in result.coords]
        # a,b (low) should fall on the opposite side of c,d (high) along PC1
        assert (pc1[0] < 0) == (pc1[1] < 0)          # a,b same side
        assert (pc1[2] < 0) == (pc1[3] < 0)          # c,d same side
        assert (pc1[0] < 0) != (pc1[2] < 0)          # groups opposite sides

    def test_monomorphic_returns_empty(self):
        dm = DosageMatrix(samples=("a", "b"), locus_ids=("L1",), ploidy=2, matrix=[[2], [2]])
        result = pca(dm)
        assert result.n_components == 0


class TestLabels:
    def test_labels_from_cohorts(self):
        labels = labels_from_cohorts(["t0#0#c", "t1#0#c"], ["o0#0#c"])
        assert labels == {"t0": "target", "t1": "target", "o0": "offtarget"}


class TestDapc:
    def test_dapc_runs_or_degrades(self):
        dm = _dm()
        labels = {"a": "lo", "b": "lo", "c": "hi", "d": "hi"}
        if is_available("sklearn"):
            result = dapc(dm, labels)
            assert set(result.assigned) <= {"lo", "hi"}
            assert len(result.coords) == 4
        else:
            with pytest.raises(MissingDependencyError):
                dapc(dm, labels)

    def test_dapc_needs_two_groups(self):
        if not is_available("sklearn"):
            pytest.skip("sklearn not installed")
        dm = _dm()
        with pytest.raises(ValueError):
            dapc(dm, {"a": "x", "b": "x", "c": "x", "d": "x"})
