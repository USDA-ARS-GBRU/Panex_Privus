"""Unit tests for the polyploid layer (src/privy/polyploid)."""

from __future__ import annotations

from pathlib import Path

import pytest

from privy.io.gfa import parse_gfa
from privy.microhap import detect_microhaplotypes
from privy.polyploid import (
    alt_dosage,
    group_paths_by_sample,
    is_heterozygous,
    observed_heterozygosity,
    observed_ploidy,
    sample_allele_dosage,
)
from privy.synteny.coordinates import PathCoordinateModel
from privy.synthetic import autopolyploid_pangenome


def _locus_and_paths(tmp_path: Path, ploidy: int = 2):
    graph = parse_gfa(autopolyploid_pangenome(ploidy=ploidy, seg_len=10).write(tmp_path / "g.gfa"))
    model = PathCoordinateModel.from_graph(graph)
    loci = detect_microhaplotypes(graph, model, "sample0#0#chr1")
    return loci[0], model.path_ids()


# ---------------------------------------------------------------------------
# Grouping
# ---------------------------------------------------------------------------


class TestGrouping:
    def test_group_paths_by_sample(self, tmp_path):
        _mh, paths = _locus_and_paths(tmp_path, ploidy=2)
        groups = group_paths_by_sample(paths)
        # samples 0,1,2 each with 2 haplotype paths
        assert set(groups) == {"sample0", "sample1", "sample2"}
        assert all(observed_ploidy(p) == 2 for p in groups.values())
        assert groups["sample1"] == ["sample1#0#chr1", "sample1#1#chr1"]


# ---------------------------------------------------------------------------
# Dosage gradient (AA / AB / BB)
# ---------------------------------------------------------------------------


class TestDosage:
    def test_alt_dosage_gradient(self, tmp_path):
        mh, paths = _locus_and_paths(tmp_path, ploidy=2)
        groups = group_paths_by_sample(paths)
        # reference is sample0#0 -> allele A; sample_k has k copies of B
        assert alt_dosage(mh, groups["sample0"]) == 0   # AA
        assert alt_dosage(mh, groups["sample1"]) == 1   # AB
        assert alt_dosage(mh, groups["sample2"]) == 2   # BB

    def test_sample_allele_dosage_counts(self, tmp_path):
        mh, paths = _locus_and_paths(tmp_path, ploidy=2)
        groups = group_paths_by_sample(paths)
        d1 = sample_allele_dosage(mh, groups["sample1"])
        assert sorted(d1.values()) == [1, 1]   # one A, one B
        d0 = sample_allele_dosage(mh, groups["sample0"])
        assert list(d0.values()) == [2]        # both A

    def test_tetraploid_dosage_range(self, tmp_path):
        mh, paths = _locus_and_paths(tmp_path, ploidy=4)
        groups = group_paths_by_sample(paths)
        dosages = sorted(alt_dosage(mh, groups[s]) for s in groups)
        assert dosages == [0, 1, 2, 3, 4]      # full 0..ploidy gradient

    def test_alt_dosage_none_when_no_call(self, tmp_path):
        mh, _paths = _locus_and_paths(tmp_path, ploidy=2)
        assert alt_dosage(mh, ["ghost#0#chr1"]) is None


# ---------------------------------------------------------------------------
# Heterozygosity
# ---------------------------------------------------------------------------


class TestHeterozygosity:
    def test_is_heterozygous(self, tmp_path):
        mh, paths = _locus_and_paths(tmp_path, ploidy=2)
        groups = group_paths_by_sample(paths)
        assert is_heterozygous(mh, groups["sample1"]) is True    # AB
        assert is_heterozygous(mh, groups["sample0"]) is False   # AA
        assert is_heterozygous(mh, groups["sample2"]) is False   # BB

    def test_observed_heterozygosity(self, tmp_path):
        mh, paths = _locus_and_paths(tmp_path, ploidy=2)
        groups = group_paths_by_sample(paths)
        # 1 of 3 samples (sample1) is heterozygous
        assert observed_heterozygosity(mh, groups) == pytest.approx(1 / 3)
