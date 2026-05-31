"""Unit tests for the population-genetics layer (src/privy/popgen)."""

from __future__ import annotations

from pathlib import Path

import pytest

from privy.io.gfa import parse_gfa
from privy.microhap import detect_microhaplotypes
from privy.popgen import (
    allele_frequencies,
    effective_n_alleles,
    genome_wide_fst,
    locus_differentiation,
    locus_diversity,
    nei_gene_diversity,
)
from privy.synteny.coordinates import PathCoordinateModel
from privy.synthetic import microhaplotype_pangenome


def _loci(tmp_path: Path):
    graph = parse_gfa(microhaplotype_pangenome(seg_len=10).write(tmp_path / "g.gfa"))
    model = PathCoordinateModel.from_graph(graph)
    return detect_microhaplotypes(graph, model, "sample0#0#chr1")


TARGETS = ["sample0#0#chr1", "sample1#0#chr1"]
OFFTARGETS = ["sample2#0#chr1", "sample3#0#chr1"]


# ---------------------------------------------------------------------------
# Diversity primitives (closed-form checks)
# ---------------------------------------------------------------------------


class TestDiversityPrimitives:
    def test_nei_gene_diversity(self):
        assert nei_gene_diversity({"a": 0.5, "b": 0.5}) == pytest.approx(0.5)
        assert nei_gene_diversity({"a": 1.0}) == pytest.approx(0.0)
        quarters = {"a": 0.25, "b": 0.25, "c": 0.25, "d": 0.25}
        assert nei_gene_diversity(quarters) == pytest.approx(0.75)

    def test_effective_n_alleles(self):
        assert effective_n_alleles({"a": 0.5, "b": 0.5}) == pytest.approx(2.0)
        assert effective_n_alleles({"a": 1.0}) == pytest.approx(1.0)

    def test_allele_frequencies_cohort_subset(self, tmp_path):
        mh = _loci(tmp_path)[0]
        # within targets only -> a single allele at frequency 1.0
        ft = allele_frequencies(mh.alleles, TARGETS)
        assert set(ft.values()) == {1.0}
        assert len(ft) == 1


# ---------------------------------------------------------------------------
# Per-locus diversity
# ---------------------------------------------------------------------------


class TestLocusDiversity:
    def test_overall_diversity(self, tmp_path):
        mh = _loci(tmp_path)[0]
        div = locus_diversity(mh)
        assert div.n_alleles == 2
        assert div.n_genomes == 4
        assert div.gene_diversity == pytest.approx(0.5)   # 1 - 2*0.5^2
        assert div.effective_alleles == pytest.approx(2.0)

    def test_within_target_is_monomorphic(self, tmp_path):
        mh = _loci(tmp_path)[0]
        div = locus_diversity(mh, genomes=TARGETS)
        assert div.n_alleles == 1
        assert div.gene_diversity == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# Differentiation (the breeder signal)
# ---------------------------------------------------------------------------


class TestDifferentiation:
    def test_private_locus_is_fully_differentiated(self, tmp_path):
        mh = _loci(tmp_path)[0]
        d = locus_differentiation(mh, TARGETS, OFFTARGETS)
        assert d is not None
        assert d.h_s == pytest.approx(0.0)    # each cohort monomorphic
        assert d.h_t == pytest.approx(0.5)
        assert d.gst == pytest.approx(1.0)    # complete differentiation
        assert d.jost_d == pytest.approx(1.0)
        assert d.is_diagnostic is True        # no shared allele
        assert d.target_n == 2 and d.offtarget_n == 2

    def test_no_differentiation_when_alleles_shared_equally(self, tmp_path):
        mh = _loci(tmp_path)[0]
        # split cohorts so each has one T and one R -> identical freqs -> GST 0
        mixed_a = ["sample0#0#chr1", "sample2#0#chr1"]
        mixed_b = ["sample1#0#chr1", "sample3#0#chr1"]
        d = locus_differentiation(mh, mixed_a, mixed_b)
        assert d is not None
        assert d.gst == pytest.approx(0.0)
        assert d.is_diagnostic is False

    def test_returns_none_when_cohort_empty(self, tmp_path):
        mh = _loci(tmp_path)[0]
        assert locus_differentiation(mh, TARGETS, ["ghost#0#chr1"]) is None

    def test_genome_wide_fst(self, tmp_path):
        loci = _loci(tmp_path)
        fst = genome_wide_fst(loci, TARGETS, OFFTARGETS)
        assert fst == pytest.approx(1.0)   # the one locus is fully diagnostic
