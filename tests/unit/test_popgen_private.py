"""Unit tests for popgen private-allele metrics and F_IS."""

from __future__ import annotations

from pathlib import Path

import pytest

from privy.io.gfa import parse_gfa
from privy.microhap import detect_microhaplotypes
from privy.polyploid import group_paths_by_sample
from privy.popgen.diversity import inbreeding_fis
from privy.popgen.private import (
    private_allele_counts,
    rarefied_private_allelic_richness,
)
from privy.synteny.coordinates import PathCoordinateModel
from privy.synthetic import autopolyploid_pangenome, microhaplotype_pangenome


def _loci(pg, tmp_path: Path):
    graph = parse_gfa(pg.write(tmp_path / "g.gfa"))
    model = PathCoordinateModel.from_graph(graph)
    return graph, model, detect_microhaplotypes(graph, model, model.path_ids()[0])


TARGETS = ["sample0#0#chr1", "sample1#0#chr1"]
OFFTARGETS = ["sample2#0#chr1", "sample3#0#chr1"]


class TestPrivateAlleleCounts:
    def test_one_private_allele_each_cohort(self, tmp_path):
        _g, _m, loci = _loci(microhaplotype_pangenome(seg_len=10), tmp_path)
        counts = private_allele_counts(loci, {"target": TARGETS, "offtarget": OFFTARGETS})
        assert counts == {"target": 1, "offtarget": 1}

    def test_no_private_when_shared(self, tmp_path):
        _g, _m, loci = _loci(microhaplotype_pangenome(seg_len=10), tmp_path)
        mixed_a = ["sample0#0#chr1", "sample2#0#chr1"]
        mixed_b = ["sample1#0#chr1", "sample3#0#chr1"]
        counts = private_allele_counts(loci, {"a": mixed_a, "b": mixed_b})
        assert counts == {"a": 0, "b": 0}


class TestRarefiedPrivateRichness:
    def test_rarefied_richness_matches_hand_calc(self, tmp_path):
        _g, _m, loci = _loci(microhaplotype_pangenome(seg_len=10), tmp_path)
        rich = rarefied_private_allelic_richness(
            loci, {"target": TARGETS, "offtarget": OFFTARGETS}, g=2
        )
        # each cohort fixes its own allele (present prob 1, absent-elsewhere prob 1)
        assert rich["target"] == pytest.approx(1.0)
        assert rich["offtarget"] == pytest.approx(1.0)


class TestFis:
    def test_fis_positive_with_het_deficit(self, tmp_path):
        # autopolyploid AA/AB/BB: He=0.5, Ho=1/3 -> FIS = 1 - (1/3)/0.5 = 1/3
        _g, _m, loci = _loci(autopolyploid_pangenome(ploidy=2, seg_len=10), tmp_path)
        graph = parse_gfa(autopolyploid_pangenome(ploidy=2, seg_len=10).write(tmp_path / "h.gfa"))
        model = PathCoordinateModel.from_graph(graph)
        groups = group_paths_by_sample(model.path_ids())
        fis = inbreeding_fis(loci[0], groups)
        assert fis == pytest.approx(1 / 3)

    def test_fis_none_when_monomorphic(self, tmp_path):
        from privy.microhap.model import Microhaplotype

        mono = Microhaplotype(
            locus_id="L", contig="c", start=0, end=10,
            alleles={"s#0#c": "a", "s#1#c": "a"}, ref_allele="a",
        )
        assert inbreeding_fis(mono, {"s": ["s#0#c", "s#1#c"]}) is None
