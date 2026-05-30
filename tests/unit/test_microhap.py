"""Unit tests for the microhaplotype layer (src/privy/microhap)."""

from __future__ import annotations

from pathlib import Path

from privy.io.gfa import parse_gfa
from privy.microhap import Microhaplotype, detect_microhaplotypes
from privy.synteny.coordinates import PathCoordinateModel
from privy.synthetic import (
    collinear_pangenome,
    microhaplotype_pangenome,
    presence_absence_pangenome,
)


def _graph_and_model(pg, tmp_path: Path):
    graph = parse_gfa(pg.write(tmp_path / "g.gfa"))
    return graph, PathCoordinateModel.from_graph(graph)


# ---------------------------------------------------------------------------
# Detection
# ---------------------------------------------------------------------------


class TestDetect:
    def test_multiallelic_bubble_detected(self, tmp_path):
        graph, model = _graph_and_model(microhaplotype_pangenome(seg_len=10), tmp_path)
        loci = detect_microhaplotypes(graph, model, "sample0#0#chr1")
        assert len(loci) == 1
        mh = loci[0]
        assert mh.is_multiallelic
        assert mh.n_alleles == 2          # sT vs sR
        assert mh.n_genomes == 4
        assert mh.contig == "chr1"

    def test_collinear_has_no_multiallelic_loci(self, tmp_path):
        graph, model = _graph_and_model(
            collinear_pangenome(n_genomes=4, n_segments=6), tmp_path
        )
        loci = detect_microhaplotypes(graph, model, "sample0#0#chr1")
        assert loci == []   # everyone identical -> no variation

    def test_presence_absence_is_multiallelic(self, tmp_path):
        # targets carry [s2,s3]; off-targets deleted them -> 2 alleles ("" vs s2s3)
        graph, model = _graph_and_model(presence_absence_pangenome(seg_len=10), tmp_path)
        loci = detect_microhaplotypes(graph, model, "sample0#0#chr1")
        assert len(loci) == 1
        assert loci[0].n_alleles == 2

    def test_unknown_reference_raises(self, tmp_path):
        graph, model = _graph_and_model(microhaplotype_pangenome(), tmp_path)
        try:
            detect_microhaplotypes(graph, model, "ghost#0#chr1")
        except KeyError as exc:
            assert "ghost" in str(exc)
        else:  # pragma: no cover
            raise AssertionError("expected KeyError")


# ---------------------------------------------------------------------------
# Allele identity / frequencies
# ---------------------------------------------------------------------------


class TestAlleles:
    def test_shared_allele_has_same_md5(self, tmp_path):
        graph, model = _graph_and_model(microhaplotype_pangenome(seg_len=10), tmp_path)
        mh = detect_microhaplotypes(graph, model, "sample0#0#chr1")[0]
        # the two targets share an allele; the two off-targets share the other
        assert mh.alleles["sample0#0#chr1"] == mh.alleles["sample1#0#chr1"]
        assert mh.alleles["sample2#0#chr1"] == mh.alleles["sample3#0#chr1"]
        assert mh.alleles["sample0#0#chr1"] != mh.alleles["sample2#0#chr1"]

    def test_allele_frequencies_and_aaf(self, tmp_path):
        graph, model = _graph_and_model(microhaplotype_pangenome(seg_len=10), tmp_path)
        mh = detect_microhaplotypes(graph, model, "sample0#0#chr1")[0]
        freqs = mh.allele_frequencies()
        assert set(freqs.values()) == {0.5}            # 2 + 2 of 4
        # ref allele (sample0's sT) has freq 0.5 -> AAF = 0.5
        assert mh.aaf() == 0.5

    def test_allele_counts(self, tmp_path):
        graph, model = _graph_and_model(microhaplotype_pangenome(), tmp_path)
        mh = detect_microhaplotypes(graph, model, "sample0#0#chr1")[0]
        assert sorted(mh.allele_counts().values()) == [2, 2]


# ---------------------------------------------------------------------------
# Private microhaplotypes (the core tie-in)
# ---------------------------------------------------------------------------


class TestPrivateMicrohaplotypes:
    def test_target_private_allele_detected(self, tmp_path):
        graph, model = _graph_and_model(microhaplotype_pangenome(seg_len=10), tmp_path)
        mh = detect_microhaplotypes(graph, model, "sample0#0#chr1")[0]
        targets = ["sample0#0#chr1", "sample1#0#chr1"]
        off_targets = ["sample2#0#chr1", "sample3#0#chr1"]
        private = mh.private_alleles(targets, off_targets)
        assert len(private) == 1
        # the private allele is the targets' shared sT allele
        assert private[0] == mh.alleles["sample0#0#chr1"]
        assert mh.is_target_private(targets, off_targets)

    def test_not_private_when_allele_shared_across_cohorts(self, tmp_path):
        graph, model = _graph_and_model(microhaplotype_pangenome(), tmp_path)
        mh = detect_microhaplotypes(graph, model, "sample0#0#chr1")[0]
        # put one target and one off-target in the same cohort split -> no private allele
        targets = ["sample0#0#chr1", "sample2#0#chr1"]   # one sT, one sR
        off_targets = ["sample1#0#chr1", "sample3#0#chr1"]
        assert mh.is_target_private(targets, off_targets) is False

    def test_microhaplotype_type(self, tmp_path):
        graph, model = _graph_and_model(microhaplotype_pangenome(), tmp_path)
        mh = detect_microhaplotypes(graph, model, "sample0#0#chr1")[0]
        assert isinstance(mh, Microhaplotype)
