"""Unit tests for src/privy/io/hvcf.py (PHG hVCF read/write)."""

from __future__ import annotations

from pathlib import Path

from privy.io.gfa import parse_gfa
from privy.io.hvcf import read_hvcf, write_hvcf
from privy.microhap import detect_microhaplotypes
from privy.microhap.model import Microhaplotype
from privy.synteny.coordinates import PathCoordinateModel
from privy.synthetic import microhaplotype_pangenome


def _loci(tmp_path: Path) -> list[Microhaplotype]:
    graph = parse_gfa(microhaplotype_pangenome(seg_len=10).write(tmp_path / "g.gfa"))
    model = PathCoordinateModel.from_graph(graph)
    return detect_microhaplotypes(graph, model, "sample0#0#chr1")


class TestWriteHvcf:
    def test_header_and_record(self, tmp_path):
        loci = _loci(tmp_path)
        out = write_hvcf(loci, tmp_path / "h.vcf")
        text = out.read_text(encoding="utf-8")
        assert "##fileformat=VCFv4.2" in text
        assert "##ALT=<ID=" in text                 # symbolic ALT (MD5) alleles
        assert "END=" in text
        assert text.count("\n") >= 5

    def test_md5_allele_ids_in_alt(self, tmp_path):
        loci = _loci(tmp_path)
        out = write_hvcf(loci, tmp_path / "h.vcf")
        # the locus's distinct allele ids appear as symbolic ALTs
        distinct = set(loci[0].alleles.values())
        text = out.read_text(encoding="utf-8")
        for allele in distinct:
            assert f"<{allele}>" in text


class TestRoundTrip:
    def test_roundtrip_preserves_alleles(self, tmp_path):
        loci = _loci(tmp_path)
        out = write_hvcf(loci, tmp_path / "h.vcf")
        recs = read_hvcf(out)
        assert len(recs) == 1
        rec = recs[0]
        mh = loci[0]
        assert rec.locus_id == mh.locus_id
        assert rec.contig == mh.contig
        assert rec.start == mh.start
        assert rec.end == mh.end
        # per-genome allele assignment recovered exactly
        assert rec.alleles == dict(mh.alleles)

    def test_explicit_sample_order(self, tmp_path):
        loci = _loci(tmp_path)
        samples = ["sample0#0#chr1", "sample1#0#chr1", "sample2#0#chr1", "sample3#0#chr1"]
        out = write_hvcf(loci, tmp_path / "h.vcf", samples=samples)
        header = next(
            ln for ln in out.read_text().splitlines() if ln.startswith("#CHROM")
        )
        assert header.split("\t")[9:] == samples
