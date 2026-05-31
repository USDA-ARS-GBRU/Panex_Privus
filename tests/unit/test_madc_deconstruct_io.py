"""Unit tests for src/privy/io/madc.py and src/privy/io/deconstruct.py."""

from __future__ import annotations

import gzip
from pathlib import Path

from privy.io.deconstruct import read_deconstruct_sites
from privy.io.madc import call_alleles, read_madc

# ---------------------------------------------------------------------------
# MADC
# ---------------------------------------------------------------------------

_MADC = "\n".join([
    "MarkerID\tAlleleID\tAlleleType\ts1\ts2\ts3",
    "chr1_000000100\tchr1_000000100_Ref\tRef\t20\t0\t10",
    "chr1_000000100\tchr1_000000100_Alt\tAlt\t0\t18\t9",
    "chr1_000000500\tchr1_000000500_Ref\tRef\t30\t1\t0",
    "chr1_000000500\tchr1_000000500_Alt\tAlt\t1\t25\t0",
]) + "\n"


def _write(tmp_path: Path, text: str, name: str = "report.madc", gz: bool = False) -> Path:
    if gz:
        p = tmp_path / (name + ".gz")
        with gzip.open(p, "wt", encoding="utf-8") as fh:
            fh.write(text)
        return p
    p = tmp_path / name
    p.write_text(text, encoding="utf-8")
    return p


class TestReadMadc:
    def test_parses_loci_and_counts(self, tmp_path):
        loci = read_madc(_write(tmp_path, _MADC))
        assert len(loci) == 2
        first = next(loc for loc in loci if loc.locus_id == "chr1_000000100")
        assert first.samples == ("s1", "s2", "s3")
        assert len(first.alleles) == 2
        ref = next(a for a in first.alleles if a.allele_class == "Ref")
        assert ref.counts == {"s1": 20, "s2": 0, "s3": 10}

    def test_reads_gzip(self, tmp_path):
        loci = read_madc(_write(tmp_path, _MADC, gz=True))
        assert len(loci) == 2

    def test_call_alleles_dominant(self, tmp_path):
        loci = read_madc(_write(tmp_path, _MADC))
        locus = next(loc for loc in loci if loc.locus_id == "chr1_000000100")
        calls = call_alleles(locus, min_reads=2)
        assert calls["s1"] == "chr1_000000100_Ref"   # 20 vs 0
        assert calls["s2"] == "chr1_000000100_Alt"   # 0 vs 18
        # s3 is 10 vs 9 -> dominant Ref (no tie)
        assert calls["s3"] == "chr1_000000100_Ref"

    def test_call_alleles_threshold_and_tie(self, tmp_path):
        # s3 at locus 500 is 0/0 -> no call; tie handling: equal nonzero -> no call
        tie = "\n".join([
            "MarkerID\tAlleleID\tAlleleType\tsX",
            "L1\tL1_Ref\tRef\t5",
            "L1\tL1_Alt\tAlt\t5",
        ]) + "\n"
        loci = read_madc(_write(tmp_path, tie))
        assert call_alleles(loci[0]) == {}   # tie -> no call


# ---------------------------------------------------------------------------
# vg deconstruct
# ---------------------------------------------------------------------------

_DECONSTRUCT = "\n".join([
    "##fileformat=VCFv4.2",
    '##INFO=<ID=LV,Number=1,Type=Integer,Description="Level">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ">s1>s3\t100\t.\tA\tG\t.\t.\tLV=0;AT=>s1>s2>s3,>s1>s3",
    ">s4>s6\t250\t.\tACG\tA,ACGT\t.\t.\tLV=1;PS=>s1>s3;AT=>s4>s5>s6",
]) + "\n"


class TestReadDeconstruct:
    def test_parses_sites(self, tmp_path):
        p = tmp_path / "d.vcf"
        p.write_text(_DECONSTRUCT, encoding="utf-8")
        sites = list(read_deconstruct_sites(p))
        assert len(sites) == 2
        s0 = sites[0]
        assert s0.pos == 99                 # 1-based 100 -> 0-based 99
        assert s0.ref == "A" and s0.alts == ["G"]
        assert s0.level == 0 and s0.is_top_level
        assert s0.n_alleles == 2
        assert s0.allele_traversals == [">s1>s2>s3", ">s1>s3"]
        s1 = sites[1]
        assert s1.n_alleles == 3            # REF + 2 ALTs
        assert s1.level == 1 and not s1.is_top_level
        assert s1.parent_snarl == ">s1>s3"

    def test_top_level_only_filter(self, tmp_path):
        p = tmp_path / "d.vcf"
        p.write_text(_DECONSTRUCT, encoding="utf-8")
        sites = list(read_deconstruct_sites(p, top_level_only=True))
        assert len(sites) == 1
        assert sites[0].level == 0
