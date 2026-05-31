"""Unit tests for src/privy/structure/karyotype.py."""

from __future__ import annotations

from privy.structure import (
    bin_chromosome,
    find_centromere,
    find_telomeres,
    kmer_blocks,
)
from privy.structure.karyotype import PLANT_TELOMERE_MOTIFS

TELO = "CCCTAAA"
SAT = "ATATATGG"   # toy satellite motif


def _chromosome() -> str:
    """A toy chromosome: telomere arrays at both ends, a satellite array in the middle."""
    left_telo = TELO * 60                  # ~420 bp telomere
    right_telo = TELO * 60
    filler = "ACGTACGTAC" * 500            # 5000 bp arm filler (no motifs)
    satellite = SAT * 400                  # ~3200 bp centromeric satellite
    return left_telo + filler + satellite + filler + right_telo


class TestKmerBlocks:
    def test_finds_and_merges(self):
        seq = ("N" * 100) + (TELO * 20) + ("N" * 100)
        blocks = kmer_blocks(seq, [TELO], max_gap=10)
        assert len(blocks) == 1
        start, end = blocks[0]
        assert start == 100
        assert end == 100 + 7 * 20

    def test_no_motif_empty(self):
        assert kmer_blocks("ACGTACGT", ["TTTTT"]) == []

    def test_gap_splits_blocks(self):
        seq = (TELO * 5) + ("N" * 500) + (TELO * 5)
        blocks = kmer_blocks(seq, [TELO], max_gap=50)
        assert len(blocks) == 2


class TestTelomeres:
    def test_both_ends_capped(self):
        result = find_telomeres(_chromosome(), end_window=1000, min_array_len=100)
        assert result.five_prime_capped is True
        assert result.three_prime_capped is True
        assert len(result.blocks) >= 2

    def test_no_telomere_uncapped(self):
        seq = "ACGTACGTAC" * 1000
        result = find_telomeres(seq)
        assert result.five_prime_capped is False
        assert result.three_prime_capped is False
        assert result.blocks == []


class TestCentromere:
    def test_finds_satellite_array(self):
        centro = find_centromere(_chromosome(), [SAT], min_array_len=500)
        assert centro is not None
        start, end = centro
        assert end - start >= 3000   # the ~3200 bp satellite array

    def test_none_without_satellite(self):
        assert find_centromere(_chromosome(), ["GGGGGGGG"]) is None


class TestBinChromosome:
    def test_full_partition(self):
        seq = _chromosome()
        bins = bin_chromosome(
            seq, satellite_motifs=[SAT], end_window=1000,
            telomere_min_len=100, centromere_min_len=500, pericentromere_flank=800,
        )
        types = [b.bin_type for b in bins]
        assert "telomere" in types
        assert "centromere" in types
        assert "pericentromere" in types
        assert "arm" in types
        # contiguous cover of the whole sequence
        assert bins[0].start == 0
        assert bins[-1].end == len(seq)
        for a, b in zip(bins, bins[1:], strict=False):
            assert a.end == b.start

    def test_empty_sequence(self):
        assert bin_chromosome("") == []

    def test_default_motifs_are_plant(self):
        assert PLANT_TELOMERE_MOTIFS == ("CCCTAAA", "TTTAGGG")
