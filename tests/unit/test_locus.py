"""Unit tests for privy.core.locus.Locus."""

import pytest

from privy.core.locus import Locus, LocusType, PrimarySource


def make_locus(
    locus_id: str = "L001",
    contig: str = "chr1",
    start: int = 100,
    end: int = 101,
    locus_type: LocusType = LocusType.SNP,
) -> Locus:
    return Locus(
        locus_id=locus_id,
        contig=contig,
        start=start,
        end=end,
        locus_type=locus_type,
        primary_source=PrimarySource.VCF,
    )


class TestLocusConstruction:
    def test_basic(self) -> None:
        loc = make_locus(start=0, end=1)
        assert loc.length == 1

    def test_negative_start_raises(self) -> None:
        with pytest.raises(ValueError, match="start must be >= 0"):
            make_locus(start=-1, end=0)

    def test_end_before_start_raises(self) -> None:
        with pytest.raises(ValueError, match="end must be >= start"):
            make_locus(start=10, end=5)

    def test_zero_length_is_valid(self) -> None:
        loc = make_locus(start=100, end=100)
        assert loc.length == 0

    def test_length(self) -> None:
        loc = make_locus(start=100, end=200)
        assert loc.length == 100

    def test_region_string(self) -> None:
        loc = make_locus(start=99, end=200)  # 0-based 99 → 1-based 100
        assert loc.region_string == "chr1:100-200"

    def test_repr(self) -> None:
        loc = make_locus(locus_id="PPX000001", start=100, end=200)
        assert "PPX000001" in repr(loc)
        assert "chr1:100-200" in repr(loc)


class TestLocusOverlap:
    def test_overlapping(self) -> None:
        a = make_locus(start=100, end=200)
        b = make_locus(start=150, end=250)
        assert a.overlaps(b)
        assert b.overlaps(a)

    def test_adjacent_not_overlapping(self) -> None:
        a = make_locus(start=100, end=200)
        b = make_locus(start=200, end=300)
        assert not a.overlaps(b)

    def test_contained(self) -> None:
        outer = make_locus(start=100, end=300)
        inner = make_locus(start=150, end=200)
        assert outer.overlaps(inner)
        assert outer.contains(inner)
        assert not inner.contains(outer)

    def test_different_contig(self) -> None:
        a = make_locus(contig="chr1", start=100, end=200)
        b = make_locus(contig="chr2", start=100, end=200)
        assert not a.overlaps(b)
        assert not a.contains(b)


class TestLocusDistance:
    def test_non_overlapping(self) -> None:
        a = make_locus(start=100, end=200)
        b = make_locus(start=300, end=400)
        assert a.distance_to(b) == 100
        assert b.distance_to(a) == 100

    def test_overlapping_returns_zero(self) -> None:
        a = make_locus(start=100, end=200)
        b = make_locus(start=150, end=250)
        assert a.distance_to(b) == 0

    def test_different_contig_returns_none(self) -> None:
        a = make_locus(contig="chr1", start=100, end=200)
        b = make_locus(contig="chr2", start=100, end=200)
        assert a.distance_to(b) is None

    def test_adjacent_distance_zero(self) -> None:
        a = make_locus(start=100, end=200)
        b = make_locus(start=200, end=300)
        assert a.distance_to(b) == 0


class TestLocusMerge:
    def test_merge_two_loci(self) -> None:
        a = make_locus(locus_id="L1", start=100, end=200)
        b = make_locus(locus_id="L2", start=300, end=400)
        merged = a.merge_with(b)
        assert merged.start == 100
        assert merged.end == 400
        assert merged.locus_type == LocusType.REGION
        assert merged.contig == "chr1"

    def test_merge_with_custom_id(self) -> None:
        a = make_locus(locus_id="L1", start=100, end=200)
        b = make_locus(locus_id="L2", start=300, end=400)
        merged = a.merge_with(b, merged_id="CUSTOM001")
        assert merged.locus_id == "CUSTOM001"

    def test_merge_different_contigs_raises(self) -> None:
        a = make_locus(contig="chr1", start=100, end=200)
        b = make_locus(contig="chr2", start=100, end=200)
        with pytest.raises(ValueError, match="different contigs"):
            a.merge_with(b)

    def test_merge_accumulates_source_ids(self) -> None:
        a = Locus("L1", "chr1", 100, 200, LocusType.SNP, PrimarySource.VCF, source_ids=["rs1"])
        b = Locus("L2", "chr1", 300, 400, LocusType.SNP, PrimarySource.VCF, source_ids=["rs2"])
        merged = a.merge_with(b)
        assert "rs1" in merged.source_ids
        assert "rs2" in merged.source_ids
