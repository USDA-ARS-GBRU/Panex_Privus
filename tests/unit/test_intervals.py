"""Unit tests for privy.core.intervals — locus merging and overlap logic."""

import pytest

from privy.core.intervals import (
    iter_contig_chunks,
    merge_loci_to_regions,
    reciprocal_overlap,
)
from privy.core.locus import Locus, LocusType, PrimarySource


def snp(locus_id: str, contig: str, pos: int) -> Locus:
    """Convenience: single-base SNP locus."""
    return Locus(
        locus_id=locus_id,
        contig=contig,
        start=pos,
        end=pos + 1,
        locus_type=LocusType.SNP,
        primary_source=PrimarySource.VCF,
    )


def interval(locus_id: str, contig: str, start: int, end: int,
             ltype: LocusType = LocusType.SNP) -> Locus:
    return Locus(
        locus_id=locus_id,
        contig=contig,
        start=start,
        end=end,
        locus_type=ltype,
        primary_source=PrimarySource.VCF,
    )


# ---------------------------------------------------------------------------
# merge_loci_to_regions
# ---------------------------------------------------------------------------

class TestMergeLociToRegions:
    def test_empty_input(self) -> None:
        assert merge_loci_to_regions([]) == []

    def test_single_locus_becomes_region(self) -> None:
        loci = [snp("L1", "chr1", 100)]
        result = merge_loci_to_regions(loci)
        assert len(result) == 1
        assert result[0].locus_type == LocusType.REGION

    def test_two_adjacent_loci_merged(self) -> None:
        loci = [snp("L1", "chr1", 100), snp("L2", "chr1", 101)]
        result = merge_loci_to_regions(loci, merge_distance=0)
        assert len(result) == 1
        assert result[0].start == 100
        assert result[0].end == 102

    def test_gap_within_distance_merged(self) -> None:
        loci = [snp("L1", "chr1", 100), snp("L2", "chr1", 200)]
        result = merge_loci_to_regions(loci, merge_distance=100)
        assert len(result) == 1
        assert result[0].start == 100
        assert result[0].end == 201

    def test_gap_exceeds_distance_not_merged(self) -> None:
        loci = [snp("L1", "chr1", 100), snp("L2", "chr1", 500)]
        result = merge_loci_to_regions(loci, merge_distance=100)
        assert len(result) == 2

    def test_different_contigs_not_merged(self) -> None:
        loci = [snp("L1", "chr1", 100), snp("L2", "chr2", 100)]
        result = merge_loci_to_regions(loci, merge_distance=10_000)
        assert len(result) == 2
        contigs = {r.contig for r in result}
        assert contigs == {"chr1", "chr2"}

    def test_unsorted_input_sorted_internally(self) -> None:
        loci = [snp("L2", "chr1", 200), snp("L1", "chr1", 100)]
        result = merge_loci_to_regions(loci, merge_distance=200)
        assert len(result) == 1
        assert result[0].start == 100

    def test_three_loci_two_groups(self) -> None:
        loci = [
            snp("L1", "chr1", 100),
            snp("L2", "chr1", 110),
            snp("L3", "chr1", 500),
        ]
        result = merge_loci_to_regions(loci, merge_distance=50)
        assert len(result) == 2

    def test_region_id_prefix(self) -> None:
        loci = [snp("L1", "chr1", 100)]
        result = merge_loci_to_regions(loci, region_id_prefix="PPX")
        assert result[0].locus_id.startswith("PPX")

    def test_n_constituent_loci_metadata(self) -> None:
        loci = [snp("L1", "chr1", 100), snp("L2", "chr1", 105)]
        result = merge_loci_to_regions(loci, merge_distance=10)
        assert result[0].metadata.get("n_constituent_loci") == "2"

    def test_same_variant_class_only_prevents_merge(self) -> None:
        loci = [
            interval("L1", "chr1", 100, 101, LocusType.SNP),
            interval("L2", "chr1", 105, 110, LocusType.INDEL),
        ]
        result = merge_loci_to_regions(
            loci, merge_distance=100, same_variant_class_only=True
        )
        assert len(result) == 2

    def test_same_variant_class_allows_merge(self) -> None:
        loci = [
            interval("L1", "chr1", 100, 101, LocusType.SNP),
            interval("L2", "chr1", 105, 106, LocusType.SNP),
        ]
        result = merge_loci_to_regions(
            loci, merge_distance=100, same_variant_class_only=True
        )
        assert len(result) == 1

    def test_overlapping_loci_merged(self) -> None:
        loci = [
            interval("L1", "chr1", 100, 200),
            interval("L2", "chr1", 150, 300),
        ]
        result = merge_loci_to_regions(loci, merge_distance=0)
        assert len(result) == 1
        assert result[0].end == 300

    def test_source_ids_accumulated(self) -> None:
        a = Locus("L1", "chr1", 100, 101, LocusType.SNP, PrimarySource.VCF,
                  source_ids=["rs1"])
        b = Locus("L2", "chr1", 105, 106, LocusType.SNP, PrimarySource.VCF,
                  source_ids=["rs2"])
        result = merge_loci_to_regions([a, b], merge_distance=10)
        assert "rs1" in result[0].source_ids
        assert "rs2" in result[0].source_ids


# ---------------------------------------------------------------------------
# reciprocal_overlap
# ---------------------------------------------------------------------------

class TestReciprocalOverlap:
    def test_identical_loci(self) -> None:
        a = interval("A", "chr1", 100, 200)
        assert reciprocal_overlap(a, a) == pytest.approx(1.0)

    def test_no_overlap(self) -> None:
        a = interval("A", "chr1", 100, 200)
        b = interval("B", "chr1", 300, 400)
        assert reciprocal_overlap(a, b) == pytest.approx(0.0)

    def test_50_pct_overlap(self) -> None:
        a = interval("A", "chr1", 0, 100)
        b = interval("B", "chr1", 50, 150)
        # intersection = 50, union = 150 → 50/150 ≈ 0.333
        assert reciprocal_overlap(a, b) == pytest.approx(50 / 150)

    def test_contained_locus(self) -> None:
        outer = interval("O", "chr1", 0, 200)
        inner = interval("I", "chr1", 50, 100)
        # intersection = 50, union = 200 → 0.25
        assert reciprocal_overlap(outer, inner) == pytest.approx(50 / 200)

    def test_different_contigs(self) -> None:
        a = interval("A", "chr1", 0, 100)
        b = interval("B", "chr2", 0, 100)
        assert reciprocal_overlap(a, b) == 0.0


# ---------------------------------------------------------------------------
# iter_contig_chunks
# ---------------------------------------------------------------------------

class TestIterContigChunks:
    def test_empty(self) -> None:
        assert list(iter_contig_chunks([], chunk_size=1000)) == []

    def test_single_locus(self) -> None:
        loci = [snp("L1", "chr1", 100)]
        chunks = list(iter_contig_chunks(loci, chunk_size=1000))
        assert len(chunks) == 1
        assert chunks[0] == loci

    def test_loci_in_different_windows(self) -> None:
        loci = [
            snp("L1", "chr1", 100),
            snp("L2", "chr1", 2000),  # different 1kb chunk
        ]
        chunks = list(iter_contig_chunks(loci, chunk_size=1000))
        assert len(chunks) == 2

    def test_contig_boundary_splits(self) -> None:
        loci = [
            snp("L1", "chr1", 100),
            snp("L2", "chr2", 100),
        ]
        chunks = list(iter_contig_chunks(loci, chunk_size=10_000))
        assert len(chunks) == 2
