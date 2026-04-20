"""Unit tests for src/privy/io/bam.py."""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from privy.io.bam import (
    get_bam_sample_name,
    load_bam_manifest,
    query_allele_counts_at_locus,
    query_position_depth,
    validate_bam_index,
)


# ---------------------------------------------------------------------------
# validate_bam_index
# ---------------------------------------------------------------------------

class TestValidateBamIndex:
    def test_raises_when_no_index(self, tmp_path: Path) -> None:
        bam = tmp_path / "sample.bam"
        bam.write_bytes(b"")
        with pytest.raises(FileNotFoundError, match="BAM index not found"):
            validate_bam_index(bam)

    def test_passes_with_bai_index(self, tmp_path: Path) -> None:
        bam = tmp_path / "sample.bam"
        bai = tmp_path / "sample.bam.bai"
        bam.write_bytes(b"")
        bai.write_bytes(b"")
        validate_bam_index(bam)  # no exception

    def test_passes_with_csi_index(self, tmp_path: Path) -> None:
        bam = tmp_path / "sample.bam"
        csi = tmp_path / "sample.bam.csi"
        bam.write_bytes(b"")
        csi.write_bytes(b"")
        validate_bam_index(bam)  # no exception

    def test_error_message_includes_path(self, tmp_path: Path) -> None:
        bam = tmp_path / "mysample.bam"
        bam.write_bytes(b"")
        with pytest.raises(FileNotFoundError, match="mysample.bam"):
            validate_bam_index(bam)

    def test_error_mentions_samtools(self, tmp_path: Path) -> None:
        bam = tmp_path / "sample.bam"
        bam.write_bytes(b"")
        with pytest.raises(FileNotFoundError, match="samtools index"):
            validate_bam_index(bam)


# ---------------------------------------------------------------------------
# get_bam_sample_name
# ---------------------------------------------------------------------------

class TestGetBamSampleName:
    def test_returns_sm_tag(self, bam_target_t1: Path) -> None:
        name = get_bam_sample_name(bam_target_t1)
        assert name == "T1"

    def test_returns_sm_tag_offtarget(self, bam_offtarget_o1: Path) -> None:
        name = get_bam_sample_name(bam_offtarget_o1)
        assert name == "O1"

    def test_returns_none_for_no_rg(self, tmp_path: Path) -> None:
        import pysam  # noqa: PLC0415

        bam_path = tmp_path / "norg.bam"
        header = pysam.AlignmentHeader.from_dict({
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": "chr1", "LN": 1000}],
        })
        with pysam.AlignmentFile(str(bam_path), "wb", header=header):
            pass
        pysam.index(str(bam_path))
        assert get_bam_sample_name(bam_path) is None

    def test_returns_none_for_rg_without_sm(self, tmp_path: Path) -> None:
        import pysam  # noqa: PLC0415

        bam_path = tmp_path / "nosm.bam"
        header = pysam.AlignmentHeader.from_dict({
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": "chr1", "LN": 1000}],
            "RG": [{"ID": "lane1"}],
        })
        with pysam.AlignmentFile(str(bam_path), "wb", header=header):
            pass
        pysam.index(str(bam_path))
        assert get_bam_sample_name(bam_path) is None


# ---------------------------------------------------------------------------
# load_bam_manifest
# ---------------------------------------------------------------------------

class TestLoadBamManifest:
    def _write_manifest(self, path: Path, rows: list[dict]) -> Path:
        fieldnames = list(rows[0].keys()) if rows else ["bam_path", "sample_id"]
        with open(path, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        return path

    def test_parses_valid_manifest(self, tmp_path: Path) -> None:
        manifest = self._write_manifest(tmp_path / "manifest.tsv", [
            {"bam_path": "/data/T1.bam", "sample_id": "T1"},
            {"bam_path": "/data/O1.bam", "sample_id": "O1"},
        ])
        rows = load_bam_manifest(manifest)
        assert len(rows) == 2
        assert rows[0]["sample_id"] == "T1"
        assert rows[1]["bam_path"] == "/data/O1.bam"

    def test_preserves_optional_group_column(self, tmp_path: Path) -> None:
        manifest = self._write_manifest(tmp_path / "manifest.tsv", [
            {"bam_path": "/data/T1.bam", "sample_id": "T1", "group": "focal"},
        ])
        rows = load_bam_manifest(manifest)
        assert rows[0]["group"] == "focal"

    def test_skips_comment_lines(self, tmp_path: Path) -> None:
        path = tmp_path / "manifest.tsv"
        path.write_text(
            "# This is a comment\nbam_path\tsample_id\n/d/T1.bam\tT1\n"
        )
        rows = load_bam_manifest(path)
        assert len(rows) == 1
        assert rows[0]["sample_id"] == "T1"

    def test_raises_for_missing_file(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError, match="not found"):
            load_bam_manifest(tmp_path / "missing.tsv")

    def test_raises_for_missing_bam_path_column(self, tmp_path: Path) -> None:
        path = tmp_path / "bad.tsv"
        path.write_text("sample_id\nT1\n")
        with pytest.raises(ValueError, match="bam_path"):
            load_bam_manifest(path)

    def test_raises_for_missing_sample_id_column(self, tmp_path: Path) -> None:
        path = tmp_path / "bad.tsv"
        path.write_text("bam_path\n/data/T1.bam\n")
        with pytest.raises(ValueError, match="sample_id"):
            load_bam_manifest(path)

    def test_empty_manifest_returns_empty_list(self, tmp_path: Path) -> None:
        path = tmp_path / "empty.tsv"
        path.write_text("bam_path\tsample_id\n")
        rows = load_bam_manifest(path)
        assert rows == []


# ---------------------------------------------------------------------------
# query_position_depth
# ---------------------------------------------------------------------------

class TestQueryPositionDepth:
    def test_returns_nonzero_depth_at_covered_positions(
        self, bam_target_t1: Path
    ) -> None:
        # 12 reads starting at 85, length 30 → cover [85, 115)
        depths = query_position_depth(bam_target_t1, "chr1", 85, 115)
        assert len(depths) == 30
        assert all(d == 12 for d in depths)

    def test_returns_zeros_for_empty_region(self, bam_target_t1: Path) -> None:
        depths = query_position_depth(bam_target_t1, "chr1", 0, 50)
        assert len(depths) == 50
        assert all(d == 0 for d in depths)

    def test_returns_zeros_for_unknown_contig(self, bam_target_t1: Path) -> None:
        depths = query_position_depth(bam_target_t1, "chrX", 0, 10)
        assert depths == [0] * 10

    def test_returns_empty_for_zero_width_region(self, bam_target_t1: Path) -> None:
        depths = query_position_depth(bam_target_t1, "chr1", 100, 100)
        assert depths == []

    def test_mapq_filter_excludes_low_mapq_reads(self, tmp_path: Path) -> None:
        import pysam  # noqa: PLC0415

        bam_path = tmp_path / "lowmapq.bam"
        header = pysam.AlignmentHeader.from_dict({
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": "chr1", "LN": 10000}],
            "RG": [{"ID": "X", "SM": "X"}],
        })
        with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
            for i in range(10):
                seg = pysam.AlignedSegment(header)
                seg.query_name = f"r{i}"
                seg.query_sequence = "A" * 30
                seg.flag = 0
                seg.reference_id = 0
                seg.reference_start = 85
                seg.mapping_quality = 5  # below default min_mapq=20
                seg.cigar = [(0, 30)]
                seg.query_qualities = pysam.qualitystring_to_array("I" * 30)
                seg.set_tag("RG", "X")
                bam.write(seg)
        pysam.index(str(bam_path))
        depths = query_position_depth(bam_path, "chr1", 85, 115, min_mapq=20)
        assert all(d == 0 for d in depths)

    def test_mapq_filter_passes_high_mapq_reads(self, bam_target_t1: Path) -> None:
        depths = query_position_depth(bam_target_t1, "chr1", 85, 115, min_mapq=20)
        assert all(d == 12 for d in depths)


# ---------------------------------------------------------------------------
# query_allele_counts_at_locus
# ---------------------------------------------------------------------------

class TestQueryAlleleCountsAtLocus:
    def test_snp_counts_alt_reads(self, bam_target_t1: Path) -> None:
        ref_c, alt_c, other_c = query_allele_counts_at_locus(
            bam_target_t1, "chr1", 99, "A", "T"
        )
        assert alt_c == 12
        assert ref_c == 0
        assert other_c == 0

    def test_snp_counts_ref_reads(self, bam_offtarget_o1: Path) -> None:
        ref_c, alt_c, other_c = query_allele_counts_at_locus(
            bam_offtarget_o1, "chr1", 99, "A", "T"
        )
        assert ref_c == 12
        assert alt_c == 0
        assert other_c == 0

    def test_snp_sum_equals_depth(self, bam_target_t1: Path) -> None:
        ref_c, alt_c, other_c = query_allele_counts_at_locus(
            bam_target_t1, "chr1", 99, "A", "T"
        )
        assert ref_c + alt_c + other_c == 12

    def test_indel_returns_depth_only(self, bam_target_t1: Path) -> None:
        ref_c, alt_c, other_c = query_allele_counts_at_locus(
            bam_target_t1, "chr1", 98, "AGG", "A"
        )
        assert alt_c == 0
        assert other_c == 0
        assert ref_c > 0  # depth-only mode: total depth in ref_count

    def test_returns_zeros_for_uncovered_position(
        self, bam_target_t1: Path
    ) -> None:
        ref_c, alt_c, other_c = query_allele_counts_at_locus(
            bam_target_t1, "chr1", 0, "A", "T"
        )
        assert ref_c == 0 and alt_c == 0 and other_c == 0

    def test_returns_zeros_for_unknown_contig(
        self, bam_target_t1: Path
    ) -> None:
        ref_c, alt_c, other_c = query_allele_counts_at_locus(
            bam_target_t1, "chrX", 99, "A", "T"
        )
        assert ref_c == 0 and alt_c == 0 and other_c == 0

    def test_case_insensitive_allele_matching(self, bam_target_t1: Path) -> None:
        ref_c_lower, alt_c_lower, _ = query_allele_counts_at_locus(
            bam_target_t1, "chr1", 99, "a", "t"
        )
        assert alt_c_lower == 12

    def test_baseq_filter_at_zero_passes_all(self, bam_target_t1: Path) -> None:
        ref_c, alt_c, _ = query_allele_counts_at_locus(
            bam_target_t1, "chr1", 99, "A", "T", min_baseq=0
        )
        assert alt_c == 12

    def test_high_baseq_excludes_reads(self, tmp_path: Path) -> None:
        import pysam  # noqa: PLC0415

        bam_path = tmp_path / "lowbq.bam"
        header = pysam.AlignmentHeader.from_dict({
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": "chr1", "LN": 10000}],
            "RG": [{"ID": "X", "SM": "X"}],
        })
        seq = "A" * 14 + "T" + "A" * 15
        with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
            for i in range(10):
                seg = pysam.AlignedSegment(header)
                seg.query_name = f"r{i}"
                seg.query_sequence = seq
                seg.flag = 0
                seg.reference_id = 0
                seg.reference_start = 85
                seg.mapping_quality = 60
                seg.cigar = [(0, 30)]
                # Base quality = 10 (low) at offset 14
                quals = [40] * 30
                quals[14] = 10
                seg.query_qualities = quals
                seg.set_tag("RG", "X")
                bam.write(seg)
        pysam.index(str(bam_path))
        _, alt_c, _ = query_allele_counts_at_locus(
            bam_path, "chr1", 99, "A", "T", min_baseq=20
        )
        assert alt_c == 0
