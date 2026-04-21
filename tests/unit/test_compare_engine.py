"""Unit tests for src/privy/backends/compare.py."""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from privy.backends.compare import (
    HitsRow,
    classify_match,
    compute_comparison_score,
    find_best_match,
    infer_source_label,
    is_state_compatible,
    load_hits_tsv,
    reciprocal_overlap_rows,
)
from privy.core.config import CompareConfig, default_config
from privy.core.evidence import MatchClass


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_row(
    locus_id: str = "PPX000001",
    contig: str = "chr1",
    start: int = 100,
    end: int = 200,
    strictness_class: str = "strict_complete",
    final_score: float = 1.0,
) -> HitsRow:
    return HitsRow(
        locus_id=locus_id,
        contig=contig,
        start=start,
        end=end,
        variant_type="snp",
        allele_key="chr1:100:A:T",
        strictness_class=strictness_class,
        final_score=final_score,
    )


def _make_hits_tsv(tmp_path: Path, rows: list[dict[str, str]]) -> Path:
    path = tmp_path / "hits.tsv"
    if not rows:
        path.write_text(
            "locus_id\tcontig\tstart\tend\tvariant_type\tallele_key\t"
            "target_support_n\ttarget_total_n\tofftarget_support_n\t"
            "offtarget_total_n\ttarget_missing_n\tofftarget_missing_n\t"
            "strictness_class\tdiscovery_score\tsupport_score\t"
            "penalty_score\tfinal_score\n",
            encoding="utf-8",
        )
        return path
    fieldnames = list(rows[0].keys())
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return path


_HITS_ROW_DEFAULTS: dict[str, str] = {
    "variant_type": "snp",
    "allele_key": "chr1:100:A:T",
    "target_support_n": "3",
    "target_total_n": "3",
    "offtarget_support_n": "0",
    "offtarget_total_n": "2",
    "target_missing_n": "0",
    "offtarget_missing_n": "0",
    "discovery_score": "1.0",
    "support_score": "0.0",
    "penalty_score": "0.0",
    "final_score": "1.0",
}


def _row_dict(**kwargs: str) -> dict[str, str]:
    d = dict(_HITS_ROW_DEFAULTS)
    d.update(kwargs)
    return d


# ---------------------------------------------------------------------------
# TestLoadHitsTsv
# ---------------------------------------------------------------------------

class TestLoadHitsTsv:
    def test_loads_valid_file(self, tmp_path: Path) -> None:
        path = _make_hits_tsv(tmp_path, [
            _row_dict(locus_id="PPX000001", contig="chr1", start="100", end="200",
                      strictness_class="strict_complete"),
        ])
        rows = load_hits_tsv(path)
        assert len(rows) == 1
        assert rows[0].locus_id == "PPX000001"
        assert rows[0].start == 100
        assert rows[0].end == 200
        assert rows[0].final_score == 1.0

    def test_empty_file_returns_empty_list(self, tmp_path: Path) -> None:
        path = _make_hits_tsv(tmp_path, [])
        rows = load_hits_tsv(path)
        assert rows == []

    def test_missing_file_raises(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            load_hits_tsv(tmp_path / "nonexistent.tsv")

    def test_malformed_start_raises(self, tmp_path: Path) -> None:
        path = _make_hits_tsv(tmp_path, [
            _row_dict(locus_id="PPX000001", contig="chr1", start="notanint", end="200",
                      strictness_class="strict_complete"),
        ])
        with pytest.raises(ValueError, match="Malformed hits.tsv row"):
            load_hits_tsv(path)

    def test_multiple_rows_parsed(self, tmp_path: Path) -> None:
        path = _make_hits_tsv(tmp_path, [
            _row_dict(locus_id="PPX000001", contig="chr1", start="100", end="200",
                      strictness_class="strict_complete"),
            _row_dict(locus_id="PPX000002", contig="chr2", start="500", end="600",
                      strictness_class="strict_target_missing", final_score="0.7"),
        ])
        rows = load_hits_tsv(path)
        assert len(rows) == 2
        assert rows[1].contig == "chr2"
        assert rows[1].final_score == pytest.approx(0.7)


# ---------------------------------------------------------------------------
# TestInferSourceLabel
# ---------------------------------------------------------------------------

class TestInferSourceLabel:
    def test_explicit_label_wins(self) -> None:
        rows = [_make_row(locus_id="GPX000001")]
        assert infer_source_label(rows, "my_source") == "my_source"

    def test_gpx_prefix_infers_gfa(self) -> None:
        rows = [_make_row(locus_id="GPX000001")]
        assert infer_source_label(rows, None) == "gfa"

    def test_ppx_prefix_infers_vcf(self) -> None:
        rows = [_make_row(locus_id="PPX000001")]
        assert infer_source_label(rows, None) == "vcf"

    def test_empty_rows_defaults_vcf(self) -> None:
        assert infer_source_label([], None) == "vcf"

    def test_unknown_prefix_defaults_vcf(self) -> None:
        rows = [_make_row(locus_id="XYZ000001")]
        assert infer_source_label(rows, None) == "vcf"


# ---------------------------------------------------------------------------
# TestReciprocalOverlapRows
# ---------------------------------------------------------------------------

class TestReciprocalOverlapRows:
    def test_perfect_overlap(self) -> None:
        a = _make_row(start=100, end=200)
        b = _make_row(start=100, end=200)
        assert reciprocal_overlap_rows(a, b) == pytest.approx(1.0)

    def test_no_overlap(self) -> None:
        a = _make_row(start=100, end=200)
        b = _make_row(start=300, end=400)
        assert reciprocal_overlap_rows(a, b) == pytest.approx(0.0)

    def test_different_contig(self) -> None:
        a = _make_row(contig="chr1", start=100, end=200)
        b = _make_row(contig="chr2", start=100, end=200)
        assert reciprocal_overlap_rows(a, b) == pytest.approx(0.0)

    def test_partial_overlap(self) -> None:
        # a=[100,200), b=[150,250) → intersection=50, union=150
        a = _make_row(start=100, end=200)
        b = _make_row(start=150, end=250)
        assert reciprocal_overlap_rows(a, b) == pytest.approx(50 / 150)

    def test_contained_locus(self) -> None:
        # a=[100,400), b=[200,300) → intersection=100, union=300
        a = _make_row(start=100, end=400)
        b = _make_row(start=200, end=300)
        assert reciprocal_overlap_rows(a, b) == pytest.approx(100 / 300)

    def test_adjacent_no_overlap(self) -> None:
        a = _make_row(start=100, end=200)
        b = _make_row(start=200, end=300)
        assert reciprocal_overlap_rows(a, b) == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# TestIsStateCompatible
# ---------------------------------------------------------------------------

class TestIsStateCompatible:
    def test_a_contradicted_returns_false(self) -> None:
        assert is_state_compatible("contradicted", "strict_complete", require=False) is False

    def test_b_contradicted_returns_false(self) -> None:
        assert is_state_compatible("strict_complete", "contradicted", require=False) is False

    def test_both_contradicted_returns_false(self) -> None:
        assert is_state_compatible("contradicted", "contradicted", require=False) is False

    def test_both_strict_no_require(self) -> None:
        assert is_state_compatible("strict_complete", "strict_target_missing", require=False) is True

    def test_both_strict_require(self) -> None:
        assert is_state_compatible("strict_complete", "strict_target_missing", require=True) is True

    def test_strict_vs_relaxed_no_require(self) -> None:
        assert is_state_compatible("strict_complete", "relaxed_threshold", require=False) is True

    def test_strict_vs_relaxed_require(self) -> None:
        assert is_state_compatible("strict_complete", "relaxed_threshold", require=True) is False

    def test_both_relaxed_no_require(self) -> None:
        assert is_state_compatible("relaxed_threshold", "relaxed_threshold", require=False) is True

    def test_both_relaxed_require(self) -> None:
        assert is_state_compatible("relaxed_threshold", "relaxed_threshold", require=True) is True


# ---------------------------------------------------------------------------
# TestClassifyMatch
# ---------------------------------------------------------------------------

class TestClassifyMatch:
    def _cfg(self, **kwargs: object) -> CompareConfig:
        return CompareConfig(**kwargs)  # type: ignore[arg-type]

    def test_contradicted_a(self) -> None:
        cfg = self._cfg(min_reciprocal_overlap=0.5, require_state_compatibility=False)
        mc = classify_match(0.8, True, "contradicted", "strict_complete", cfg)
        assert mc == MatchClass.CONTRADICTED

    def test_contradicted_b(self) -> None:
        cfg = self._cfg(min_reciprocal_overlap=0.5, require_state_compatibility=False)
        mc = classify_match(0.8, True, "strict_complete", "contradicted", cfg)
        assert mc == MatchClass.CONTRADICTED

    def test_zero_overlap_source_specific(self) -> None:
        cfg = self._cfg(min_reciprocal_overlap=0.5, require_state_compatibility=False)
        mc = classify_match(0.0, True, "strict_complete", "strict_complete", cfg)
        assert mc == MatchClass.SOURCE_SPECIFIC

    def test_overlap_state_compat_supported(self) -> None:
        cfg = self._cfg(min_reciprocal_overlap=0.5, require_state_compatibility=False)
        mc = classify_match(0.8, True, "strict_complete", "strict_complete", cfg)
        assert mc == MatchClass.SUPPORTED

    def test_overlap_state_incompat_no_require(self) -> None:
        cfg = self._cfg(min_reciprocal_overlap=0.5, require_state_compatibility=False)
        mc = classify_match(0.8, False, "strict_complete", "relaxed_threshold", cfg)
        assert mc == MatchClass.PARTIALLY_SUPPORTED

    def test_overlap_state_incompat_require(self) -> None:
        cfg = self._cfg(min_reciprocal_overlap=0.5, require_state_compatibility=True)
        mc = classify_match(0.8, False, "strict_complete", "relaxed_threshold", cfg)
        assert mc == MatchClass.PARTIALLY_SUPPORTED

    def test_overlap_state_compat_require_supported(self) -> None:
        cfg = self._cfg(min_reciprocal_overlap=0.5, require_state_compatibility=True)
        mc = classify_match(0.8, True, "strict_complete", "strict_target_missing", cfg)
        assert mc == MatchClass.SUPPORTED


# ---------------------------------------------------------------------------
# TestComputeComparisonScore
# ---------------------------------------------------------------------------

class TestComputeComparisonScore:
    def test_supported_perfect_overlap(self) -> None:
        score = compute_comparison_score(MatchClass.SUPPORTED, 1.0)
        assert score == pytest.approx(1.0)

    def test_supported_half_overlap(self) -> None:
        score = compute_comparison_score(MatchClass.SUPPORTED, 0.5)
        assert score == pytest.approx(0.5)

    def test_partially_supported_half_overlap(self) -> None:
        score = compute_comparison_score(MatchClass.PARTIALLY_SUPPORTED, 0.5)
        assert score == pytest.approx(0.25)

    def test_contradicted_zero(self) -> None:
        assert compute_comparison_score(MatchClass.CONTRADICTED, 0.9) == pytest.approx(0.0)

    def test_source_specific_fixed(self) -> None:
        assert compute_comparison_score(MatchClass.SOURCE_SPECIFIC, 0.0) == pytest.approx(0.3)

    def test_uninformative_fixed(self) -> None:
        assert compute_comparison_score(MatchClass.UNINFORMATIVE, 0.0) == pytest.approx(0.1)

    def test_missing_data_zero(self) -> None:
        assert compute_comparison_score(MatchClass.MISSING_DATA, 0.0) == pytest.approx(0.0)

    def test_low_overlap_clamps_to_minimum(self) -> None:
        score = compute_comparison_score(MatchClass.SUPPORTED, 0.0)
        assert score == pytest.approx(1.0 * 0.01)


# ---------------------------------------------------------------------------
# TestFindBestMatch
# ---------------------------------------------------------------------------

class TestFindBestMatch:
    def _cfg(self, min_overlap: float = 0.5, tol: int = 0) -> CompareConfig:
        return CompareConfig(min_reciprocal_overlap=min_overlap, breakpoint_tolerance_bp=tol)

    def test_exact_match(self) -> None:
        query = _make_row(start=100, end=200)
        cands = [_make_row(start=100, end=200, locus_id="GPX000001")]
        best, ov = find_best_match(query, cands, self._cfg())
        assert best is not None
        assert best.locus_id == "GPX000001"
        assert ov == pytest.approx(1.0)

    def test_no_overlap_returns_none(self) -> None:
        query = _make_row(start=100, end=200)
        cands = [_make_row(start=500, end=600, locus_id="GPX000001")]
        best, ov = find_best_match(query, cands, self._cfg())
        assert best is None
        assert ov == pytest.approx(0.0)

    def test_picks_best_overlap(self) -> None:
        query = _make_row(start=100, end=300)
        cand_a = _make_row(start=100, end=300, locus_id="GPX000001")  # perfect
        cand_b = _make_row(start=200, end=400, locus_id="GPX000002")  # partial
        best, ov = find_best_match(query, [cand_b, cand_a], self._cfg())
        assert best is not None
        assert best.locus_id == "GPX000001"

    def test_breakpoint_fallback(self) -> None:
        query = _make_row(start=100, end=200)
        near = _make_row(start=210, end=310, locus_id="GPX000001")  # gap=10 bp
        best, ov = find_best_match(query, [near], self._cfg(min_overlap=0.5, tol=50))
        assert best is not None
        assert best.locus_id == "GPX000001"

    def test_breakpoint_fallback_disabled(self) -> None:
        query = _make_row(start=100, end=200)
        near = _make_row(start=210, end=310, locus_id="GPX000001")  # gap=10 bp
        best, _ = find_best_match(query, [near], self._cfg(min_overlap=0.5, tol=0))
        assert best is None

    def test_empty_candidates(self) -> None:
        query = _make_row()
        best, ov = find_best_match(query, [], self._cfg())
        assert best is None
        assert ov == pytest.approx(0.0)

    def test_different_contig_no_match(self) -> None:
        query = _make_row(contig="chr1", start=100, end=200)
        cand = _make_row(contig="chr2", start=100, end=200, locus_id="GPX000001")
        best, _ = find_best_match(query, [cand], self._cfg())
        assert best is None
