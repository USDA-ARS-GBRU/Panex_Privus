"""Unit tests for privy.report.summary helper functions.

Tests cover:
    - _rank_hits: sort order and top-N truncation
    - _compute_strictness_summary: counts, percentages, canonical ordering
    - _compute_support_summary: grouping by source_type × evidence_class
    - _compute_contradiction_summary: extraction from qc_rows and compare_rows
    - _compute_run_summary: key metrics assembled correctly
    - _write_* helpers: correct column schemas and row contents
"""

from __future__ import annotations

from pathlib import Path

from privy.core.config import default_config
from privy.io.tsv import (
    HITS_COLUMNS,
    QC_COLUMNS,
    RANKED_HITS_COLUMNS,
    STRICTNESS_SUMMARY_COLUMNS,
    read_tsv,
)
from privy.report.summary import (
    _compute_contradiction_summary,
    _compute_run_summary,
    _compute_strictness_summary,
    _compute_support_summary,
    _rank_hits,
    _write_ranked_hits_tsv,
    _write_strictness_summary_tsv,
    _write_summary_tsv,
    run_report,
)

# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

def _make_hit(
    locus_id: str = "PPX00000001",
    contig: str = "chr1",
    start: str = "99",
    end: str = "100",
    variant_type: str = "snp",
    allele_key: str = "chr1:100:A:T",
    target_support_n: str = "2",
    target_total_n: str = "2",
    offtarget_support_n: str = "0",
    offtarget_total_n: str = "3",
    target_missing_n: str = "0",
    offtarget_missing_n: str = "0",
    strictness_class: str = "strict_complete",
    discovery_score: str = "1.0",
    support_score: str = "0.0",
    penalty_score: str = "0.0",
    final_score: str = "1.0",
) -> dict[str, str]:
    return {
        "locus_id": locus_id,
        "contig": contig,
        "start": start,
        "end": end,
        "variant_type": variant_type,
        "allele_key": allele_key,
        "target_support_n": target_support_n,
        "target_total_n": target_total_n,
        "offtarget_support_n": offtarget_support_n,
        "offtarget_total_n": offtarget_total_n,
        "target_missing_n": target_missing_n,
        "offtarget_missing_n": offtarget_missing_n,
        "strictness_class": strictness_class,
        "discovery_score": discovery_score,
        "support_score": support_score,
        "penalty_score": penalty_score,
        "final_score": final_score,
    }


SAMPLE_HITS = [
    _make_hit("PPX00000001", final_score="0.85", strictness_class="strict_complete"),
    _make_hit("PPX00000002", final_score="1.0",  strictness_class="strict_complete"),
    _make_hit("PPX00000003", final_score="0.55", strictness_class="strict_target_missing"),
    _make_hit("PPX00000004", final_score="0.75", strictness_class="strict_offtarget_missing"),
]


# ---------------------------------------------------------------------------
# _rank_hits
# ---------------------------------------------------------------------------

class TestRankHits:
    def test_sorts_descending_by_final_score(self) -> None:
        ranked = _rank_hits(SAMPLE_HITS, top_n=10)
        scores = [float(r["final_score"]) for r in ranked]
        assert scores == sorted(scores, reverse=True)

    def test_top_n_truncation(self) -> None:
        ranked = _rank_hits(SAMPLE_HITS, top_n=2)
        assert len(ranked) == 2
        assert ranked[0]["locus_id"] == "PPX00000002"  # score 1.0

    def test_returns_all_when_top_n_exceeds_count(self) -> None:
        ranked = _rank_hits(SAMPLE_HITS, top_n=100)
        assert len(ranked) == len(SAMPLE_HITS)

    def test_empty_input(self) -> None:
        assert _rank_hits([], top_n=10) == []

    def test_handles_missing_final_score(self) -> None:
        rows = [{"locus_id": "X", "final_score": ""}, {"locus_id": "Y", "final_score": "0.5"}]
        ranked = _rank_hits(rows, top_n=10)
        assert ranked[0]["locus_id"] == "Y"

    def test_handles_non_numeric_score(self) -> None:
        rows = [{"locus_id": "X", "final_score": "NA"}, {"locus_id": "Y", "final_score": "0.3"}]
        ranked = _rank_hits(rows, top_n=10)
        assert ranked[0]["locus_id"] == "Y"


# ---------------------------------------------------------------------------
# _compute_strictness_summary
# ---------------------------------------------------------------------------

class TestComputeStrictnessSummary:
    def test_counts_classes_correctly(self) -> None:
        rows = _compute_strictness_summary(SAMPLE_HITS)
        by_class = {r["strictness_class"]: r["n_loci"] for r in rows}
        assert by_class["strict_complete"] == 2
        assert by_class["strict_target_missing"] == 1
        assert by_class["strict_offtarget_missing"] == 1
        assert by_class["strict_both_missing"] == 0

    def test_percentages_sum_to_100(self) -> None:
        rows = _compute_strictness_summary(SAMPLE_HITS)
        total_pct = sum(r["pct_hits"] for r in rows if isinstance(r["pct_hits"], float))
        assert abs(total_pct - 100.0) < 0.2  # allow rounding error

    def test_canonical_order(self) -> None:
        rows = _compute_strictness_summary(SAMPLE_HITS)
        classes = [r["strictness_class"] for r in rows]
        assert classes[0] == "strict_complete"
        assert classes[1] == "strict_target_missing"

    def test_empty_input(self) -> None:
        rows = _compute_strictness_summary([])
        assert all(r["n_loci"] == 0 for r in rows)
        assert all(r["pct_hits"] == 0.0 for r in rows)

    def test_zero_pct_for_absent_class(self) -> None:
        rows = _compute_strictness_summary(
            [_make_hit(strictness_class="strict_complete")]
        )
        by_class = {r["strictness_class"]: r["pct_hits"] for r in rows}
        assert by_class["strict_target_missing"] == 0.0


# ---------------------------------------------------------------------------
# _compute_support_summary
# ---------------------------------------------------------------------------

class TestComputeSupportSummary:
    def test_groups_by_source_and_class(self) -> None:
        rows = [
            {"source_type": "vcf", "evidence_class": "support"},
            {"source_type": "vcf", "evidence_class": "support"},
            {"source_type": "vcf", "evidence_class": "absence"},
            {"source_type": "bam", "evidence_class": "support"},
        ]
        result = _compute_support_summary(rows)
        by_key = {(r["source_type"], r["evidence_class"]): r["n_records"] for r in result}
        assert by_key[("vcf", "support")] == 2
        assert by_key[("vcf", "absence")] == 1
        assert by_key[("bam", "support")] == 1

    def test_pct_of_source_sums_to_100_per_source(self) -> None:
        rows = [
            {"source_type": "vcf", "evidence_class": "support"},
            {"source_type": "vcf", "evidence_class": "absence"},
        ]
        result = _compute_support_summary(rows)
        vcf_rows = [r for r in result if r["source_type"] == "vcf"]
        total_pct = sum(r["pct_of_source"] for r in vcf_rows)
        assert abs(total_pct - 100.0) < 0.2

    def test_empty_input(self) -> None:
        assert _compute_support_summary([]) == []


# ---------------------------------------------------------------------------
# _compute_contradiction_summary
# ---------------------------------------------------------------------------

class TestComputeContradictionSummary:
    def test_extracts_from_qc_rows(self) -> None:
        qc = [
            {"metric": "alleles_contradicted", "value": "5", "description": "..."},
            {"metric": "n_hits", "value": "20", "description": "..."},
        ]
        result = _compute_contradiction_summary(qc, None)
        by_metric = {r["metric"]: r["value"] for r in result}
        assert by_metric["alleles_contradicted"] == "5"

    def test_counts_compare_contradicted_loci(self) -> None:
        compare = [
            {"match_class": "contradicted"},
            {"match_class": "supported"},
            {"match_class": "contradicted"},
        ]
        result = _compute_contradiction_summary(None, compare)
        by_metric = {r["metric"]: r["value"] for r in result}
        assert by_metric["compare_contradicted_loci"] == "2"

    def test_none_inputs(self) -> None:
        assert _compute_contradiction_summary(None, None) == []


# ---------------------------------------------------------------------------
# _compute_run_summary
# ---------------------------------------------------------------------------

class TestComputeRunSummary:
    def test_n_hits_matches_input(self) -> None:
        cfg = default_config()
        rows = _compute_run_summary(SAMPLE_HITS, None, None, None, cfg)
        by_metric = {r["metric"]: r["value"] for r in rows}
        assert by_metric["n_hits"] == str(len(SAMPLE_HITS))

    def test_top_locus_is_highest_score(self) -> None:
        cfg = default_config()
        rows = _compute_run_summary(SAMPLE_HITS, None, None, None, cfg)
        by_metric = {r["metric"]: r["value"] for r in rows}
        assert by_metric.get("top_final_score") == "1.0"
        assert by_metric.get("top_locus_id") == "PPX00000002"

    def test_project_name_from_config(self) -> None:
        cfg = default_config()
        cfg.project_name = "my_project"
        rows = _compute_run_summary([], None, None, None, cfg)
        by_metric = {r["metric"]: r["value"] for r in rows}
        assert by_metric["project_name"] == "my_project"

    def test_n_regions_included_when_provided(self) -> None:
        cfg = default_config()
        regions = [{"region_id": "R1"}, {"region_id": "R2"}]
        rows = _compute_run_summary(SAMPLE_HITS, regions, None, None, cfg)  # type: ignore[arg-type]
        by_metric = {r["metric"]: r["value"] for r in rows}
        assert by_metric.get("n_regions") == "2"

    def test_empty_hits(self) -> None:
        cfg = default_config()
        rows = _compute_run_summary([], None, None, None, cfg)
        by_metric = {r["metric"]: r["value"] for r in rows}
        assert by_metric["n_hits"] == "0"
        assert "top_locus_id" not in by_metric


# ---------------------------------------------------------------------------
# TSV writer helpers
# ---------------------------------------------------------------------------

class TestWriteSummaryTsv:
    def test_creates_file_with_correct_columns(self, tmp_path: Path) -> None:
        rows = [{"metric": "n_hits", "value": "5", "description": "hits"}]
        _write_summary_tsv(rows, tmp_path)
        result = read_tsv(tmp_path / "summary.tsv")
        assert result[0]["metric"] == "n_hits"
        assert list(result[0].keys()) == QC_COLUMNS

    def test_empty_rows_writes_header_only(self, tmp_path: Path) -> None:
        _write_summary_tsv([], tmp_path)
        content = (tmp_path / "summary.tsv").read_text()
        assert content.startswith("\t".join(QC_COLUMNS))


class TestWriteRankedHitsTsv:
    def test_creates_file_with_rank_column(self, tmp_path: Path) -> None:
        _write_ranked_hits_tsv(SAMPLE_HITS, tmp_path)
        result = read_tsv(tmp_path / "ranked_hits.tsv")
        assert "rank" in result[0]
        assert result[0]["rank"] == "1"

    def test_rank_increments_correctly(self, tmp_path: Path) -> None:
        _write_ranked_hits_tsv(SAMPLE_HITS[:3], tmp_path)
        result = read_tsv(tmp_path / "ranked_hits.tsv")
        assert [r["rank"] for r in result] == ["1", "2", "3"]

    def test_all_hits_columns_present(self, tmp_path: Path) -> None:
        _write_ranked_hits_tsv(SAMPLE_HITS, tmp_path)
        result = read_tsv(tmp_path / "ranked_hits.tsv")
        for col in HITS_COLUMNS:
            assert col in result[0], f"missing column: {col}"

    def test_column_order_matches_schema(self, tmp_path: Path) -> None:
        _write_ranked_hits_tsv(SAMPLE_HITS, tmp_path)
        header = (tmp_path / "ranked_hits.tsv").read_text().splitlines()[0]
        assert header == "\t".join(RANKED_HITS_COLUMNS)


class TestWriteStrictnessSummaryTsv:
    def test_columns_match_schema(self, tmp_path: Path) -> None:
        rows = _compute_strictness_summary(SAMPLE_HITS)
        _write_strictness_summary_tsv(rows, tmp_path)
        result = read_tsv(tmp_path / "strictness_summary.tsv")
        assert list(result[0].keys()) == STRICTNESS_SUMMARY_COLUMNS

    def test_all_canonical_classes_present(self, tmp_path: Path) -> None:
        rows = _compute_strictness_summary(SAMPLE_HITS)
        _write_strictness_summary_tsv(rows, tmp_path)
        result = read_tsv(tmp_path / "strictness_summary.tsv")
        classes = {r["strictness_class"] for r in result}
        assert "strict_complete" in classes
        assert "strict_target_missing" in classes


# ---------------------------------------------------------------------------
# run_report — smoke test with minimal inputs
# ---------------------------------------------------------------------------

class TestRunReportSmoke:
    def _write_hits(self, path: Path, rows: list[dict[str, str]]) -> None:
        from privy.io.tsv import TsvWriter
        with TsvWriter(path, HITS_COLUMNS) as w:
            w.write_rows(rows)

    def test_creates_all_output_files_markdown(self, tmp_path: Path) -> None:
        hits_path = tmp_path / "hits.tsv"
        self._write_hits(hits_path, SAMPLE_HITS)
        outdir = tmp_path / "report"
        outdir.mkdir()

        run_report(
            hits=hits_path,
            regions=None,
            evidence=None,
            compare=None,
            qc=None,
            run_json=None,
            cfg=default_config(),
            fmt="markdown",
            top_n=20,
            include_qc=True,
            include_strictness=True,
            include_compare=True,
            include_regions=True,
            title="Test Report",
            outdir=outdir,
        )

        assert (outdir / "summary.tsv").exists()
        assert (outdir / "ranked_hits.tsv").exists()
        assert (outdir / "strictness_summary.tsv").exists()
        assert (outdir / "report.md").exists()
        assert not (outdir / "report.html").exists()

    def test_html_format_creates_html(self, tmp_path: Path) -> None:
        hits_path = tmp_path / "hits.tsv"
        self._write_hits(hits_path, SAMPLE_HITS)
        outdir = tmp_path / "report"
        outdir.mkdir()

        run_report(
            hits=hits_path,
            regions=None,
            evidence=None,
            compare=None,
            qc=None,
            run_json=None,
            cfg=default_config(),
            fmt="html",
            top_n=20,
            include_qc=True,
            include_strictness=True,
            include_compare=True,
            include_regions=True,
            title="Test Report",
            outdir=outdir,
        )

        assert (outdir / "report.md").exists()
        assert (outdir / "report.html").exists()

    def test_both_format_creates_both(self, tmp_path: Path) -> None:
        hits_path = tmp_path / "hits.tsv"
        self._write_hits(hits_path, SAMPLE_HITS)
        outdir = tmp_path / "report"
        outdir.mkdir()

        run_report(
            hits=hits_path,
            regions=None,
            evidence=None,
            compare=None,
            qc=None,
            run_json=None,
            cfg=default_config(),
            fmt="both",
            top_n=20,
            include_qc=True,
            include_strictness=True,
            include_compare=True,
            include_regions=True,
            title="Test Report",
            outdir=outdir,
        )

        assert (outdir / "report.md").exists()
        assert (outdir / "report.html").exists()

    def test_empty_hits_does_not_raise(self, tmp_path: Path) -> None:
        hits_path = tmp_path / "hits.tsv"
        self._write_hits(hits_path, [])
        outdir = tmp_path / "report"
        outdir.mkdir()

        run_report(
            hits=hits_path,
            regions=None,
            evidence=None,
            compare=None,
            qc=None,
            run_json=None,
            cfg=default_config(),
            fmt="markdown",
            top_n=20,
            include_qc=True,
            include_strictness=True,
            include_compare=True,
            include_regions=True,
            title="Empty Report",
            outdir=outdir,
        )

        assert (outdir / "report.md").exists()

    def test_top_n_truncates_ranked_hits(self, tmp_path: Path) -> None:
        hits_path = tmp_path / "hits.tsv"
        self._write_hits(hits_path, SAMPLE_HITS)
        outdir = tmp_path / "report"
        outdir.mkdir()

        run_report(
            hits=hits_path,
            regions=None,
            evidence=None,
            compare=None,
            qc=None,
            run_json=None,
            cfg=default_config(),
            fmt="markdown",
            top_n=2,
            include_qc=True,
            include_strictness=True,
            include_compare=True,
            include_regions=True,
            title="Truncated Report",
            outdir=outdir,
        )

        ranked = read_tsv(outdir / "ranked_hits.tsv")
        assert len(ranked) == 2
        assert ranked[0]["final_score"] == "1.0"
