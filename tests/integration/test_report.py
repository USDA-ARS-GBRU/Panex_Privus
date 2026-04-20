"""Integration tests for ``privy report``.

Tests the full ``run_report()`` pipeline using synthetic TSV inputs that
mirror real ``privy scan`` output.  No actual scan is executed here — inputs
are written directly so the tests are fast and self-contained.

Test classes:
    TestOutputFiles        — all expected output files are created
    TestSummaryTsv         — summary.tsv content is correct
    TestRankedHitsTsv      — ranked_hits.tsv ordering and top-N
    TestStrictnessSummary  — strictness_summary.tsv content
    TestMarkdownReport     — report.md structure and key sections
    TestHtmlReport         — report.html created with correct format
    TestWithOptionalInputs — regions.tsv, qc.tsv, evidence.tsv wired through
    TestCli                — privy report CLI command works end-to-end
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterator

import pytest

from privy.core.config import default_config
from privy.io.tsv import (
    EVIDENCE_COLUMNS,
    HITS_COLUMNS,
    QC_COLUMNS,
    REGIONS_COLUMNS,
    STRICTNESS_SUMMARY_COLUMNS,
    TsvWriter,
    read_tsv,
)
from privy.report.summary import run_report


# ---------------------------------------------------------------------------
# Synthetic input data
# ---------------------------------------------------------------------------

_HITS: list[dict[str, str]] = [
    {
        "locus_id": "PPX00000001", "contig": "chr1", "start": "99", "end": "100",
        "variant_type": "snp", "allele_key": "chr1:100:A:T",
        "target_support_n": "2", "target_total_n": "2",
        "offtarget_support_n": "0", "offtarget_total_n": "3",
        "target_missing_n": "0", "offtarget_missing_n": "0",
        "strictness_class": "strict_complete",
        "discovery_score": "1.0", "support_score": "0.0",
        "penalty_score": "0.0", "final_score": "1.0",
    },
    {
        "locus_id": "PPX00000002", "contig": "chr1", "start": "199", "end": "200",
        "variant_type": "snp", "allele_key": "chr1:200:A:T",
        "target_support_n": "1", "target_total_n": "2",
        "offtarget_support_n": "0", "offtarget_total_n": "3",
        "target_missing_n": "1", "offtarget_missing_n": "0",
        "strictness_class": "strict_target_missing",
        "discovery_score": "0.85", "support_score": "0.0",
        "penalty_score": "0.1", "final_score": "0.75",
    },
    {
        "locus_id": "PPX00000003", "contig": "chr1", "start": "299", "end": "300",
        "variant_type": "snp", "allele_key": "chr1:300:A:T",
        "target_support_n": "2", "target_total_n": "2",
        "offtarget_support_n": "0", "offtarget_total_n": "2",
        "target_missing_n": "0", "offtarget_missing_n": "1",
        "strictness_class": "strict_offtarget_missing",
        "discovery_score": "0.9", "support_score": "0.0",
        "penalty_score": "0.05", "final_score": "0.85",
    },
    {
        "locus_id": "PPX00000004", "contig": "chr1", "start": "399", "end": "400",
        "variant_type": "snp", "allele_key": "chr1:400:A:T",
        "target_support_n": "1", "target_total_n": "2",
        "offtarget_support_n": "0", "offtarget_total_n": "2",
        "target_missing_n": "1", "offtarget_missing_n": "1",
        "strictness_class": "strict_both_missing",
        "discovery_score": "0.7", "support_score": "0.0",
        "penalty_score": "0.15", "final_score": "0.55",
    },
]

_QC_ROWS: list[dict[str, str]] = [
    {"metric": "records_evaluated",    "value": "50", "description": "VCF records processed"},
    {"metric": "alleles_passed",        "value": "4",  "description": "alleles passing"},
    {"metric": "alleles_contradicted",  "value": "2",  "description": "alleles contradicted"},
    {"metric": "loci_emitted",          "value": "4",  "description": "loci written"},
    {"metric": "regions_emitted",       "value": "1",  "description": "regions written"},
]

_REGIONS: list[dict[str, str]] = [
    {
        "region_id": "REGION_00000001", "contig": "chr1",
        "start": "99", "end": "400",
        "n_loci": "4", "variant_types": "snp",
        "dominant_strictness_class": "strict_complete",
        "target_consistency": "0.875", "offtarget_exclusion": "1.0",
        "final_score": "0.79",
    },
]

_EVIDENCE: list[dict[str, str]] = [
    {
        "locus_id": "PPX00000001", "source_type": "vcf",
        "sample_id": "T1", "evidence_class": "support",
        "metric_name": "gt", "metric_value": "0/1", "details": "",
    },
    {
        "locus_id": "PPX00000001", "source_type": "vcf",
        "sample_id": "O1", "evidence_class": "absence",
        "metric_name": "gt", "metric_value": "0/0", "details": "",
    },
]


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def scan_dir(tmp_path: Path) -> Path:
    """Write synthetic scan outputs to a temp directory."""
    d = tmp_path / "scan"
    d.mkdir()

    with TsvWriter(d / "hits.tsv", HITS_COLUMNS) as w:
        w.write_rows(_HITS)
    with TsvWriter(d / "qc.tsv", QC_COLUMNS) as w:
        w.write_rows(_QC_ROWS)
    with TsvWriter(d / "regions.tsv", REGIONS_COLUMNS) as w:
        w.write_rows(_REGIONS)
    with TsvWriter(d / "evidence.tsv", EVIDENCE_COLUMNS) as w:
        w.write_rows(_EVIDENCE)

    return d


@pytest.fixture
def report_dir(tmp_path: Path) -> Path:
    d = tmp_path / "report"
    d.mkdir()
    return d


@pytest.fixture
def base_report(scan_dir: Path, report_dir: Path) -> Path:
    """Run report with all inputs and return the report output directory."""
    run_report(
        hits=scan_dir / "hits.tsv",
        regions=scan_dir / "regions.tsv",
        evidence=scan_dir / "evidence.tsv",
        compare=None,
        qc=scan_dir / "qc.tsv",
        run_json=None,
        cfg=default_config(),
        fmt="both",
        top_n=20,
        include_qc=True,
        include_strictness=True,
        include_compare=True,
        include_regions=True,
        title="Integration Test Report",
        outdir=report_dir,
    )
    return report_dir


# ---------------------------------------------------------------------------
# TestOutputFiles
# ---------------------------------------------------------------------------

class TestOutputFiles:
    def test_summary_tsv_created(self, base_report: Path) -> None:
        assert (base_report / "summary.tsv").exists()

    def test_ranked_hits_tsv_created(self, base_report: Path) -> None:
        assert (base_report / "ranked_hits.tsv").exists()

    def test_strictness_summary_tsv_created(self, base_report: Path) -> None:
        assert (base_report / "strictness_summary.tsv").exists()

    def test_report_md_created(self, base_report: Path) -> None:
        assert (base_report / "report.md").exists()

    def test_report_html_created(self, base_report: Path) -> None:
        assert (base_report / "report.html").exists()

    def test_support_summary_created_when_evidence_provided(self, base_report: Path) -> None:
        assert (base_report / "support_summary.tsv").exists()

    def test_contradiction_summary_created_when_qc_provided(self, base_report: Path) -> None:
        assert (base_report / "contradiction_summary.tsv").exists()

    def test_no_html_for_markdown_only(self, scan_dir: Path, tmp_path: Path) -> None:
        out = tmp_path / "md_only"
        out.mkdir()
        run_report(
            hits=scan_dir / "hits.tsv",
            regions=None, evidence=None, compare=None,
            qc=None, run_json=None,
            cfg=default_config(), fmt="markdown",
            top_n=20, include_qc=True, include_strictness=True,
            include_compare=True, include_regions=True,
            title="MD Only", outdir=out,
        )
        assert (out / "report.md").exists()
        assert not (out / "report.html").exists()


# ---------------------------------------------------------------------------
# TestSummaryTsv
# ---------------------------------------------------------------------------

class TestSummaryTsv:
    def test_has_correct_columns(self, base_report: Path) -> None:
        rows = read_tsv(base_report / "summary.tsv")
        assert list(rows[0].keys()) == QC_COLUMNS

    def test_n_hits_is_correct(self, base_report: Path) -> None:
        rows = {r["metric"]: r["value"] for r in read_tsv(base_report / "summary.tsv")}
        assert rows["n_hits"] == str(len(_HITS))

    def test_n_regions_is_correct(self, base_report: Path) -> None:
        rows = {r["metric"]: r["value"] for r in read_tsv(base_report / "summary.tsv")}
        assert rows["n_regions"] == str(len(_REGIONS))

    def test_top_locus_id_is_highest_score(self, base_report: Path) -> None:
        rows = {r["metric"]: r["value"] for r in read_tsv(base_report / "summary.tsv")}
        assert rows["top_locus_id"] == "PPX00000001"
        assert rows["top_final_score"] == "1.0"

    def test_records_evaluated_from_qc(self, base_report: Path) -> None:
        rows = {r["metric"]: r["value"] for r in read_tsv(base_report / "summary.tsv")}
        assert rows.get("records_evaluated") == "50"


# ---------------------------------------------------------------------------
# TestRankedHitsTsv
# ---------------------------------------------------------------------------

class TestRankedHitsTsv:
    def test_rank_column_present(self, base_report: Path) -> None:
        rows = read_tsv(base_report / "ranked_hits.tsv")
        assert "rank" in rows[0]

    def test_rank_1_has_highest_score(self, base_report: Path) -> None:
        rows = read_tsv(base_report / "ranked_hits.tsv")
        assert rows[0]["rank"] == "1"
        assert rows[0]["final_score"] == "1.0"
        assert rows[0]["locus_id"] == "PPX00000001"

    def test_rows_sorted_descending_by_score(self, base_report: Path) -> None:
        rows = read_tsv(base_report / "ranked_hits.tsv")
        scores = [float(r["final_score"]) for r in rows]
        assert scores == sorted(scores, reverse=True)

    def test_all_hits_columns_present(self, base_report: Path) -> None:
        rows = read_tsv(base_report / "ranked_hits.tsv")
        for col in HITS_COLUMNS:
            assert col in rows[0]

    def test_top_n_limits_rows(self, scan_dir: Path, tmp_path: Path) -> None:
        out = tmp_path / "top2"
        out.mkdir()
        run_report(
            hits=scan_dir / "hits.tsv",
            regions=None, evidence=None, compare=None,
            qc=None, run_json=None,
            cfg=default_config(), fmt="markdown",
            top_n=2, include_qc=True, include_strictness=True,
            include_compare=True, include_regions=True,
            title="Top2", outdir=out,
        )
        rows = read_tsv(out / "ranked_hits.tsv")
        assert len(rows) == 2


# ---------------------------------------------------------------------------
# TestStrictnessSummary
# ---------------------------------------------------------------------------

class TestStrictnessSummary:
    def test_columns_match_schema(self, base_report: Path) -> None:
        rows = read_tsv(base_report / "strictness_summary.tsv")
        assert list(rows[0].keys()) == STRICTNESS_SUMMARY_COLUMNS

    def test_strict_complete_count(self, base_report: Path) -> None:
        rows = read_tsv(base_report / "strictness_summary.tsv")
        by_class = {r["strictness_class"]: r["n_loci"] for r in rows}
        assert by_class["strict_complete"] == "1"

    def test_all_canonical_classes_listed(self, base_report: Path) -> None:
        rows = read_tsv(base_report / "strictness_summary.tsv")
        classes = {r["strictness_class"] for r in rows}
        for expected in (
            "strict_complete", "strict_target_missing",
            "strict_offtarget_missing", "strict_both_missing",
        ):
            assert expected in classes

    def test_not_written_when_include_strictness_false(
        self, scan_dir: Path, tmp_path: Path
    ) -> None:
        out = tmp_path / "no_strict"
        out.mkdir()
        run_report(
            hits=scan_dir / "hits.tsv",
            regions=None, evidence=None, compare=None,
            qc=None, run_json=None,
            cfg=default_config(), fmt="markdown",
            top_n=20, include_qc=True, include_strictness=False,
            include_compare=True, include_regions=True,
            title="No Strictness", outdir=out,
        )
        assert not (out / "strictness_summary.tsv").exists()


# ---------------------------------------------------------------------------
# TestMarkdownReport
# ---------------------------------------------------------------------------

class TestMarkdownReport:
    def test_starts_with_h1_title(self, base_report: Path) -> None:
        text = (base_report / "report.md").read_text()
        assert text.startswith("# Integration Test Report")

    def test_contains_run_summary_section(self, base_report: Path) -> None:
        text = (base_report / "report.md").read_text()
        assert "## Run Summary" in text

    def test_contains_top_hits_section(self, base_report: Path) -> None:
        text = (base_report / "report.md").read_text()
        assert "## Top" in text
        assert "Hits" in text

    def test_contains_strictness_section(self, base_report: Path) -> None:
        text = (base_report / "report.md").read_text()
        assert "## Strictness Class Distribution" in text

    def test_contains_qc_section_when_enabled(self, base_report: Path) -> None:
        text = (base_report / "report.md").read_text()
        assert "## Filtering and QC" in text

    def test_contains_regions_section(self, base_report: Path) -> None:
        text = (base_report / "report.md").read_text()
        assert "## Candidate Regions" in text

    def test_contains_caveats_section(self, base_report: Path) -> None:
        text = (base_report / "report.md").read_text()
        assert "## Caveats" in text

    def test_top_locus_id_appears_in_hits_table(self, base_report: Path) -> None:
        text = (base_report / "report.md").read_text()
        assert "PPX00000001" in text

    def test_qc_omitted_when_include_qc_false(
        self, scan_dir: Path, tmp_path: Path
    ) -> None:
        out = tmp_path / "no_qc"
        out.mkdir()
        run_report(
            hits=scan_dir / "hits.tsv",
            regions=None, evidence=None, compare=None,
            qc=scan_dir / "qc.tsv", run_json=None,
            cfg=default_config(), fmt="markdown",
            top_n=20, include_qc=False, include_strictness=True,
            include_compare=True, include_regions=True,
            title="No QC", outdir=out,
        )
        text = (out / "report.md").read_text()
        assert "## Filtering and QC" not in text


# ---------------------------------------------------------------------------
# TestHtmlReport
# ---------------------------------------------------------------------------

class TestHtmlReport:
    def test_html_is_valid_document(self, base_report: Path) -> None:
        text = (base_report / "report.html").read_text()
        assert "<!DOCTYPE html>" in text
        assert "<html" in text
        assert "</html>" in text

    def test_html_contains_title(self, base_report: Path) -> None:
        text = (base_report / "report.html").read_text()
        assert "Integration Test Report" in text

    def test_html_has_table_elements(self, base_report: Path) -> None:
        text = (base_report / "report.html").read_text()
        assert "<table>" in text

    def test_html_contains_locus_id(self, base_report: Path) -> None:
        text = (base_report / "report.html").read_text()
        assert "PPX00000001" in text


# ---------------------------------------------------------------------------
# TestWithOptionalInputs
# ---------------------------------------------------------------------------

class TestWithOptionalInputs:
    def test_no_optional_inputs_succeeds(self, scan_dir: Path, tmp_path: Path) -> None:
        out = tmp_path / "minimal"
        out.mkdir()
        run_report(
            hits=scan_dir / "hits.tsv",
            regions=None, evidence=None, compare=None,
            qc=None, run_json=None,
            cfg=default_config(), fmt="markdown",
            top_n=20, include_qc=True, include_strictness=True,
            include_compare=True, include_regions=True,
            title="Minimal", outdir=out,
        )
        assert (out / "report.md").exists()

    def test_support_summary_not_written_without_evidence(
        self, scan_dir: Path, tmp_path: Path
    ) -> None:
        out = tmp_path / "no_evid"
        out.mkdir()
        run_report(
            hits=scan_dir / "hits.tsv",
            regions=None, evidence=None, compare=None,
            qc=None, run_json=None,
            cfg=default_config(), fmt="markdown",
            top_n=20, include_qc=True, include_strictness=True,
            include_compare=True, include_regions=True,
            title="No Evidence", outdir=out,
        )
        assert not (out / "support_summary.tsv").exists()

    def test_regions_section_in_report_when_provided(self, base_report: Path) -> None:
        text = (base_report / "report.md").read_text()
        assert "REGION_00000001" in text

    def test_contradiction_summary_has_qc_metric(self, base_report: Path) -> None:
        rows = read_tsv(base_report / "contradiction_summary.tsv")
        metrics = {r["metric"] for r in rows}
        assert "alleles_contradicted" in metrics

    def test_support_summary_columns_correct(self, base_report: Path) -> None:
        from privy.io.tsv import SUPPORT_SUMMARY_COLUMNS
        rows = read_tsv(base_report / "support_summary.tsv")
        assert list(rows[0].keys()) == SUPPORT_SUMMARY_COLUMNS


# ---------------------------------------------------------------------------
# TestCli
# ---------------------------------------------------------------------------

class TestCli:
    def test_privy_report_runs_via_cli(self, scan_dir: Path, tmp_path: Path) -> None:
        from typer.testing import CliRunner
        from privy.cli.main import app

        out = tmp_path / "cli_report"
        out.mkdir()

        runner = CliRunner()
        result = runner.invoke(app, [
            "report",
            "--hits", str(scan_dir / "hits.tsv"),
            "--qc", str(scan_dir / "qc.tsv"),
            "--outdir", str(out),
            "--title", "CLI Test",
        ])

        assert result.exit_code == 0, result.output
        assert (out / "report.md").exists()

    def test_privy_report_fails_without_hits(self, tmp_path: Path) -> None:
        from typer.testing import CliRunner
        from privy.cli.main import app

        runner = CliRunner()
        result = runner.invoke(app, ["report", "--outdir", str(tmp_path)])
        assert result.exit_code != 0

    def test_privy_report_fails_with_missing_hits_file(self, tmp_path: Path) -> None:
        from typer.testing import CliRunner
        from privy.cli.main import app

        runner = CliRunner()
        result = runner.invoke(app, [
            "report",
            "--hits", str(tmp_path / "nonexistent_hits.tsv"),
            "--outdir", str(tmp_path),
        ])
        assert result.exit_code != 0
