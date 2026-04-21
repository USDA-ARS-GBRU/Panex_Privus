"""Integration tests for privy compare — run_compare() and CLI."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Generator

import pytest
from typer.testing import CliRunner

from privy.backends.compare import run_compare
from privy.cli.main import app
from privy.core.config import default_config
from privy.core.evidence import MatchClass
from privy.io.tsv import COMPARE_COLUMNS, COMPARE_SUMMARY_COLUMNS, read_tsv


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

_VCF_HITS_HEADER = (
    "locus_id\tcontig\tstart\tend\tvariant_type\tallele_key\t"
    "target_support_n\ttarget_total_n\tofftarget_support_n\t"
    "offtarget_total_n\ttarget_missing_n\tofftarget_missing_n\t"
    "strictness_class\tdiscovery_score\tsupport_score\tpenalty_score\tfinal_score\n"
)


def _write_hits(path: Path, rows: list[dict[str, str]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        fh.write(_VCF_HITS_HEADER)
        writer = csv.DictWriter(
            fh,
            fieldnames=_VCF_HITS_HEADER.strip().split("\t"),
            delimiter="\t",
            extrasaction="ignore",
        )
        for row in rows:
            writer.writerow(row)
    return path


def _hit(
    locus_id: str,
    contig: str,
    start: int,
    end: int,
    strictness: str = "strict_complete",
    final_score: float = 1.0,
) -> dict[str, str]:
    return {
        "locus_id": locus_id,
        "contig": contig,
        "start": str(start),
        "end": str(end),
        "variant_type": "snp",
        "allele_key": f"{contig}:{start}:A:T",
        "target_support_n": "3",
        "target_total_n": "3",
        "offtarget_support_n": "0",
        "offtarget_total_n": "2",
        "target_missing_n": "0",
        "offtarget_missing_n": "0",
        "strictness_class": strictness,
        "discovery_score": "1.0",
        "support_score": "0.0",
        "penalty_score": "0.0",
        "final_score": str(final_score),
    }


@pytest.fixture()
def vcf_hits(tmp_path: Path) -> Path:
    """VCF scan results with three loci."""
    return _write_hits(tmp_path / "vcf" / "hits.tsv", [
        _hit("PPX000001", "chr1", 100, 200, "strict_complete"),
        _hit("PPX000002", "chr1", 500, 600, "strict_target_missing", 0.7),
        _hit("PPX000003", "chr2", 1000, 1100, "strict_complete"),
    ])


@pytest.fixture()
def gfa_hits(tmp_path: Path) -> Path:
    """GFA scan results: matches PPX000001 and PPX000003; PPX000002 not present."""
    return _write_hits(tmp_path / "gfa" / "hits.tsv", [
        _hit("GPX000001", "chr1", 100, 200, "strict_complete"),      # matches PPX000001
        _hit("GPX000002", "chr2", 1000, 1100, "strict_complete"),    # matches PPX000003
        _hit("GPX000003", "chr3", 2000, 2100, "strict_complete"),    # no VCF counterpart
    ])


@pytest.fixture()
def runner() -> CliRunner:
    return CliRunner()


# ---------------------------------------------------------------------------
# TestRunCompareOutputFiles
# ---------------------------------------------------------------------------

class TestRunCompareOutputFiles:
    def test_compare_tsv_created(
        self, tmp_path: Path, vcf_hits: Path, gfa_hits: Path
    ) -> None:
        outdir = tmp_path / "out"
        run_compare(vcf_hits, gfa_hits, outdir, default_config())
        assert (outdir / "compare.tsv").exists()

    def test_compare_summary_tsv_created(
        self, tmp_path: Path, vcf_hits: Path, gfa_hits: Path
    ) -> None:
        outdir = tmp_path / "out"
        run_compare(vcf_hits, gfa_hits, outdir, default_config())
        assert (outdir / "compare_summary.tsv").exists()

    def test_compare_json_created(
        self, tmp_path: Path, vcf_hits: Path, gfa_hits: Path
    ) -> None:
        outdir = tmp_path / "out"
        run_compare(vcf_hits, gfa_hits, outdir, default_config())
        assert (outdir / "compare.json").exists()

    def test_write_flags_suppress_files(
        self, tmp_path: Path, vcf_hits: Path, gfa_hits: Path
    ) -> None:
        outdir = tmp_path / "out"
        run_compare(
            vcf_hits, gfa_hits, outdir, default_config(),
            write_compare_tsv=False,
            write_summary_tsv=False,
            write_json=False,
        )
        assert not (outdir / "compare.tsv").exists()
        assert not (outdir / "compare_summary.tsv").exists()
        assert not (outdir / "compare.json").exists()


# ---------------------------------------------------------------------------
# TestRunCompareColumns
# ---------------------------------------------------------------------------

class TestRunCompareColumns:
    def test_compare_tsv_has_correct_columns(
        self, tmp_path: Path, vcf_hits: Path, gfa_hits: Path
    ) -> None:
        outdir = tmp_path / "out"
        run_compare(vcf_hits, gfa_hits, outdir, default_config())
        rows = read_tsv(outdir / "compare.tsv")
        assert list(rows[0].keys()) == COMPARE_COLUMNS

    def test_compare_summary_tsv_has_correct_columns(
        self, tmp_path: Path, vcf_hits: Path, gfa_hits: Path
    ) -> None:
        outdir = tmp_path / "out"
        run_compare(vcf_hits, gfa_hits, outdir, default_config())
        rows = read_tsv(outdir / "compare_summary.tsv")
        assert list(rows[0].keys()) == COMPARE_SUMMARY_COLUMNS


# ---------------------------------------------------------------------------
# TestRunCompareMatchClasses
# ---------------------------------------------------------------------------

class TestRunCompareMatchClasses:
    def test_perfect_overlap_is_supported(
        self, tmp_path: Path, vcf_hits: Path, gfa_hits: Path
    ) -> None:
        outdir = tmp_path / "out"
        rows = run_compare(vcf_hits, gfa_hits, outdir, default_config())
        supported = [r for r in rows if r["match_class"] == MatchClass.SUPPORTED.value]
        assert len(supported) >= 2  # PPX000001 and PPX000003 both overlap GFA loci

    def test_vcf_only_locus_is_source_specific(
        self, tmp_path: Path, vcf_hits: Path, gfa_hits: Path
    ) -> None:
        outdir = tmp_path / "out"
        rows = run_compare(vcf_hits, gfa_hits, outdir, default_config())
        specific = [r for r in rows if r["match_class"] == MatchClass.SOURCE_SPECIFIC.value]
        # PPX000002 (chr1:500-600) has no GFA counterpart; GPX000003 (chr3) has no VCF counterpart
        assert len(specific) >= 2

    def test_vcf_only_locus_has_na_locus_id_b(
        self, tmp_path: Path, vcf_hits: Path, gfa_hits: Path
    ) -> None:
        outdir = tmp_path / "out"
        rows = run_compare(vcf_hits, gfa_hits, outdir, default_config())
        vcf_specific = [
            r for r in rows
            if r["match_class"] == MatchClass.SOURCE_SPECIFIC.value
            and r["locus_id_a"] != "NA"
        ]
        assert all(r["locus_id_b"] == "NA" for r in vcf_specific)

    def test_gfa_only_locus_has_na_locus_id_a(
        self, tmp_path: Path, vcf_hits: Path, gfa_hits: Path
    ) -> None:
        outdir = tmp_path / "out"
        rows = run_compare(vcf_hits, gfa_hits, outdir, default_config())
        gfa_specific = [
            r for r in rows
            if r["match_class"] == MatchClass.SOURCE_SPECIFIC.value
            and r["locus_id_b"] != "NA"
        ]
        assert all(r["locus_id_a"] == "NA" for r in gfa_specific)

    def test_total_row_count(
        self, tmp_path: Path, vcf_hits: Path, gfa_hits: Path
    ) -> None:
        # 3 rows in A (2 matched + 1 source-specific A) + 1 unmatched B = 4 total
        outdir = tmp_path / "out"
        rows = run_compare(vcf_hits, gfa_hits, outdir, default_config())
        assert len(rows) == 4

    def test_contradicted_classification(self, tmp_path: Path) -> None:
        vcf = _write_hits(tmp_path / "v" / "hits.tsv", [
            _hit("PPX000001", "chr1", 100, 200, "contradicted"),
        ])
        gfa = _write_hits(tmp_path / "g" / "hits.tsv", [
            _hit("GPX000001", "chr1", 100, 200, "strict_complete"),
        ])
        rows = run_compare(vcf, gfa, tmp_path / "out", default_config())
        assert rows[0]["match_class"] == MatchClass.CONTRADICTED.value

    def test_partially_supported_mixed_strictness(self, tmp_path: Path) -> None:
        cfg = default_config().model_copy(
            update={"compare": default_config().compare.model_copy(
                update={"require_state_compatibility": True}
            )}
        )
        vcf = _write_hits(tmp_path / "v" / "hits.tsv", [
            _hit("PPX000001", "chr1", 100, 200, "strict_complete"),
        ])
        gfa = _write_hits(tmp_path / "g" / "hits.tsv", [
            _hit("GPX000001", "chr1", 100, 200, "relaxed_threshold"),
        ])
        rows = run_compare(vcf, gfa, tmp_path / "out", cfg)
        assert rows[0]["match_class"] == MatchClass.PARTIALLY_SUPPORTED.value


# ---------------------------------------------------------------------------
# TestRunCompareSourceLabels
# ---------------------------------------------------------------------------

class TestRunCompareSourceLabels:
    def test_infers_vcf_label(
        self, tmp_path: Path, vcf_hits: Path, gfa_hits: Path
    ) -> None:
        outdir = tmp_path / "out"
        rows = run_compare(vcf_hits, gfa_hits, outdir, default_config())
        matched = [r for r in rows if r["source_a"] != "NA"]
        assert all(r["source_a"] == "vcf" for r in matched)

    def test_infers_gfa_label(
        self, tmp_path: Path, vcf_hits: Path, gfa_hits: Path
    ) -> None:
        outdir = tmp_path / "out"
        rows = run_compare(vcf_hits, gfa_hits, outdir, default_config())
        matched = [r for r in rows if r["source_b"] != "NA"]
        assert all(r["source_b"] == "gfa" for r in matched)

    def test_explicit_labels_used(
        self, tmp_path: Path, vcf_hits: Path, gfa_hits: Path
    ) -> None:
        outdir = tmp_path / "out"
        rows = run_compare(
            vcf_hits, gfa_hits, outdir, default_config(),
            source_label_a="my_vcf", source_label_b="my_gfa",
        )
        a_labels = {r["source_a"] for r in rows if r["source_a"] != "NA"}
        b_labels = {r["source_b"] for r in rows if r["source_b"] != "NA"}
        assert a_labels == {"my_vcf"}
        assert b_labels == {"my_gfa"}


# ---------------------------------------------------------------------------
# TestRunCompareJson
# ---------------------------------------------------------------------------

class TestRunCompareJson:
    def test_json_has_expected_keys(
        self, tmp_path: Path, vcf_hits: Path, gfa_hits: Path
    ) -> None:
        outdir = tmp_path / "out"
        run_compare(vcf_hits, gfa_hits, outdir, default_config())
        meta = json.loads((outdir / "compare.json").read_text())
        for key in ("tool", "hits_a", "hits_b", "source_a", "source_b",
                    "n_rows_a", "n_rows_b", "n_compare_rows", "config", "timestamp"):
            assert key in meta

    def test_json_row_counts_match(
        self, tmp_path: Path, vcf_hits: Path, gfa_hits: Path
    ) -> None:
        outdir = tmp_path / "out"
        rows = run_compare(vcf_hits, gfa_hits, outdir, default_config())
        meta = json.loads((outdir / "compare.json").read_text())
        assert meta["n_rows_a"] == 3
        assert meta["n_rows_b"] == 3
        assert meta["n_compare_rows"] == len(rows)


# ---------------------------------------------------------------------------
# TestRunCompareSummary
# ---------------------------------------------------------------------------

class TestRunCompareSummary:
    def test_summary_covers_all_match_classes(
        self, tmp_path: Path, vcf_hits: Path, gfa_hits: Path
    ) -> None:
        outdir = tmp_path / "out"
        run_compare(vcf_hits, gfa_hits, outdir, default_config())
        summary_rows = read_tsv(outdir / "compare_summary.tsv")
        classes = {r["match_class"] for r in summary_rows}
        assert classes == {mc.value for mc in MatchClass}

    def test_summary_pct_sums_to_100(
        self, tmp_path: Path, vcf_hits: Path, gfa_hits: Path
    ) -> None:
        outdir = tmp_path / "out"
        run_compare(vcf_hits, gfa_hits, outdir, default_config())
        summary_rows = read_tsv(outdir / "compare_summary.tsv")
        total = sum(float(r["pct_total"]) for r in summary_rows)
        assert total == pytest.approx(100.0, abs=0.5)


# ---------------------------------------------------------------------------
# TestRunCompareConfig
# ---------------------------------------------------------------------------

class TestRunCompareConfig:
    def test_low_overlap_threshold_increases_matches(self, tmp_path: Path) -> None:
        vcf = _write_hits(tmp_path / "v" / "hits.tsv", [
            _hit("PPX000001", "chr1", 100, 200),
        ])
        gfa = _write_hits(tmp_path / "g" / "hits.tsv", [
            _hit("GPX000001", "chr1", 150, 250),  # overlap = 50/150 ≈ 0.33
        ])
        # Strict config: threshold=0.5, no breakpoint fallback → no match
        cfg_strict = default_config().model_copy(
            update={"compare": default_config().compare.model_copy(
                update={"min_reciprocal_overlap": 0.5, "breakpoint_tolerance_bp": 0}
            )}
        )
        rows_strict = run_compare(
            vcf, gfa, tmp_path / "out1", cfg_strict,
            write_compare_tsv=False, write_summary_tsv=False, write_json=False,
        )
        assert any(r["match_class"] == MatchClass.SOURCE_SPECIFIC.value for r in rows_strict)

        # Low threshold: threshold=0.2 → should match (overlap ≈ 0.33 > 0.2)
        cfg_low = default_config().model_copy(
            update={"compare": default_config().compare.model_copy(
                update={"min_reciprocal_overlap": 0.2, "breakpoint_tolerance_bp": 0}
            )}
        )
        rows_low = run_compare(
            vcf, gfa, tmp_path / "out2", cfg_low,
            write_compare_tsv=False, write_summary_tsv=False, write_json=False,
        )
        assert any(r["match_class"] in (
            MatchClass.SUPPORTED.value, MatchClass.PARTIALLY_SUPPORTED.value
        ) for r in rows_low)

    def test_breakpoint_tolerance_enables_near_miss_match(self, tmp_path: Path) -> None:
        vcf = _write_hits(tmp_path / "v" / "hits.tsv", [
            _hit("PPX000001", "chr1", 100, 200),
        ])
        gfa = _write_hits(tmp_path / "g" / "hits.tsv", [
            _hit("GPX000001", "chr1", 210, 310),  # gap = 10 bp
        ])
        cfg = default_config().model_copy(
            update={"compare": default_config().compare.model_copy(
                update={"breakpoint_tolerance_bp": 50}
            )}
        )
        rows = run_compare(
            vcf, gfa, tmp_path / "out", cfg,
            write_compare_tsv=False, write_summary_tsv=False, write_json=False,
        )
        matched = [r for r in rows if r["locus_id_a"] != "NA" and r["locus_id_b"] != "NA"]
        assert len(matched) == 1


# ---------------------------------------------------------------------------
# TestCompareCli
# ---------------------------------------------------------------------------

class TestCompareCli:
    def test_successful_run(
        self,
        tmp_path: Path,
        vcf_hits: Path,
        gfa_hits: Path,
        runner: CliRunner,
    ) -> None:
        outdir = tmp_path / "cli_out"
        result = runner.invoke(app, [
            "compare",
            "--hits-a", str(vcf_hits),
            "--hits-b", str(gfa_hits),
            "--outdir", str(outdir),
        ])
        assert result.exit_code == 0, result.output
        assert (outdir / "compare.tsv").exists()

    def test_missing_hits_a_exits_nonzero(
        self,
        tmp_path: Path,
        gfa_hits: Path,
        runner: CliRunner,
    ) -> None:
        result = runner.invoke(app, [
            "compare",
            "--hits-a", str(tmp_path / "nonexistent.tsv"),
            "--hits-b", str(gfa_hits),
            "--outdir", str(tmp_path / "out"),
        ])
        assert result.exit_code != 0

    def test_missing_hits_b_exits_nonzero(
        self,
        tmp_path: Path,
        vcf_hits: Path,
        runner: CliRunner,
    ) -> None:
        result = runner.invoke(app, [
            "compare",
            "--hits-a", str(vcf_hits),
            "--hits-b", str(tmp_path / "nonexistent.tsv"),
            "--outdir", str(tmp_path / "out"),
        ])
        assert result.exit_code != 0

    def test_min_reciprocal_overlap_override(
        self,
        tmp_path: Path,
        vcf_hits: Path,
        gfa_hits: Path,
        runner: CliRunner,
    ) -> None:
        outdir = tmp_path / "out"
        result = runner.invoke(app, [
            "compare",
            "--hits-a", str(vcf_hits),
            "--hits-b", str(gfa_hits),
            "--min-reciprocal-overlap", "0.1",
            "--outdir", str(outdir),
        ])
        assert result.exit_code == 0, result.output
        meta = json.loads((outdir / "compare.json").read_text())
        assert meta["config"]["min_reciprocal_overlap"] == pytest.approx(0.1)

    def test_explicit_source_labels(
        self,
        tmp_path: Path,
        vcf_hits: Path,
        gfa_hits: Path,
        runner: CliRunner,
    ) -> None:
        outdir = tmp_path / "out"
        result = runner.invoke(app, [
            "compare",
            "--hits-a", str(vcf_hits),
            "--hits-b", str(gfa_hits),
            "--source-a", "myvcf",
            "--source-b", "mygfa",
            "--outdir", str(outdir),
        ])
        assert result.exit_code == 0, result.output
        rows = read_tsv(outdir / "compare.tsv")
        a_labels = {r["source_a"] for r in rows if r["source_a"] != "NA"}
        assert a_labels == {"myvcf"}

    def test_no_write_flags(
        self,
        tmp_path: Path,
        vcf_hits: Path,
        gfa_hits: Path,
        runner: CliRunner,
    ) -> None:
        outdir = tmp_path / "out"
        result = runner.invoke(app, [
            "compare",
            "--hits-a", str(vcf_hits),
            "--hits-b", str(gfa_hits),
            "--no-write-compare-tsv",
            "--no-write-summary-tsv",
            "--no-write-json",
            "--outdir", str(outdir),
        ])
        assert result.exit_code == 0, result.output
        assert not (outdir / "compare.tsv").exists()
