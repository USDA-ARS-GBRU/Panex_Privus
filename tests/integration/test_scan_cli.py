"""CLI integration tests for ``privy scan``."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

from typer.testing import CliRunner

from privy.cli.main import app
from privy.io.gfa import default_gfa_index_path

runner = CliRunner()
GFA_PATH = Path(__file__).parent.parent / "data" / "small_cohort.gfa"


def _copy_gfa(tmp_path: Path) -> Path:
    gfa = tmp_path / "small_cohort.gfa"
    gfa.write_bytes(GFA_PATH.read_bytes())
    return gfa


def test_scan_cli_runs_end_to_end(indexed_vcf: Path, tmp_path: Path) -> None:
    outdir = tmp_path / "cli-out"
    result = runner.invoke(
        app,
        [
            "scan",
            "--vcf",
            str(indexed_vcf),
            "--targets",
            "T1",
            "--targets",
            "T2",
            "--off-targets",
            "O1",
            "--off-targets",
            "O2",
            "--off-targets",
            "O3",
            "--outdir",
            str(outdir),
        ],
    )

    assert result.exit_code == 0, result.stdout
    assert (outdir / "vcf" / "hits.tsv").exists()
    assert (outdir / "vcf" / "run.json").exists()


def test_scan_cli_accepts_grouped_cohort_values(
    indexed_vcf: Path, tmp_path: Path
) -> None:
    outdir = tmp_path / "cli-grouped-cohort-out"
    result = runner.invoke(
        app,
        [
            "scan",
            "--vcf",
            str(indexed_vcf),
            "--targets",
            "T1",
            "T2",
            "--off-targets",
            "O1",
            "O2",
            "O3",
            "--outdir",
            str(outdir),
        ],
    )

    assert result.exit_code == 0, result.output
    run_data = json.loads((outdir / "vcf" / "run.json").read_text())
    assert run_data["cohort"]["targets"] == ["T1", "T2"]
    assert run_data["cohort"]["off_targets"] == ["O1", "O2", "O3"]


def test_scan_cli_applies_boolean_scan_override(indexed_vcf: Path, tmp_path: Path) -> None:
    outdir = tmp_path / "vcf-no-multiallelic-out"
    result = runner.invoke(
        app,
        [
            "scan",
            "--vcf",
            str(indexed_vcf),
            "--targets",
            "T1",
            "--targets",
            "T2",
            "--off-targets",
            "O1",
            "--off-targets",
            "O2",
            "--off-targets",
            "O3",
            "--no-allow-multiallelic",
            "--outdir",
            str(outdir),
        ],
    )

    assert result.exit_code == 0, result.output
    run_data = json.loads((outdir / "vcf" / "run.json").read_text())
    assert run_data["config"]["scan"]["allow_multiallelic"] is False
    hits_lines = (outdir / "vcf" / "hits.tsv").read_text().strip().splitlines()
    assert len(hits_lines) == 7  # header + 6 hits


def test_scan_cli_requires_a_complete_cohort(indexed_vcf: Path, tmp_path: Path) -> None:
    result = runner.invoke(
        app,
        [
            "scan",
            "--vcf",
            str(indexed_vcf),
            "--targets",
            "T1",
            "--outdir",
            str(tmp_path / "cli-out"),
        ],
    )

    assert result.exit_code == 1
    assert "Cohort is incomplete" in result.output


def test_scan_cli_runs_gfa_backend_end_to_end(tmp_path: Path) -> None:
    outdir = tmp_path / "gfa-cli-out"
    result = runner.invoke(
        app,
        [
            "scan",
            "--gfa",
            str(GFA_PATH),
            "--targets",
            "T1",
            "--targets",
            "T2",
            "--off-targets",
            "O1",
            "--off-targets",
            "O2",
            "--off-targets",
            "O3",
            "--outdir",
            str(outdir),
        ],
    )

    assert result.exit_code == 0, result.output
    assert (outdir / "gfa" / "hits.tsv").exists()
    assert (outdir / "gfa" / "run.json").exists()


def test_scan_cli_runs_gfa_gz_backend_end_to_end(tmp_path: Path) -> None:
    gfa_gz = tmp_path / "small_cohort.gfa.gz"
    with gzip.open(gfa_gz, "wb") as fh:
        fh.write(GFA_PATH.read_bytes())

    outdir = tmp_path / "gfa-gz-cli-out"
    result = runner.invoke(
        app,
        [
            "scan",
            "--gfa",
            str(gfa_gz),
            "--targets",
            "T1",
            "T2",
            "--off-targets",
            "O1",
            "O2",
            "O3",
            "--outdir",
            str(outdir),
        ],
    )

    assert result.exit_code == 0, result.output
    assert (outdir / "gfa" / "hits.tsv").exists()
    assert (outdir / "gfa" / "run.json").exists()


def test_index_gfa_cli_writes_default_sidecar(tmp_path: Path) -> None:
    gfa = _copy_gfa(tmp_path)
    index_path = default_gfa_index_path(gfa)

    result = runner.invoke(app, ["index", "gfa", "--gfa", str(gfa)])

    assert result.exit_code == 0, result.output
    assert index_path.exists()
    assert "coordinate segments" in result.output


def test_scan_cli_uses_explicit_gfa_index(tmp_path: Path) -> None:
    gfa = _copy_gfa(tmp_path)
    index_path = tmp_path / "custom.privy.gfaidx"
    index_result = runner.invoke(
        app,
        ["index", "gfa", "--gfa", str(gfa), "--out", str(index_path)],
    )
    assert index_result.exit_code == 0, index_result.output

    outdir = tmp_path / "gfa-index-cli-out"
    result = runner.invoke(
        app,
        [
            "scan",
            "--gfa",
            str(gfa),
            "--gfa-index",
            str(index_path),
            "--targets",
            "T1",
            "T2",
            "--off-targets",
            "O1",
            "O2",
            "O3",
            "--outdir",
            str(outdir),
        ],
    )

    assert result.exit_code == 0, result.output
    run_data = json.loads((outdir / "gfa" / "run.json").read_text())
    assert run_data["inputs"]["gfa_index"] == str(index_path)


def test_scan_cli_auto_detects_gfa_index(tmp_path: Path) -> None:
    gfa = _copy_gfa(tmp_path)
    index_path = default_gfa_index_path(gfa)
    index_result = runner.invoke(app, ["index", "gfa", "--gfa", str(gfa)])
    assert index_result.exit_code == 0, index_result.output

    outdir = tmp_path / "gfa-auto-index-cli-out"
    result = runner.invoke(
        app,
        [
            "scan",
            "--gfa",
            str(gfa),
            "--targets",
            "T1",
            "T2",
            "--off-targets",
            "O1",
            "O2",
            "O3",
            "--outdir",
            str(outdir),
        ],
    )

    assert result.exit_code == 0, result.output
    run_data = json.loads((outdir / "gfa" / "run.json").read_text())
    assert run_data["inputs"]["gfa_index"] == str(index_path)


def test_scan_cli_applies_gfa_min_segment_length_override(tmp_path: Path) -> None:
    outdir = tmp_path / "gfa-minlen-out"
    result = runner.invoke(
        app,
        [
            "scan",
            "--gfa",
            str(GFA_PATH),
            "--targets",
            "T1",
            "--targets",
            "T2",
            "--off-targets",
            "O1",
            "--off-targets",
            "O2",
            "--off-targets",
            "O3",
            "--min-segment-length",
            "11",
            "--outdir",
            str(outdir),
        ],
    )

    assert result.exit_code == 0, result.output
    run_data = json.loads((outdir / "gfa" / "run.json").read_text())
    assert run_data["config"]["gfa"]["min_segment_length"] == 11
    hits_lines = (outdir / "gfa" / "hits.tsv").read_text().strip().splitlines()
    assert len(hits_lines) == 1  # header only


def test_scan_cli_runs_vcf_and_gfa_then_compares(
    indexed_vcf: Path, tmp_path: Path
) -> None:
    outdir = tmp_path / "combined-out"
    result = runner.invoke(
        app,
        [
            "scan",
            "--vcf",
            str(indexed_vcf),
            "--gfa",
            str(GFA_PATH),
            "--targets",
            "T1",
            "--targets",
            "T2",
            "--off-targets",
            "O1",
            "--off-targets",
            "O2",
            "--off-targets",
            "O3",
            "--outdir",
            str(outdir),
        ],
    )

    assert result.exit_code == 0, result.output
    assert (outdir / "vcf" / "hits.tsv").exists()
    assert (outdir / "gfa" / "hits.tsv").exists()
    assert (outdir / "compare" / "compare.tsv").exists()
    assert (outdir / "compare" / "compare_summary.tsv").exists()
    compare_data = json.loads((outdir / "compare" / "compare.json").read_text())
    assert compare_data["source_a"] == "vcf"
    assert compare_data["source_b"] == "gfa"


def test_scan_cli_loads_cohort_from_yaml_file(tmp_path: Path) -> None:
    outdir = tmp_path / "gfa-yaml-cohort-out"
    cohort_yaml = tmp_path / "cohort.yaml"
    cohort_yaml.write_text(
        "targets: [T1, T2]\n"
        "off_targets: [O1, O2, O3]\n"
    )

    result = runner.invoke(
        app,
        [
            "scan",
            "--gfa",
            str(GFA_PATH),
            "--cohort-file",
            str(cohort_yaml),
            "--outdir",
            str(outdir),
        ],
    )

    assert result.exit_code == 0, result.output
    assert (outdir / "gfa" / "hits.tsv").exists()


def test_scan_cli_loads_cohort_from_tsv_file(tmp_path: Path) -> None:
    outdir = tmp_path / "gfa-tsv-cohort-out"
    cohort_tsv = tmp_path / "cohort.tsv"
    cohort_tsv.write_text(
        "sample_id\tcohort_role\n"
        "T1\ttarget\n"
        "T2\ttarget\n"
        "O1\toff_target\n"
        "O2\toff_target\n"
        "O3\toff_target\n"
    )

    result = runner.invoke(
        app,
        [
            "scan",
            "--gfa",
            str(GFA_PATH),
            "--cohort-file",
            str(cohort_tsv),
            "--outdir",
            str(outdir),
        ],
    )

    assert result.exit_code == 0, result.output
    assert (outdir / "gfa" / "hits.tsv").exists()
