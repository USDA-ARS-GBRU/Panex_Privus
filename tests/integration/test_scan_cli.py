"""CLI integration tests for ``privy scan``."""

from __future__ import annotations

from pathlib import Path

from typer.testing import CliRunner

from privy.cli.main import app

runner = CliRunner()


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
    assert (outdir / "hits.tsv").exists()
    assert (outdir / "run.json").exists()


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
