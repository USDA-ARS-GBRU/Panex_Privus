"""CLI integration tests for ``privy landscape``."""

from __future__ import annotations

import json
from pathlib import Path

from typer.testing import CliRunner

from privy.cli.main import app
from privy.io.tsv import read_tsv

runner = CliRunner()


def test_landscape_cli_runs_vcf_and_writes_tables_and_plots(
    indexed_vcf: Path, tmp_path: Path
) -> None:
    outdir = tmp_path / "landscape-out"
    result = runner.invoke(
        app,
        [
            "landscape",
            "--vcf",
            str(indexed_vcf),
            "--targets",
            "T1",
            "--targets",
            "T2",
            "--window-records",
            "3",
            "--step-records",
            "3",
            "--min-called-for-freq",
            "0",
            "--min-freq-values",
            "0",
            "--outdir",
            str(outdir),
        ],
    )

    assert result.exit_code == 0, result.output
    assert (outdir / "sample_windows.tsv").exists()
    assert (outdir / "windows.tsv").exists()
    assert (outdir / "background_blocks.tsv").exists()
    assert (outdir / "candidate_introgression_blocks.tsv").exists()
    assert (outdir / "similarity.tsv").exists()
    assert (outdir / "missingness_heatmap.png").exists()
    assert (outdir / "private_burden_heatmap.png").exists()
    assert (outdir / "local_background_map.png").exists()
    assert (outdir / "similarity_cluster_map.png").exists()

    rows = read_tsv(outdir / "sample_windows.tsv")
    assert len(rows) == 15
    t1_first = next(row for row in rows if row["window_id"] == "LW00000001"
                    and row["sample"] == "T1")
    assert t1_first["private_alt_n"] == "3"
    data = json.loads((outdir / "landscape.json").read_text())
    assert data["analysis"] == "landscape"
    assert data["parameters"]["window_mode"] == "records"
    assert "candidate_introgression_blocks.tsv" in data["outputs"]
    assert "local_background_map.png" in data["outputs"]


def test_landscape_cli_accepts_grouped_cohort_sample_lists(
    indexed_vcf: Path, tmp_path: Path
) -> None:
    outdir = tmp_path / "landscape-grouped-out"
    result = runner.invoke(
        app,
        [
            "landscape",
            "--vcf",
            str(indexed_vcf),
            "--targets",
            "T1",
            "T2",
            "--off-targets",
            "O1",
            "O2",
            "O3",
            "--window-records",
            "3",
            "--step-records",
            "3",
            "--min-called-for-freq",
            "0",
            "--min-freq-values",
            "0",
            "--no-plots",
            "--outdir",
            str(outdir),
        ],
    )

    assert result.exit_code == 0, result.output
    data = json.loads((outdir / "landscape.json").read_text())
    assert data["samples"]["target"] == ["T1", "T2"]
    assert data["samples"]["off_target"] == ["O1", "O2", "O3"]
    assert data["summary"]["n_sample_window_rows"] == 15


def test_landscape_cli_accepts_cohort_file(
    indexed_vcf: Path, tmp_path: Path
) -> None:
    outdir = tmp_path / "landscape-cohort-file-out"
    cohort_yaml = tmp_path / "cohort.yaml"
    cohort_yaml.write_text(
        "targets: [T1, T2]\n"
        "off_targets: [O1, O2, O3]\n",
        encoding="utf-8",
    )

    result = runner.invoke(
        app,
        [
            "landscape",
            "--vcf",
            str(indexed_vcf),
            "--cohort-file",
            str(cohort_yaml),
            "--window-records",
            "3",
            "--step-records",
            "3",
            "--min-called-for-freq",
            "0",
            "--min-freq-values",
            "0",
            "--no-plots",
            "--outdir",
            str(outdir),
        ],
    )

    assert result.exit_code == 0, result.output
    data = json.loads((outdir / "landscape.json").read_text())
    assert data["samples"]["target"] == ["T1", "T2"]
    assert data["samples"]["off_target"] == ["O1", "O2", "O3"]


def test_landscape_cli_can_skip_plots(indexed_vcf: Path, tmp_path: Path) -> None:
    outdir = tmp_path / "landscape-no-plots"
    result = runner.invoke(
        app,
        [
            "landscape",
            "--vcf",
            str(indexed_vcf),
            "--targets",
            "T1",
            "--targets",
            "T2",
            "--window-records",
            "3",
            "--step-records",
            "3",
            "--no-plots",
            "--outdir",
            str(outdir),
        ],
    )

    assert result.exit_code == 0, result.output
    assert (outdir / "sample_windows.tsv").exists()
    assert not (outdir / "local_background_map.png").exists()
    data = json.loads((outdir / "landscape.json").read_text())
    assert data["parameters"]["write_plots"] is False
