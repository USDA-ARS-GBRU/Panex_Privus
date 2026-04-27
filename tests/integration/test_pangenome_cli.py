"""CLI integration tests for ``privy pangenome``."""

from __future__ import annotations

import json
from pathlib import Path

from typer.testing import CliRunner

from privy.cli.main import app
from privy.io.tsv import read_tsv

runner = CliRunner()
GFA_PATH = Path(__file__).parent.parent / "data" / "small_cohort.gfa"


def test_pangenome_cli_runs_gfa_with_inferred_offtargets(tmp_path: Path) -> None:
    outdir = tmp_path / "pangenome-out"
    result = runner.invoke(
        app,
        [
            "pangenome",
            "--gfa",
            str(GFA_PATH),
            "--targets",
            "T1",
            "--targets",
            "T2",
            "--permutations",
            "2",
            "--outdir",
            str(outdir),
        ],
    )

    assert result.exit_code == 0, result.output
    assert (outdir / "feature_summary.tsv").exists()
    assert (outdir / "coverage_histogram.tsv").exists()
    assert (outdir / "composition.tsv").exists()
    assert (outdir / "growth_curves.tsv").exists()
    assert (outdir / "pangenome_growth.png").exists()
    assert (outdir / "pangenome_coverage.png").exists()
    assert (outdir / "pangenome_composition.png").exists()
    data = json.loads((outdir / "pangenome.json").read_text())
    assert data["samples"]["off_target"] == ["O1", "O2", "O3"]
    assert "pangenome_growth.png" in data["outputs"]


def test_pangenome_cli_accepts_target_and_offtarget_list_files(tmp_path: Path) -> None:
    targets_file = tmp_path / "targets.txt"
    offtargets_file = tmp_path / "offtargets.txt"
    targets_file.write_text("T1\nT2\n", encoding="utf-8")
    offtargets_file.write_text("O1\nO2\nO3\n", encoding="utf-8")
    outdir = tmp_path / "pangenome-files-out"

    result = runner.invoke(
        app,
        [
            "pangenome",
            "--gfa",
            str(GFA_PATH),
            "--targets-file",
            str(targets_file),
            "--off-targets-file",
            str(offtargets_file),
            "--permutations",
            "1",
            "--outdir",
            str(outdir),
        ],
    )

    assert result.exit_code == 0, result.output
    rows = read_tsv(outdir / "feature_summary.tsv")
    s2_target = next(row for row in rows if row["feature_id"] == "s2_target")
    assert s2_target["target_private"] == "True"


def test_pangenome_cli_can_skip_plots(tmp_path: Path) -> None:
    outdir = tmp_path / "pangenome-no-plots"
    result = runner.invoke(
        app,
        [
            "pangenome",
            "--gfa",
            str(GFA_PATH),
            "--targets",
            "T1",
            "--targets",
            "T2",
            "--permutations",
            "1",
            "--no-plots",
            "--outdir",
            str(outdir),
        ],
    )

    assert result.exit_code == 0, result.output
    assert not (outdir / "pangenome_growth.png").exists()
    data = json.loads((outdir / "pangenome.json").read_text())
    assert data["parameters"]["write_plots"] is False
