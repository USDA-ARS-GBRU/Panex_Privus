"""Integration tests for ``privy dashboard``."""

from __future__ import annotations

from pathlib import Path

from typer.testing import CliRunner

from privy.cli.main import app
from privy.synthetic import presence_absence_pangenome

runner = CliRunner()


def _synteny_dir(tmp_path: Path) -> Path:
    gfa = presence_absence_pangenome(seg_len=10).write(tmp_path / "g.gfa")
    out = tmp_path / "syn"
    r = runner.invoke(
        app,
        ["synteny", "--gfa", str(gfa), "--reference", "sample0#0#chr1",
         "--targets", "sample1", "--off-targets", "sample2,sample3", "--outdir", str(out)],
    )
    assert r.exit_code == 0, r.output
    return out


class TestDashboardCli:
    def test_builds_dashboard(self, tmp_path):
        syn = _synteny_dir(tmp_path)
        result = runner.invoke(app, ["dashboard", "--synteny", str(syn)])
        assert result.exit_code == 0, result.output
        out = syn / "synteny_dashboard.html"
        assert out.exists()
        assert '{"__privy_placeholder__": true}' not in out.read_text(encoding="utf-8")

    def test_outdir(self, tmp_path):
        syn = _synteny_dir(tmp_path)
        dash = tmp_path / "dash"
        result = runner.invoke(
            app, ["dashboard", "--synteny", str(syn), "--outdir", str(dash)]
        )
        assert result.exit_code == 0, result.output
        assert (dash / "synteny_dashboard.html").exists()

    def test_missing_dir_errors(self, tmp_path):
        result = runner.invoke(app, ["dashboard", "--synteny", str(tmp_path / "nope")])
        assert result.exit_code == 1
        assert "not found" in result.output

    def test_missing_blocks_errors(self, tmp_path):
        empty = tmp_path / "empty"
        empty.mkdir()
        result = runner.invoke(app, ["dashboard", "--synteny", str(empty)])
        assert result.exit_code == 1
        assert "not found" in result.output
