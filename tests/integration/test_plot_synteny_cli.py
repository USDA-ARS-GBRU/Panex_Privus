"""Integration tests for ``privy plot --plot-set synteny``."""

from __future__ import annotations

from pathlib import Path

from typer.testing import CliRunner

from privy.cli.main import app
from privy.synthetic import inversion_pangenome

runner = CliRunner()


def _synteny_dir(tmp_path: Path) -> Path:
    gfa = inversion_pangenome(seg_len=10).write(tmp_path / "g.gfa")
    out = tmp_path / "syn"
    result = runner.invoke(
        app, ["synteny", "--gfa", str(gfa), "--reference", "sample0#0#chr1",
               "--outdir", str(out)]
    )
    assert result.exit_code == 0, result.output
    return out


class TestPlotSyntenyCli:
    def test_renders_figures_from_synteny_dir(self, tmp_path):
        syn = _synteny_dir(tmp_path)
        result = runner.invoke(
            app, ["plot", "--plot-set", "synteny", "--input-dir", str(syn)]
        )
        assert result.exit_code == 0, result.output
        assert (syn / "riparian.png").exists()
        assert (syn / "dotplot.png").exists()

    def test_pdf_to_separate_outdir(self, tmp_path):
        syn = _synteny_dir(tmp_path)
        figs = tmp_path / "figs"
        result = runner.invoke(
            app,
            ["plot", "--plot-set", "synteny", "--input-dir", str(syn),
             "--output-format", "pdf", "--outdir", str(figs)],
        )
        assert result.exit_code == 0, result.output
        assert (figs / "riparian.pdf").exists()
        assert (figs / "dotplot.pdf").exists()

    def test_missing_blocks_tsv_errors(self, tmp_path):
        empty = tmp_path / "empty"
        empty.mkdir()
        result = runner.invoke(
            app, ["plot", "--plot-set", "synteny", "--input-dir", str(empty)]
        )
        assert result.exit_code == 1
        assert "not found" in result.output
