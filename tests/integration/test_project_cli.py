"""Integration tests for ``privy project``."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from typer.testing import CliRunner

from privy.cli.main import app
from privy.synthetic import collinear_pangenome

runner = CliRunner()


def _gfa(tmp_path: Path) -> Path:
    return collinear_pangenome(n_genomes=3, n_segments=6, seg_len=10).write(tmp_path / "g.gfa")


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


class TestProjectRegion:
    def test_region_projection_writes_outputs(self, tmp_path):
        gfa = _gfa(tmp_path)
        out = tmp_path / "out"
        result = runner.invoke(
            app,
            ["project", "--gfa", str(gfa), "--region", "sample0#0#chr1:15-35",
             "--outdir", str(out)],
        )
        assert result.exit_code == 0, result.output
        tsv = out / "projection.tsv"
        meta = out / "project.json"
        assert tsv.exists() and meta.exists()
        rows = _read_tsv(tsv)
        # all three collinear genomes present, spanning chr1:10-40 (s2,s3,s4)
        assert len(rows) == 3
        for row in rows:
            assert row["present"] == "True"
            assert row["contig"] == "chr1"
            assert (row["start"], row["end"]) == ("10", "40")
        assert json.loads(meta.read_text())["n_present"] == 3

    def test_to_genomes_subset(self, tmp_path):
        gfa = _gfa(tmp_path)
        out = tmp_path / "out"
        result = runner.invoke(
            app,
            [
                "project", "--gfa", str(gfa),
                "--region", "sample0#0#chr1:0-20",
                "--to-genomes", "sample2#0#chr1",
                "--outdir", str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        rows = _read_tsv(out / "projection.tsv")
        assert [r["target_path"] for r in rows] == ["sample2#0#chr1"]


class TestProjectNodeSet:
    def test_node_set_projection(self, tmp_path):
        gfa = _gfa(tmp_path)
        out = tmp_path / "out"
        result = runner.invoke(
            app,
            ["project", "--gfa", str(gfa), "--node-set", "s3,s4", "--outdir", str(out)],
        )
        assert result.exit_code == 0, result.output
        rows = _read_tsv(out / "projection.tsv")
        for row in rows:
            assert (row["start"], row["end"]) == ("20", "40")   # s3,s4 -> chr1:20-40


class TestProjectErrors:
    def test_missing_gfa(self, tmp_path):
        result = runner.invoke(
            app, ["project", "--gfa", str(tmp_path / "nope.gfa"), "--node-set", "s1"]
        )
        assert result.exit_code == 1
        assert "not found" in result.output

    def test_both_region_and_node_set(self, tmp_path):
        gfa = _gfa(tmp_path)
        result = runner.invoke(
            app,
            ["project", "--gfa", str(gfa), "--region", "sample0#0#chr1:0-10", "--node-set", "s1"],
        )
        assert result.exit_code == 1
        assert "exactly one" in result.output

    def test_neither_region_nor_node_set(self, tmp_path):
        gfa = _gfa(tmp_path)
        result = runner.invoke(app, ["project", "--gfa", str(gfa)])
        assert result.exit_code == 1
        assert "exactly one" in result.output

    def test_bad_region_format(self, tmp_path):
        gfa = _gfa(tmp_path)
        result = runner.invoke(
            app, ["project", "--gfa", str(gfa), "--region", "sample0#0#chr1-15-35"]
        )
        assert result.exit_code == 1

    def test_unknown_source_path(self, tmp_path):
        gfa = _gfa(tmp_path)
        result = runner.invoke(
            app, ["project", "--gfa", str(gfa), "--region", "ghost#0#chr1:0-10"]
        )
        assert result.exit_code == 1
        assert "not found" in result.output
