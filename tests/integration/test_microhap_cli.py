"""Integration tests for ``privy microhap``."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from typer.testing import CliRunner

from privy.cli.main import app
from privy.synthetic import collinear_pangenome, microhaplotype_pangenome

runner = CliRunner()


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


class TestMicrohapCli:
    def test_detects_and_writes_tables(self, tmp_path):
        gfa = microhaplotype_pangenome(seg_len=10).write(tmp_path / "g.gfa")
        out = tmp_path / "out"
        result = runner.invoke(
            app,
            ["microhap", "--gfa", str(gfa), "--reference", "sample0#0#chr1",
             "--outdir", str(out)],
        )
        assert result.exit_code == 0, result.output
        loci = _read_tsv(out / "microhaplotypes.tsv")
        assert len(loci) == 1
        assert loci[0]["n_alleles"] == "2"
        # allele matrix has one column per genome
        matrix = _read_tsv(out / "allele_matrix.tsv")
        assert "sample0#0#chr1" in matrix[0]
        meta = json.loads((out / "microhap.json").read_text())
        assert meta["n_loci"] == 1

    def test_flags_private_allele(self, tmp_path):
        gfa = microhaplotype_pangenome(seg_len=10).write(tmp_path / "g.gfa")
        out = tmp_path / "out"
        result = runner.invoke(
            app,
            [
                "microhap", "--gfa", str(gfa), "--reference", "sample0#0#chr1",
                "--targets", "sample0,sample1", "--off-targets", "sample2,sample3",
                "--outdir", str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        loci = _read_tsv(out / "microhaplotypes.tsv")
        assert loci[0]["target_private"] == "True"
        assert loci[0]["n_private_alleles"] == "1"
        meta = json.loads((out / "microhap.json").read_text())
        assert meta["n_target_private_loci"] == 1

    def test_collinear_yields_no_loci(self, tmp_path):
        gfa = collinear_pangenome(n_genomes=4, n_segments=5).write(tmp_path / "g.gfa")
        out = tmp_path / "out"
        result = runner.invoke(
            app,
            ["microhap", "--gfa", str(gfa), "--reference", "sample0#0#chr1",
             "--outdir", str(out)],
        )
        assert result.exit_code == 0, result.output
        assert _read_tsv(out / "microhaplotypes.tsv") == []

    def test_unknown_reference_errors(self, tmp_path):
        gfa = microhaplotype_pangenome().write(tmp_path / "g.gfa")
        result = runner.invoke(
            app, ["microhap", "--gfa", str(gfa), "--reference", "ghost#0#chr1"]
        )
        assert result.exit_code == 1
        assert "not found" in result.output

    def test_missing_gfa_errors(self, tmp_path):
        result = runner.invoke(
            app, ["microhap", "--gfa", str(tmp_path / "x.gfa"), "--reference", "r"]
        )
        assert result.exit_code == 1
        assert "not found" in result.output
