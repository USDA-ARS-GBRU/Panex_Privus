"""Integration tests for ``privy popgen``."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from typer.testing import CliRunner

from privy.cli.main import app
from privy.synthetic import microhaplotype_pangenome

runner = CliRunner()


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


class TestPopgenCli:
    def test_writes_tables_with_diagnostic_marker(self, tmp_path):
        gfa = microhaplotype_pangenome(seg_len=10).write(tmp_path / "g.gfa")
        out = tmp_path / "out"
        result = runner.invoke(
            app,
            [
                "popgen", "--gfa", str(gfa), "--reference", "sample0#0#chr1",
                "--targets", "sample0,sample1", "--off-targets", "sample2,sample3",
                "--outdir", str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        loci = _read_tsv(out / "popgen_loci.tsv")
        assert len(loci) == 1
        row = loci[0]
        assert row["n_alleles"] == "2"
        assert row["gst"] == "1.0000"          # fully diagnostic target-private locus
        assert row["is_diagnostic"] == "True"
        assert "fis" in row   # F_IS column present
        meta = json.loads((out / "popgen.json").read_text())
        assert meta["genome_wide_fst"] == 1.0
        assert meta["n_diagnostic_loci"] == 1
        # private-allele metrics
        assert meta["private_allele_counts"] == {"target": 1, "offtarget": 1}
        assert meta["private_allelic_richness"]["target"] >= 0.0
        # GP-ready exports also written
        assert (out / "dosage_matrix.tsv").exists()
        assert (out / "grm.tsv").exists()
        grm = _read_tsv(out / "grm.tsv")
        assert grm[0]["sample"]   # labelled square matrix with a sample column

    def test_requires_both_cohorts_present_in_graph(self, tmp_path):
        gfa = microhaplotype_pangenome().write(tmp_path / "g.gfa")
        out = tmp_path / "out"
        result = runner.invoke(
            app,
            ["popgen", "--gfa", str(gfa), "--reference", "sample0#0#chr1",
             "--targets", "sample0", "--off-targets", "ghost", "--outdir", str(out)],
        )
        assert result.exit_code == 1
        assert "off-targets" in result.output or "resolve" in result.output

    def test_unknown_reference_errors(self, tmp_path):
        gfa = microhaplotype_pangenome().write(tmp_path / "g.gfa")
        result = runner.invoke(
            app,
            ["popgen", "--gfa", str(gfa), "--reference", "ghost#0#chr1",
             "--targets", "sample0", "--off-targets", "sample2"],
        )
        assert result.exit_code == 1
        assert "not found" in result.output

    def test_missing_gfa_errors(self, tmp_path):
        result = runner.invoke(
            app,
            ["popgen", "--gfa", str(tmp_path / "x.gfa"), "--reference", "r",
             "--targets", "a", "--off-targets", "b"],
        )
        assert result.exit_code == 1
        assert "not found" in result.output
