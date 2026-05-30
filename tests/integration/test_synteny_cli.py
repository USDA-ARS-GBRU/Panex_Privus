"""Integration tests for ``privy synteny``."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from typer.testing import CliRunner

from privy.cli.main import app
from privy.io.paf import PafRecord, write_paf
from privy.synthetic import (
    inversion_pangenome,
    presence_absence_pangenome,
)

runner = CliRunner()


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


class TestSyntenyCli:
    def test_writes_blocks_and_regions(self, tmp_path):
        gfa = inversion_pangenome(seg_len=10).write(tmp_path / "g.gfa")
        out = tmp_path / "out"
        result = runner.invoke(
            app,
            ["synteny", "--gfa", str(gfa), "--reference", "sample0#0#chr1",
             "--outdir", str(out)],
        )
        assert result.exit_code == 0, result.output
        blocks = _read_tsv(out / "synteny_blocks.tsv")
        # sample3 contributes an inversion; ensure at least one inversion block exists
        assert any(b["block_type"] == "inversion" for b in blocks)
        meta = json.loads((out / "synteny.json").read_text())
        assert meta["reference"] == "sample0#0#chr1"
        assert meta["block_type_counts"].get("inversion", 0) >= 1

    def test_private_region_flagged(self, tmp_path):
        gfa = presence_absence_pangenome(seg_len=10).write(tmp_path / "g.gfa")
        out = tmp_path / "out"
        result = runner.invoke(
            app,
            [
                "synteny", "--gfa", str(gfa), "--reference", "sample0#0#chr1",
                "--targets", "sample1", "--off-targets", "sample2,sample3",
                "--outdir", str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        regions = _read_tsv(out / "synteny_regions.tsv")
        assert any(r["target_private"] == "True" for r in regions)
        meta = json.loads((out / "synteny.json").read_text())
        assert meta["n_target_private_regions"] >= 1

    def test_no_cohorts_leaves_privacy_na(self, tmp_path):
        gfa = presence_absence_pangenome().write(tmp_path / "g.gfa")
        out = tmp_path / "out"
        result = runner.invoke(
            app,
            ["synteny", "--gfa", str(gfa), "--reference", "sample0#0#chr1",
             "--outdir", str(out)],
        )
        assert result.exit_code == 0, result.output
        regions = _read_tsv(out / "synteny_regions.tsv")
        assert all(r["target_private"] == "NA" for r in regions)

    def test_unknown_reference_errors(self, tmp_path):
        gfa = inversion_pangenome().write(tmp_path / "g.gfa")
        result = runner.invoke(
            app, ["synteny", "--gfa", str(gfa), "--reference", "ghost#0#chr1"]
        )
        assert result.exit_code == 1
        assert "not found" in result.output

    def test_missing_gfa_errors(self, tmp_path):
        result = runner.invoke(
            app, ["synteny", "--gfa", str(tmp_path / "nope.gfa"), "--reference", "x"]
        )
        assert result.exit_code == 1
        assert "not found" in result.output

    def test_requires_exactly_one_source(self, tmp_path):
        gfa = inversion_pangenome().write(tmp_path / "g.gfa")
        # neither
        r1 = runner.invoke(app, ["synteny", "--outdir", str(tmp_path / "o")])
        assert r1.exit_code == 1 and "exactly one" in r1.output
        # both
        paf = tmp_path / "a.paf"
        write_paf([PafRecord("q#0#c", 100, 0, 80, "+", "t#0#c", 100, 0, 80, 80, 80, 60)], paf)
        r2 = runner.invoke(
            app, ["synteny", "--gfa", str(gfa), "--paf", str(paf), "--outdir", str(tmp_path / "o")]
        )
        assert r2.exit_code == 1 and "exactly one" in r2.output

    def test_gfa_mode_requires_reference(self, tmp_path):
        gfa = inversion_pangenome().write(tmp_path / "g.gfa")
        result = runner.invoke(app, ["synteny", "--gfa", str(gfa), "--outdir", str(tmp_path / "o")])
        assert result.exit_code == 1
        assert "--reference is required" in result.output


class TestSyntenyPafMode:
    def test_paf_chaining_writes_blocks(self, tmp_path):
        recs = [
            PafRecord("qA#0#chr1", 100_000, i * 100, i * 100 + 80, "+",
                      "tB#0#chr1", 100_000, i * 100, i * 100 + 80, 80, 80, 60)
            for i in range(5)
        ]
        paf = tmp_path / "aln.paf"
        write_paf(recs, paf)
        out = tmp_path / "out"
        result = runner.invoke(
            app, ["synteny", "--paf", str(paf), "--min-block-anchors", "3", "--outdir", str(out)]
        )
        assert result.exit_code == 0, result.output
        blocks = _read_tsv(out / "synteny_blocks.tsv")
        assert len(blocks) == 1
        assert blocks[0]["block_type"] == "collinear"
        assert blocks[0]["query_genome"] == "qA"
        # PAF mode -> privacy NA
        regions = _read_tsv(out / "synteny_regions.tsv")
        assert all(r["target_private"] == "NA" for r in regions)
