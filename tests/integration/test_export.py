"""Integration tests for privy export."""

from __future__ import annotations

import json
from pathlib import Path

from typer.testing import CliRunner

from privy.backends.export import (
    run_export,
    write_hits_bed,
    write_hits_gff3,
    write_regions_bed,
    write_regions_gff3,
)
from privy.cli.main import app
from privy.io.tsv import HITS_COLUMNS, REGIONS_COLUMNS, TsvWriter


def _hit_row(locus_id: str, start: int, end: int, score: str = "0.75") -> dict[str, str]:
    return {
        "locus_id": locus_id,
        "contig": "chr1",
        "start": str(start),
        "end": str(end),
        "variant_type": "snp",
        "allele_key": f"chr1:{start + 1}:A:T",
        "target_support_n": "2",
        "target_total_n": "2",
        "offtarget_support_n": "0",
        "offtarget_total_n": "3",
        "target_missing_n": "0",
        "offtarget_missing_n": "0",
        "strictness_class": "strict_complete",
        "discovery_score": "0.75",
        "support_score": "0.0",
        "penalty_score": "0.0",
        "final_score": score,
    }


def _region_row(region_id: str, start: int, end: int, score: str = "0.5") -> dict[str, str]:
    return {
        "region_id": region_id,
        "contig": "chr1",
        "start": str(start),
        "end": str(end),
        "n_loci": "2",
        "variant_types": "snp",
        "dominant_strictness_class": "strict_complete",
        "target_consistency": "1.0",
        "offtarget_exclusion": "1.0",
        "final_score": score,
    }


def _write_hits(path: Path) -> Path:
    with TsvWriter(path, HITS_COLUMNS) as writer:
        writer.write_rows([
            _hit_row("PPX000001", 99, 100, "0.75"),
            _hit_row("PPX000002", 199, 200, "1.25"),
        ])
    return path


def _write_regions(path: Path) -> Path:
    with TsvWriter(path, REGIONS_COLUMNS) as writer:
        writer.write_rows([
            _region_row("REGION000001", 99, 200, "0.5"),
        ])
    return path


class TestBedWriters:
    def test_write_hits_bed(self, tmp_path: Path) -> None:
        out = tmp_path / "hits.bed"
        write_hits_bed([_hit_row("PPX000001", 99, 100, "0.75")], out)

        lines = out.read_text(encoding="utf-8").splitlines()
        assert lines[0].startswith('track name="Panex Privus hits"')
        assert lines[1].split("\t") == [
            "chr1",
            "99",
            "100",
            "PPX000001",
            "750",
            ".",
            "strict_complete",
            "snp",
            "chr1:100:A:T",
            "0.75",
        ]

    def test_write_regions_bed(self, tmp_path: Path) -> None:
        out = tmp_path / "regions.bed"
        write_regions_bed([_region_row("REGION000001", 99, 200, "0.5")], out)

        lines = out.read_text(encoding="utf-8").splitlines()
        assert lines[1].split("\t") == [
            "chr1",
            "99",
            "200",
            "REGION000001",
            "500",
            ".",
            "strict_complete",
            "snp",
            "2",
            "0.5",
        ]

    def test_no_header(self, tmp_path: Path) -> None:
        out = tmp_path / "hits.bed"
        write_hits_bed([_hit_row("PPX000001", 99, 100)], out, include_header=False)
        assert out.read_text(encoding="utf-8").startswith("chr1\t99\t100")


class TestGff3Writers:
    def test_write_hits_gff3(self, tmp_path: Path) -> None:
        out = tmp_path / "hits.gff3"
        write_hits_gff3([_hit_row("PPX000001", 99, 100, "0.75")], out)

        lines = out.read_text(encoding="utf-8").splitlines()
        assert lines[0] == "##gff-version 3"
        assert lines[1].split("\t") == [
            "chr1",
            "privy",
            "sequence_variant",
            "100",
            "100",
            "0.75",
            ".",
            ".",
            "ID=PPX000001;Name=PPX000001;strictness_class=strict_complete;"
            "variant_type=snp;allele_key=chr1:100:A:T;final_score=0.75",
        ]

    def test_write_regions_gff3(self, tmp_path: Path) -> None:
        out = tmp_path / "regions.gff3"
        write_regions_gff3([_region_row("REGION000001", 99, 200, "0.5")], out)

        lines = out.read_text(encoding="utf-8").splitlines()
        assert lines[1].split("\t") == [
            "chr1",
            "privy",
            "region",
            "100",
            "200",
            "0.5",
            ".",
            ".",
            "ID=REGION000001;Name=REGION000001;"
            "dominant_strictness_class=strict_complete;variant_types=snp;"
            "n_loci=2;final_score=0.5",
        ]

    def test_gff3_attributes_are_url_escaped(self, tmp_path: Path) -> None:
        out = tmp_path / "hits.gff3"
        row = _hit_row("PPX000001", 99, 100, "0.75")
        row["allele_key"] = "chr1:100:A:<DEL with space>"

        write_hits_gff3([row], out)

        assert "allele_key=chr1:100:A:%3CDEL%20with%20space%3E" in out.read_text(
            encoding="utf-8"
        )


class TestRunExport:
    def test_exports_hits_regions_and_metadata(self, tmp_path: Path) -> None:
        hits = _write_hits(tmp_path / "hits.tsv")
        regions = _write_regions(tmp_path / "regions.tsv")
        outdir = tmp_path / "exported"

        written = run_export(hits, regions, outdir)

        assert outdir / "hits.bed" in written
        assert outdir / "regions.bed" in written
        assert outdir / "export.json" in written
        meta = json.loads((outdir / "export.json").read_text(encoding="utf-8"))
        assert meta["format"] == "bed"
        assert meta["kind"] == "both"

    def test_exports_hits_only(self, tmp_path: Path) -> None:
        hits = _write_hits(tmp_path / "hits.tsv")
        outdir = tmp_path / "exported"

        run_export(hits, None, outdir, export_kind="hits")

        assert (outdir / "hits.bed").exists()
        assert not (outdir / "regions.bed").exists()

    def test_exports_gff3(self, tmp_path: Path) -> None:
        hits = _write_hits(tmp_path / "hits.tsv")
        regions = _write_regions(tmp_path / "regions.tsv")
        outdir = tmp_path / "exported"

        written = run_export(hits, regions, outdir, export_format="gff3")

        assert outdir / "hits.gff3" in written
        assert outdir / "regions.gff3" in written
        meta = json.loads((outdir / "export.json").read_text(encoding="utf-8"))
        assert meta["format"] == "gff3"

    def test_missing_hits_for_hits_export_raises(self, tmp_path: Path) -> None:
        outdir = tmp_path / "exported"
        try:
            run_export(None, None, outdir, export_kind="hits")
        except ValueError as exc:
            assert "--hits is required" in str(exc)
        else:
            raise AssertionError("Expected ValueError")

    def test_score_is_clamped_to_bed_range(self, tmp_path: Path) -> None:
        out = tmp_path / "hits.bed"
        write_hits_bed([_hit_row("PPX000001", 99, 100, "1.25")], out)
        assert out.read_text(encoding="utf-8").splitlines()[1].split("\t")[4] == "1000"


class TestExportCli:
    def setup_method(self) -> None:
        self.runner = CliRunner()

    def test_cli_runs_successfully(self, tmp_path: Path) -> None:
        hits = _write_hits(tmp_path / "hits.tsv")
        regions = _write_regions(tmp_path / "regions.tsv")
        outdir = tmp_path / "exported"

        result = self.runner.invoke(app, [
            "export",
            "--hits", str(hits),
            "--regions", str(regions),
            "--outdir", str(outdir),
        ])

        assert result.exit_code == 0, result.output
        assert (outdir / "hits.bed").exists()
        assert (outdir / "regions.bed").exists()

    def test_cli_hits_only(self, tmp_path: Path) -> None:
        hits = _write_hits(tmp_path / "hits.tsv")
        outdir = tmp_path / "exported"

        result = self.runner.invoke(app, [
            "export",
            "--hits", str(hits),
            "--kind", "hits",
            "--no-include-header",
            "--outdir", str(outdir),
        ])

        assert result.exit_code == 0, result.output
        text = (outdir / "hits.bed").read_text(encoding="utf-8")
        assert text.startswith("chr1\t99\t100")

    def test_cli_gff3(self, tmp_path: Path) -> None:
        hits = _write_hits(tmp_path / "hits.tsv")
        outdir = tmp_path / "exported"

        result = self.runner.invoke(app, [
            "export",
            "--hits", str(hits),
            "--kind", "hits",
            "--format", "gff3",
            "--outdir", str(outdir),
        ])

        assert result.exit_code == 0, result.output
        assert (outdir / "hits.gff3").exists()

    def test_cli_missing_requested_input_exits_nonzero(self, tmp_path: Path) -> None:
        result = self.runner.invoke(app, [
            "export",
            "--hits", str(tmp_path / "missing.tsv"),
            "--kind", "hits",
            "--outdir", str(tmp_path / "exported"),
        ])

        assert result.exit_code == 1
        assert "--hits not found" in result.output
