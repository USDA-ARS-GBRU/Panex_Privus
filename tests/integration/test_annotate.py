"""Integration tests for privy annotate — run_annotate() and CLI."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

from privy.backends.annotate import classify_locus, run_annotate
from privy.cli.main import app
from privy.io.gff import build_annotation_index
from privy.io.tsv import ANNOTATED_HITS_COLUMNS, ANNOTATION_SUMMARY_COLUMNS, read_tsv

GFF3_FIXTURE = Path(__file__).parent.parent / "data" / "small_cohort.gff3"

_HITS_HEADER = (
    "locus_id\tcontig\tstart\tend\tvariant_type\tallele_key\t"
    "target_support_n\ttarget_total_n\tofftarget_support_n\t"
    "offtarget_total_n\ttarget_missing_n\tofftarget_missing_n\t"
    "strictness_class\tdiscovery_score\tsupport_score\tpenalty_score\tfinal_score\n"
)


def _write_hits(path: Path, rows: list[dict]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        fh.write(_HITS_HEADER)
        fieldnames = _HITS_HEADER.strip().split("\t")
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        for row in rows:
            writer.writerow(row)
    return path


def _hit(locus_id: str, contig: str, start: int, end: int, strictness: str = "strict_complete") -> dict:
    return {
        "locus_id": locus_id, "contig": contig,
        "start": str(start), "end": str(end),
        "variant_type": "snp", "allele_key": f"{contig}:{start}:A:T",
        "target_support_n": "2", "target_total_n": "2",
        "offtarget_support_n": "0", "offtarget_total_n": "3",
        "target_missing_n": "0", "offtarget_missing_n": "0",
        "strictness_class": strictness, "discovery_score": "1.0",
        "support_score": "0.0", "penalty_score": "0.0", "final_score": "1.0",
    }


# VCF positions → 0-based intervals (SNPs: end = start + 1)
# pos 100 → [99, 100): CDS  (GeneA CDS at [90, 110))
# pos 200 → [199, 200): UTR  (GeneA five_prime_UTR at [190, 210))
# pos 300 → [299, 300): CDS  (GeneA CDS at [290, 310))
# pos 400 → [399, 400): intergenic
# pos 500 → [499, 500): intronic  (GeneB body [450, 550), exon only [450, 470))


@pytest.fixture()
def hits_file(tmp_path: Path) -> Path:
    return _write_hits(tmp_path / "hits.tsv", [
        _hit("PPX000001", "chr1", 99, 100),   # CDS
        _hit("PPX000002", "chr1", 199, 200),  # UTR
        _hit("PPX000003", "chr1", 299, 300),  # CDS
        _hit("PPX000004", "chr1", 399, 400),  # intergenic
        _hit("PPX000005", "chr1", 499, 500),  # intronic
    ])


# ---------------------------------------------------------------------------
# Class 1: classify_locus unit checks (via index)
# ---------------------------------------------------------------------------

class TestClassifyLocus:
    def setup_method(self):
        self.idx = build_annotation_index(GFF3_FIXTURE)

    def test_cds_classification(self):
        ann_class, gene_id, *_ = classify_locus(self.idx, "chr1", 99, 100)
        assert ann_class == "CDS"
        assert gene_id == "GeneA"

    def test_utr_classification(self):
        ann_class, gene_id, *_ = classify_locus(self.idx, "chr1", 199, 200)
        assert ann_class == "UTR"
        assert gene_id == "GeneA"

    def test_second_cds(self):
        ann_class, *_ = classify_locus(self.idx, "chr1", 299, 300)
        assert ann_class == "CDS"

    def test_intronic_classification(self):
        # GeneB body [450,550), exon only [450,470); pos 499→[499,500) is intronic
        ann_class, gene_id, *_ = classify_locus(self.idx, "chr1", 499, 500)
        assert ann_class == "intronic"
        assert gene_id == "GeneB"

    def test_intergenic_classification(self):
        ann_class, gene_id, strand, g_start, g_end = classify_locus(self.idx, "chr1", 699, 700)
        assert ann_class == "intergenic"
        assert gene_id == ""
        assert g_start == -1

    def test_unknown_contig_intergenic(self):
        ann_class, *_ = classify_locus(self.idx, "chrX", 100, 200)
        assert ann_class == "intergenic"

    def test_gene_strand_returned(self):
        _, _, strand, _, _ = classify_locus(self.idx, "chr1", 99, 100)
        assert strand == "+"


# ---------------------------------------------------------------------------
# Class 2: run_annotate output files
# ---------------------------------------------------------------------------

class TestRunAnnotate:
    def test_creates_output_files(self, hits_file, tmp_path):
        outdir = tmp_path / "out"
        run_annotate(hits_path=hits_file, gff_path=GFF3_FIXTURE, outdir=outdir)
        assert (outdir / "annotated_hits.tsv").exists()
        assert (outdir / "annotation_summary.tsv").exists()
        assert (outdir / "annotate.json").exists()

    def test_annotated_hits_columns(self, hits_file, tmp_path):
        outdir = tmp_path / "out"
        run_annotate(hits_path=hits_file, gff_path=GFF3_FIXTURE, outdir=outdir)
        rows = read_tsv(outdir / "annotated_hits.tsv")
        assert len(rows) == 5
        assert set(ANNOTATED_HITS_COLUMNS).issubset(set(rows[0].keys()))

    def test_annotation_classes_correct(self, hits_file, tmp_path):
        outdir = tmp_path / "out"
        run_annotate(hits_path=hits_file, gff_path=GFF3_FIXTURE, outdir=outdir)
        rows = read_tsv(outdir / "annotated_hits.tsv")
        by_id = {r["locus_id"]: r["annotation_class"] for r in rows}
        assert by_id["PPX000001"] == "CDS"
        assert by_id["PPX000002"] == "UTR"
        assert by_id["PPX000003"] == "CDS"
        assert by_id["PPX000004"] == "intergenic"
        assert by_id["PPX000005"] == "intronic"

    def test_gene_id_populated_for_genic(self, hits_file, tmp_path):
        outdir = tmp_path / "out"
        run_annotate(hits_path=hits_file, gff_path=GFF3_FIXTURE, outdir=outdir)
        rows = read_tsv(outdir / "annotated_hits.tsv")
        by_id = {r["locus_id"]: r for r in rows}
        assert by_id["PPX000001"]["gene_id"] == "GeneA"
        assert by_id["PPX000004"]["gene_id"] == ""

    def test_summary_columns(self, hits_file, tmp_path):
        outdir = tmp_path / "out"
        run_annotate(hits_path=hits_file, gff_path=GFF3_FIXTURE, outdir=outdir)
        rows = read_tsv(outdir / "annotation_summary.tsv")
        assert set(ANNOTATION_SUMMARY_COLUMNS).issubset(set(rows[0].keys()))

    def test_summary_counts(self, hits_file, tmp_path):
        outdir = tmp_path / "out"
        run_annotate(hits_path=hits_file, gff_path=GFF3_FIXTURE, outdir=outdir)
        rows = read_tsv(outdir / "annotation_summary.tsv")
        by_cls = {r["annotation_class"]: int(r["n_loci"]) for r in rows}
        assert by_cls["CDS"] == 2
        assert by_cls["UTR"] == 1
        assert by_cls["intronic"] == 1
        assert by_cls["intergenic"] == 1
        assert by_cls.get("exonic", 0) == 0

    def test_json_metadata(self, hits_file, tmp_path):
        outdir = tmp_path / "out"
        run_annotate(hits_path=hits_file, gff_path=GFF3_FIXTURE, outdir=outdir)
        meta = json.loads((outdir / "annotate.json").read_text())
        assert meta["n_hits"] == 5
        assert "annotation_counts" in meta

    def test_creates_outdir(self, hits_file, tmp_path):
        outdir = tmp_path / "deep" / "nested" / "out"
        run_annotate(hits_path=hits_file, gff_path=GFF3_FIXTURE, outdir=outdir)
        assert outdir.is_dir()

    def test_contig_alias_hits_to_gff(self, tmp_path):
        # hits use 'Chr1', GFF3 uses 'chr1' → alias maps Chr1→chr1
        alias_path = tmp_path / "alias.tsv"
        alias_path.write_text("Chr1\tchr1\n", encoding="utf-8")
        hits = _write_hits(tmp_path / "hits.tsv", [
            _hit("PPX000001", "Chr1", 99, 100),
        ])
        outdir = tmp_path / "out"
        run_annotate(
            hits_path=hits,
            gff_path=GFF3_FIXTURE,
            outdir=outdir,
            contig_alias_path=alias_path,
            hits_contig_to_gff=True,
        )
        rows = read_tsv(outdir / "annotated_hits.tsv")
        assert rows[0]["annotation_class"] == "CDS"

    def test_empty_hits(self, tmp_path):
        empty_hits = tmp_path / "hits.tsv"
        empty_hits.write_text(_HITS_HEADER, encoding="utf-8")
        outdir = tmp_path / "out"
        run_annotate(hits_path=empty_hits, gff_path=GFF3_FIXTURE, outdir=outdir)
        rows = read_tsv(outdir / "annotated_hits.tsv")
        assert rows == []
        summary = read_tsv(outdir / "annotation_summary.tsv")
        for r in summary:
            assert int(r["n_loci"]) == 0


# ---------------------------------------------------------------------------
# Class 3: CLI integration
# ---------------------------------------------------------------------------

class TestAnnotateCli:
    def setup_method(self):
        self.runner = CliRunner()

    def test_cli_runs_successfully(self, hits_file, tmp_path):
        outdir = tmp_path / "out"
        result = self.runner.invoke(app, [
            "annotate",
            "--hits", str(hits_file),
            "--gff", str(GFF3_FIXTURE),
            "--outdir", str(outdir),
        ])
        assert result.exit_code == 0, result.output

    def test_cli_creates_outputs(self, hits_file, tmp_path):
        outdir = tmp_path / "out"
        self.runner.invoke(app, [
            "annotate",
            "--hits", str(hits_file),
            "--gff", str(GFF3_FIXTURE),
            "--outdir", str(outdir),
        ])
        assert (outdir / "annotated_hits.tsv").exists()
        assert (outdir / "annotation_summary.tsv").exists()

    def test_cli_missing_hits_exits_1(self, tmp_path):
        result = self.runner.invoke(app, [
            "annotate",
            "--hits", str(tmp_path / "no_such_hits.tsv"),
            "--gff", str(GFF3_FIXTURE),
        ])
        assert result.exit_code == 1

    def test_cli_missing_gff_exits_1(self, hits_file, tmp_path):
        result = self.runner.invoke(app, [
            "annotate",
            "--hits", str(hits_file),
            "--gff", str(tmp_path / "no_such.gff3"),
        ])
        assert result.exit_code == 1

    def test_cli_annotation_classes_correct(self, hits_file, tmp_path):
        outdir = tmp_path / "out"
        self.runner.invoke(app, [
            "annotate",
            "--hits", str(hits_file),
            "--gff", str(GFF3_FIXTURE),
            "--outdir", str(outdir),
        ])
        rows = read_tsv(outdir / "annotated_hits.tsv")
        by_id = {r["locus_id"]: r["annotation_class"] for r in rows}
        assert by_id["PPX000001"] == "CDS"
        assert by_id["PPX000004"] == "intergenic"

    def test_cli_help_available(self):
        result = self.runner.invoke(app, ["annotate", "--help"])
        assert result.exit_code == 0
        assert "annotate" in result.output.lower()

    def test_cli_contig_alias_flag(self, tmp_path):
        alias_path = tmp_path / "alias.tsv"
        alias_path.write_text("Chr1\tchr1\n", encoding="utf-8")
        hits = _write_hits(tmp_path / "hits.tsv", [
            _hit("PPX000001", "Chr1", 99, 100),
        ])
        outdir = tmp_path / "out"
        result = self.runner.invoke(app, [
            "annotate",
            "--hits", str(hits),
            "--gff", str(GFF3_FIXTURE),
            "--contig-alias", str(alias_path),
            "--hits-to-gff",
            "--outdir", str(outdir),
        ])
        assert result.exit_code == 0
        rows = read_tsv(outdir / "annotated_hits.tsv")
        assert rows[0]["annotation_class"] == "CDS"
