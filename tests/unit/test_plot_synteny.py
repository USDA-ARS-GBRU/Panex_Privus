"""Unit tests for src/privy/plot/synteny.py."""

from __future__ import annotations

from pathlib import Path

from privy.io.gfa import parse_gfa
from privy.plot.synteny import plot_dotplot, plot_riparian
from privy.synteny.build import build_synteny
from privy.synteny.coordinates import PathCoordinateModel
from privy.synthetic import inversion_pangenome


def _block_rows(tmp_path: Path) -> list[dict]:
    model = PathCoordinateModel.from_graph(parse_gfa(
        inversion_pangenome(seg_len=10).write(tmp_path / "g.gfa")
    ))
    result = build_synteny(model, "sample0#0#chr1")
    return [
        {
            "block_id": b.block_id,
            "query_genome": b.query.genome,
            "query_start": b.query.start,
            "query_end": b.query.end,
            "ref_genome": b.target.genome,
            "ref_contig": b.target.contig,
            "ref_start": b.target.start,
            "ref_end": b.target.end,
            "strand": b.strand,
            "block_type": b.block_type.value,
        }
        for b in result.blocks
    ]


class TestRiparian:
    def test_writes_nonempty_png(self, tmp_path):
        rows = _block_rows(tmp_path)
        out = plot_riparian(rows, tmp_path / "fig")
        assert out.exists()
        assert out.suffix == ".png"
        assert out.stat().st_size > 0

    def test_pdf_format(self, tmp_path):
        rows = _block_rows(tmp_path)
        out = plot_riparian(rows, tmp_path / "fig", output_format="pdf")
        assert out.suffix == ".pdf"
        assert out.stat().st_size > 0

    def test_empty_rows_still_writes_placeholder(self, tmp_path):
        out = plot_riparian([], tmp_path / "fig")
        assert out.exists()
        assert out.stat().st_size > 0


class TestDotplot:
    def test_writes_nonempty(self, tmp_path):
        rows = _block_rows(tmp_path)
        out = plot_dotplot(rows, tmp_path / "fig")
        assert out.exists()
        assert out.stat().st_size > 0

    def test_empty_rows_placeholder(self, tmp_path):
        out = plot_dotplot([], tmp_path / "fig")
        assert out.exists()


class TestRowCoercion:
    def test_string_numeric_columns_accepted(self, tmp_path):
        # rows as they'd come from a TSV (all strings)
        rows = [{
            "block_id": "B0", "query_genome": "qA", "query_start": "0", "query_end": "100",
            "ref_genome": "rB", "ref_contig": "chr1", "ref_start": "0", "ref_end": "100",
            "strand": "+", "block_type": "collinear",
        }]
        out = plot_riparian(rows, tmp_path / "fig")
        assert out.stat().st_size > 0
