"""Unit tests for src/privy/io/paf.py."""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from privy.io.paf import (
    BedpeRecord,
    PafRecord,
    PafTag,
    format_paf_record,
    parse_bedpe,
    parse_paf,
    write_paf,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

# Two minimap2-style PAF lines: a forward and a reverse alignment, with tags.
_PAF_LINES = [
    "q1\t1000\t100\t600\t+\tchrA\t5000\t1100\t1600\t480\t500\t60\t"
    "tp:A:P\tcg:Z:500M\tNM:i:20\tdv:f:0.04",
    "q2\t800\t0\t300\t-\tchrB\t4000\t200\t500\t290\t300\t30\tcs:Z::300",
]


def _write_paf(tmp_path: Path, lines: list[str], gz: bool = False) -> Path:
    content = "\n".join(lines) + "\n"
    if gz:
        p = tmp_path / "test.paf.gz"
        with gzip.open(p, "wt", encoding="utf-8") as fh:
            fh.write(content)
    else:
        p = tmp_path / "test.paf"
        p.write_text(content, encoding="utf-8")
    return p


# ---------------------------------------------------------------------------
# Class 1: parse_paf
# ---------------------------------------------------------------------------


class TestParsePaf:
    def test_parses_mandatory_columns(self, tmp_path):
        p = _write_paf(tmp_path, _PAF_LINES)
        records = list(parse_paf(p))
        assert len(records) == 2
        r0 = records[0]
        assert r0.query_name == "q1"
        assert r0.query_length == 1000
        assert r0.query_start == 100
        assert r0.query_end == 600
        assert r0.strand == "+"
        assert r0.target_name == "chrA"
        assert r0.target_length == 5000
        assert r0.target_start == 1100
        assert r0.target_end == 1600
        assert r0.residue_matches == 480
        assert r0.alignment_block_length == 500
        assert r0.mapping_quality == 60

    def test_parses_typed_tags(self, tmp_path):
        p = _write_paf(tmp_path, _PAF_LINES)
        r0 = next(parse_paf(p))
        assert r0.get_tag("tp") == "P"           # A -> str
        assert r0.cigar == "500M"                # cg:Z
        assert r0.get_tag("NM") == 20            # i -> int
        assert isinstance(r0.get_tag("NM"), int)
        assert r0.get_tag("dv") == pytest.approx(0.04)  # f -> float
        assert isinstance(r0.get_tag("dv"), float)
        assert r0.tags["cg"] == PafTag("cg", "Z", "500M")

    def test_cs_tag_and_missing_tag_default(self, tmp_path):
        p = _write_paf(tmp_path, _PAF_LINES)
        records = list(parse_paf(p))
        assert records[1].cs == ":300"
        assert records[1].cigar is None
        assert records[0].cs is None
        assert records[0].get_tag("absent", "fallback") == "fallback"

    def test_reads_gzipped(self, tmp_path):
        p = _write_paf(tmp_path, _PAF_LINES, gz=True)
        records = list(parse_paf(p))
        assert len(records) == 2
        assert records[1].query_name == "q2"

    def test_skips_blank_and_comment_lines(self, tmp_path):
        p = _write_paf(tmp_path, ["# header", "", _PAF_LINES[0], "  ", _PAF_LINES[1]])
        # "  " is non-empty after split -> would be malformed; use skip_malformed.
        records = list(parse_paf(p, skip_malformed=True))
        assert len(records) == 2

    def test_raises_on_too_few_columns(self, tmp_path):
        p = _write_paf(tmp_path, ["q1\t1000\t100\t600\t+\tchrA"])
        with pytest.raises(ValueError, match="at least 12"):
            list(parse_paf(p))

    def test_raises_on_bad_strand(self, tmp_path):
        bad = "q1\t1000\t100\t600\t*\tchrA\t5000\t1100\t1600\t480\t500\t60"
        p = _write_paf(tmp_path, [bad])
        with pytest.raises(ValueError, match="strand"):
            list(parse_paf(p))

    def test_raises_on_nonnumeric_field(self, tmp_path):
        bad = "q1\tNOTINT\t100\t600\t+\tchrA\t5000\t1100\t1600\t480\t500\t60"
        p = _write_paf(tmp_path, [bad])
        with pytest.raises(ValueError):
            list(parse_paf(p))

    def test_skip_malformed_continues(self, tmp_path):
        p = _write_paf(tmp_path, ["truncated\tline", _PAF_LINES[0]])
        records = list(parse_paf(p, skip_malformed=True))
        assert len(records) == 1
        assert records[0].query_name == "q1"


# ---------------------------------------------------------------------------
# Class 2: derived accessors
# ---------------------------------------------------------------------------


class TestPafRecordAccessors:
    def test_identity_and_lengths(self):
        r = PafRecord("q", 1000, 100, 600, "+", "t", 5000, 1100, 1600, 480, 500, 60)
        assert r.blast_identity == pytest.approx(0.96)
        assert r.query_aligned_length == 500
        assert r.target_aligned_length == 500
        assert r.is_reverse is False

    def test_reverse_strand(self):
        r = PafRecord("q", 800, 0, 300, "-", "t", 4000, 200, 500, 290, 300, 30)
        assert r.is_reverse is True

    def test_identity_zero_block_length(self):
        r = PafRecord("q", 10, 0, 0, "+", "t", 10, 0, 0, 0, 0, 0)
        assert r.blast_identity == 0.0


# ---------------------------------------------------------------------------
# Class 3: write / round-trip
# ---------------------------------------------------------------------------


class TestWritePaf:
    def test_format_record_roundtrip(self):
        original = (
            "q1\t1000\t100\t600\t+\tchrA\t5000\t1100\t1600\t480\t500\t60\t"
            "tp:A:P\tcg:Z:500M\tNM:i:20"
        )
        record = next(iter([_only(original)]))
        assert format_paf_record(record) == original

    def test_write_and_reparse_is_stable(self, tmp_path):
        records = list(parse_paf(_write_paf(tmp_path, _PAF_LINES)))
        out = tmp_path / "out.paf"
        n = write_paf(records, out)
        assert n == 2
        reparsed = list(parse_paf(out))
        assert len(reparsed) == 2
        assert reparsed[0].query_name == records[0].query_name
        assert reparsed[0].get_tag("NM") == 20
        assert reparsed[0].cigar == "500M"
        assert reparsed[1].cs == ":300"

    def test_write_gzip(self, tmp_path):
        records = list(parse_paf(_write_paf(tmp_path, _PAF_LINES)))
        out = tmp_path / "out.paf.gz"
        write_paf(records, out)
        reparsed = list(parse_paf(out))
        assert len(reparsed) == 2


def _only(line: str) -> PafRecord:
    from privy.io.paf import _parse_paf_line

    return _parse_paf_line(line, 1)


# ---------------------------------------------------------------------------
# Class 4: BEDPE
# ---------------------------------------------------------------------------


class TestParseBedpe:
    def test_parses_full_record(self, tmp_path):
        line = "chr1\t100\t200\tchr2\t1100\t1200\tblockA\t60\t+\t-\textra1"
        p = tmp_path / "t.bedpe"
        p.write_text(line + "\n", encoding="utf-8")
        rec = next(parse_bedpe(p))
        assert rec == BedpeRecord(
            "chr1", 100, 200, "chr2", 1100, 1200, "blockA", "60", "+", "-", ["extra1"]
        )

    def test_minimal_six_columns(self, tmp_path):
        p = tmp_path / "t.bedpe"
        p.write_text("chr1\t0\t50\tchr1\t900\t950\n", encoding="utf-8")
        rec = next(parse_bedpe(p))
        assert rec.chrom1 == "chr1"
        assert rec.end2 == 950
        assert rec.name is None
        assert rec.extra == []

    def test_raises_on_too_few_columns(self, tmp_path):
        p = tmp_path / "t.bedpe"
        p.write_text("chr1\t0\t50\n", encoding="utf-8")
        with pytest.raises(ValueError, match="at least 6"):
            list(parse_bedpe(p))

    def test_skips_track_and_comment_lines(self, tmp_path):
        p = tmp_path / "t.bedpe"
        p.write_text(
            "track name=x\n# c\nchr1\t0\t50\tchr1\t900\t950\n", encoding="utf-8"
        )
        assert len(list(parse_bedpe(p))) == 1
