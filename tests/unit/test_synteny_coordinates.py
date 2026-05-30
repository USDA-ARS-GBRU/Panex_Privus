"""Unit tests for src/privy/synteny/coordinates.py."""

from __future__ import annotations

from pathlib import Path

import pytest

from privy.io.gfa import parse_gfa
from privy.synteny.coordinates import (
    PathCoordinateModel,
    PathStepLocation,
    SegmentOccurrence,
)

# ---------------------------------------------------------------------------
# Fixtures: a tiny GFA with P-lines and W-lines sharing segments.
# Segment lengths: s1=10, s2=5, s3=20, s4=5
#   pA (P-line "sampleA#0#chr1"): s1+ s2+ s3+      -> starts [0,10,15], len 35
#   pC (P-line "sampleC#0#chr1"): s1+ s2+ s1-      -> s1 occurs at steps 0 and 2
#   W-line sampleB 0 chr1 100 135: >s1>s2>s3       -> stable base 100
# ---------------------------------------------------------------------------

_GFA_LINES = [
    "H\tVN:Z:1.1",
    "S\ts1\t*\tLN:i:10",
    "S\ts2\t*\tLN:i:5",
    "S\ts3\t*\tLN:i:20",
    "S\ts4\t*\tLN:i:5",
    "P\tsampleA#0#chr1\ts1+,s2+,s3+\t*",
    "P\tsampleC#0#chr1\ts1+,s2+,s1-\t*",
    "W\tsampleB\t0\tchr1\t100\t135\t>s1>s2>s3",
]


def _model(tmp_path: Path) -> PathCoordinateModel:
    p = tmp_path / "g.gfa"
    p.write_text("\n".join(_GFA_LINES) + "\n", encoding="utf-8")
    return PathCoordinateModel.from_graph(parse_gfa(p))


# ---------------------------------------------------------------------------
# Class 1: construction & container
# ---------------------------------------------------------------------------


class TestConstruction:
    def test_indexes_paths_and_walks(self, tmp_path):
        m = _model(tmp_path)
        ids = m.path_ids()
        assert "sampleA#0#chr1" in ids
        assert "sampleC#0#chr1" in ids
        # walk id built from sample#hap#seqid
        assert "sampleB#0#chr1" in ids
        assert len(m) == 3
        assert "sampleA#0#chr1" in m
        assert "nope" not in m

    def test_path_lengths(self, tmp_path):
        m = _model(tmp_path)
        assert m.path_length("sampleA#0#chr1") == 35   # 10+5+20
        assert m.path_length("sampleB#0#chr1") == 35

    def test_unknown_segment_strict_raises(self, tmp_path):
        lines = ["H\tVN:Z:1.1", "S\ts1\t*\tLN:i:10", "P\tp#0#c\ts1+,sX+\t*"]
        p = tmp_path / "bad.gfa"
        p.write_text("\n".join(lines) + "\n", encoding="utf-8")
        with pytest.raises(KeyError, match="unknown segment"):
            PathCoordinateModel.from_graph(parse_gfa(p))

    def test_unknown_segment_nonstrict_skips(self, tmp_path):
        lines = [
            "H\tVN:Z:1.1",
            "S\ts1\t*\tLN:i:10",
            "P\tgood#0#c\ts1+\t*",
            "P\tbad#0#c\ts1+,sX+\t*",
        ]
        p = tmp_path / "mixed.gfa"
        p.write_text("\n".join(lines) + "\n", encoding="utf-8")
        m = PathCoordinateModel.from_graph(parse_gfa(p), strict=False)
        assert "good#0#c" in m
        assert "bad#0#c" not in m


# ---------------------------------------------------------------------------
# Class 2: locate
# ---------------------------------------------------------------------------


class TestLocate:
    def test_locate_first_segment(self, tmp_path):
        m = _model(tmp_path)
        loc = m.locate("sampleA#0#chr1", 3)
        assert loc == PathStepLocation(
            path_id="sampleA#0#chr1",
            path_position=3,
            segment="s1",
            orientation="+",
            step_index=0,
            offset_in_segment=3,
        )

    def test_locate_boundary_into_second_segment(self, tmp_path):
        m = _model(tmp_path)
        loc = m.locate("sampleA#0#chr1", 10)   # first base of s2
        assert loc.segment == "s2"
        assert loc.step_index == 1
        assert loc.offset_in_segment == 0

    def test_locate_within_third_segment(self, tmp_path):
        m = _model(tmp_path)
        loc = m.locate("sampleA#0#chr1", 20)   # 20 -> s3 starts at 15, offset 5
        assert loc.segment == "s3"
        assert loc.offset_in_segment == 5

    def test_locate_last_base(self, tmp_path):
        m = _model(tmp_path)
        loc = m.locate("sampleA#0#chr1", 34)
        assert loc.segment == "s3"
        assert loc.offset_in_segment == 19

    @pytest.mark.parametrize("pos", [-1, 35, 999])
    def test_locate_out_of_range(self, tmp_path, pos):
        m = _model(tmp_path)
        with pytest.raises(IndexError):
            m.locate("sampleA#0#chr1", pos)

    def test_unknown_path_raises(self, tmp_path):
        m = _model(tmp_path)
        with pytest.raises(KeyError):
            m.locate("ghost", 0)


# ---------------------------------------------------------------------------
# Class 3: occurrences (CNV / repeated segment)
# ---------------------------------------------------------------------------


class TestOccurrences:
    def test_single_occurrence(self, tmp_path):
        m = _model(tmp_path)
        occ = m.occurrences("sampleA#0#chr1", "s2")
        assert occ == [SegmentOccurrence(step_index=1, start=10, end=15, orientation="+")]

    def test_repeated_segment_multiple_occurrences(self, tmp_path):
        m = _model(tmp_path)
        occ = m.occurrences("sampleC#0#chr1", "s1")
        # s1 at step 0 (0-10, +) and step 2 (15-25, -)
        assert occ == [
            SegmentOccurrence(step_index=0, start=0, end=10, orientation="+"),
            SegmentOccurrence(step_index=2, start=15, end=25, orientation="-"),
        ]

    def test_absent_segment_returns_empty(self, tmp_path):
        m = _model(tmp_path)
        assert m.occurrences("sampleA#0#chr1", "s4") == []


# ---------------------------------------------------------------------------
# Class 4: stable coordinate projection
# ---------------------------------------------------------------------------


class TestToStable:
    def test_pline_base_offset_zero(self, tmp_path):
        m = _model(tmp_path)
        assert m.to_stable("sampleA#0#chr1", 12) == ("chr1", 12)
        assert m.stable_contig("sampleA#0#chr1") == "chr1"

    def test_wline_uses_seq_start_offset(self, tmp_path):
        m = _model(tmp_path)
        # W-line seq_start=100, so path position 12 -> chr1:112
        assert m.to_stable("sampleB#0#chr1", 12) == ("chr1", 112)

    def test_to_stable_out_of_range(self, tmp_path):
        m = _model(tmp_path)
        with pytest.raises(IndexError):
            m.to_stable("sampleA#0#chr1", 35)


# ---------------------------------------------------------------------------
# Class 5: cross-path primitive (foundation for P1 projection)
# ---------------------------------------------------------------------------


class TestCrossPathPrimitive:
    def test_same_segment_maps_between_paths(self, tmp_path):
        """A shared segment lets us bridge two coordinate systems (P1 will build on this)."""
        m = _model(tmp_path)
        # Position 20 on pA -> s3 offset 5. Find s3 on the W-line path and map back.
        loc = m.locate("sampleA#0#chr1", 20)
        assert loc.segment == "s3"
        (occ,) = m.occurrences("sampleB#0#chr1", loc.segment)
        # same orientation here, so path-local position = occ.start + offset
        b_position = occ.start + loc.offset_in_segment
        assert m.to_stable("sampleB#0#chr1", b_position) == ("chr1", 100 + 20)
