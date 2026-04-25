"""Unit tests for src/privy/io/gfa.py.

Tests the GFA parser, inverted-index building, position queries,
and cohort-count extraction — all pure Python, no external tools required.

Fixture layout (small_cohort.gfa):
    Samples: T1, T2 (targets), O1, O2, O3 (off-targets)
    Bubble 1 (chr1:8-18): T1+T2 → s2_target, O1+O2+O3 → s2_offt
    Bubble 2 (chr1:60-67): T1 → s4_target, T2 MISSING, O1+O2+O3 → s4_offt
    Backbone:  s1 (0-8), s3 (18-26), s5 (80-88) — all five samples traverse
"""

from __future__ import annotations

from pathlib import Path

import pytest

from privy.io.gfa import (
    GfaGraph,
    extract_cohort_segment_counts,
    get_gfa_samples,
    get_samples_present_at_locus,
    get_samples_traversing_segment,
    parse_gfa,
    query_segments_at_locus,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

GFA_DATA = Path(__file__).parent.parent / "data" / "small_cohort.gfa"


@pytest.fixture(scope="module")
def graph() -> GfaGraph:
    return parse_gfa(GFA_DATA)


@pytest.fixture(scope="module")
def targets() -> list[str]:
    return ["T1", "T2"]


@pytest.fixture(scope="module")
def offtargets() -> list[str]:
    return ["O1", "O2", "O3"]


# ---------------------------------------------------------------------------
# TestParsing — basic file parsing
# ---------------------------------------------------------------------------


class TestParsing:
    def test_segment_count(self, graph: GfaGraph) -> None:
        assert len(graph.segments) == 7

    def test_link_count(self, graph: GfaGraph) -> None:
        assert len(graph.links) == 8

    def test_walk_count(self, graph: GfaGraph) -> None:
        # 5 walks for bubble 1 + 4 walks for bubble 2 (T2 absent)
        assert len(graph.walks) == 9

    def test_no_paths(self, graph: GfaGraph) -> None:
        # fixture uses W-lines only
        assert len(graph.paths) == 0

    def test_segment_names(self, graph: GfaGraph) -> None:
        expected = {"s1", "s2_target", "s2_offt", "s3", "s4_target", "s4_offt", "s5"}
        assert set(graph.segments) == expected

    def test_segment_length(self, graph: GfaGraph) -> None:
        assert graph.segments["s2_target"].length == 10
        assert graph.segments["s4_target"].length == 7

    def test_segment_coordinates(self, graph: GfaGraph) -> None:
        seg = graph.segments["s2_target"]
        assert seg.ref_contig == "chr1"
        assert seg.ref_start == 8
        assert seg.ref_end == 18

    def test_header_version(self, graph: GfaGraph) -> None:
        assert graph.header_tags.get("VN") == "1.1"

    def test_walk_sample_names(self, graph: GfaGraph) -> None:
        samples_in_walks = {w.sample for w in graph.walks}
        assert samples_in_walks == {"T1", "T2", "O1", "O2", "O3"}

    def test_walk_coordinates(self, graph: GfaGraph) -> None:
        t1_walks = [w for w in graph.walks if w.sample == "T1"]
        assert any(w.seq_start == 0 and w.seq_end == 26 for w in t1_walks)
        assert any(w.seq_start == 57 and w.seq_end == 88 for w in t1_walks)

    def test_walk_steps_parsed(self, graph: GfaGraph) -> None:
        # T1 bubble-1 walk: >s1>s2_target>s3
        t1_b1 = next(
            w for w in graph.walks if w.sample == "T1" and w.seq_end == 26
        )
        assert len(t1_b1.steps) == 3
        assert t1_b1.steps[0].segment == "s1"
        assert t1_b1.steps[0].orient == "+"
        assert t1_b1.steps[1].segment == "s2_target"
        assert t1_b1.steps[2].segment == "s3"

    def test_t2_absent_from_bubble2(self, graph: GfaGraph) -> None:
        # T2 has no walk with seq_id==chr1 and seq_start >= 57
        t2_walks = [w for w in graph.walks if w.sample == "T2"]
        assert all(w.seq_end <= 26 for w in t2_walks)


# ---------------------------------------------------------------------------
# TestParsePLine — P-line path fixture
# ---------------------------------------------------------------------------


@pytest.fixture
def pline_gfa(tmp_path: Path) -> GfaGraph:
    """Minimal GFA1 fixture with only P-lines (no W-lines)."""
    content = (
        "H\tVN:Z:1.0\n"
        "S\tsegA\t*\tLN:i:5\tSN:Z:chr1\tSO:i:0\n"
        "S\tsegB\t*\tLN:i:5\tSN:Z:chr1\tSO:i:5\n"
        "L\tsegA\t+\tsegB\t+\t0M\n"
        "P\tHG001#1#chr1\tsegA+,segB+\t*\n"
        "P\tHG002#1#chr1\tsegA+,segB+\t*\n"
    )
    p = tmp_path / "test.gfa"
    p.write_text(content)
    return parse_gfa(p)


class TestParsePLine:
    def test_path_sample_extraction(self, pline_gfa: GfaGraph) -> None:
        assert "HG001" in pline_gfa.sample_to_paths
        assert "HG002" in pline_gfa.sample_to_paths

    def test_segment_to_paths_index(self, pline_gfa: GfaGraph) -> None:
        assert "HG001#1#chr1" in pline_gfa.segment_to_paths["segA"]
        assert "HG002#1#chr1" in pline_gfa.segment_to_paths["segA"]

    def test_path_segment_list(self, pline_gfa: GfaGraph) -> None:
        path = pline_gfa.paths["HG001#1#chr1"]
        assert path.segment_names == ["segA", "segB"]
        assert path.orientations == ["+", "+"]

    def test_path_haplotype_extraction(self, pline_gfa: GfaGraph) -> None:
        path = pline_gfa.paths["HG001#1#chr1"]
        assert path.sample == "HG001"
        assert path.haplotype == 1

    def test_plain_path_name(self, tmp_path: Path) -> None:
        content = (
            "H\tVN:Z:1.0\n"
            "S\tX\t*\tLN:i:5\tSN:Z:chr1\tSO:i:0\n"
            "P\tMySample\tX+\t*\n"
        )
        p = tmp_path / "plain.gfa"
        p.write_text(content)
        g = parse_gfa(p)
        assert "MySample" in g.sample_to_paths


# ---------------------------------------------------------------------------
# TestInvertedIndices
# ---------------------------------------------------------------------------


class TestInvertedIndices:
    def test_sample_to_walks_populated(self, graph: GfaGraph) -> None:
        for sample in ["T1", "T2", "O1", "O2", "O3"]:
            assert sample in graph.sample_to_walks

    def test_segment_to_walks_s2_target(self, graph: GfaGraph) -> None:
        # T1 and T2 traverse s2_target; O* do not
        path_samples = {
            graph.walks[i].sample for i in graph.segment_to_walks.get("s2_target", [])
        }
        assert path_samples == {"T1", "T2"}

    def test_segment_to_walks_s2_offt(self, graph: GfaGraph) -> None:
        path_samples = {
            graph.walks[i].sample for i in graph.segment_to_walks.get("s2_offt", [])
        }
        assert path_samples == {"O1", "O2", "O3"}

    def test_segment_to_walks_backbone(self, graph: GfaGraph) -> None:
        # s1 is traversed by all five samples in bubble-1 walks
        s1_samples = {
            graph.walks[i].sample for i in graph.segment_to_walks.get("s1", [])
        }
        assert s1_samples == {"T1", "T2", "O1", "O2", "O3"}

    def test_position_index_built(self, graph: GfaGraph) -> None:
        assert "chr1" in graph._contig_segments
        assert len(graph._contig_segments["chr1"]) == 7

    def test_position_index_sorted(self, graph: GfaGraph) -> None:
        entries = graph._contig_segments["chr1"]
        starts = [e[0] for e in entries]
        assert starts == sorted(starts)


# ---------------------------------------------------------------------------
# TestGetGfaSamples
# ---------------------------------------------------------------------------


class TestGetGfaSamples:
    def test_returns_all_five_samples(self) -> None:
        samples = get_gfa_samples(GFA_DATA)
        assert samples == ["O1", "O2", "O3", "T1", "T2"]

    def test_returns_sorted_list(self) -> None:
        samples = get_gfa_samples(GFA_DATA)
        assert samples == sorted(samples)

    def test_file_not_found(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            get_gfa_samples(tmp_path / "nonexistent.gfa")


# ---------------------------------------------------------------------------
# TestQuerySegmentsAtLocus
# ---------------------------------------------------------------------------


class TestQuerySegmentsAtLocus:
    def test_bubble1_segments(self, graph: GfaGraph) -> None:
        segs = query_segments_at_locus(graph, "chr1", 8, 18)
        names = {s.name for s in segs}
        assert names == {"s2_target", "s2_offt"}

    def test_backbone_s1(self, graph: GfaGraph) -> None:
        segs = query_segments_at_locus(graph, "chr1", 0, 8)
        assert {s.name for s in segs} == {"s1"}

    def test_no_segments_past_end(self, graph: GfaGraph) -> None:
        segs = query_segments_at_locus(graph, "chr1", 200, 300)
        assert segs == []

    def test_unknown_contig(self, graph: GfaGraph) -> None:
        segs = query_segments_at_locus(graph, "chrX", 0, 1000)
        assert segs == []

    def test_partial_overlap(self, graph: GfaGraph) -> None:
        # Query starting at 5 overlaps s1 (0-8) and s2_target (8-18)
        segs = query_segments_at_locus(graph, "chr1", 5, 12)
        names = {s.name for s in segs}
        assert "s1" in names
        assert "s2_target" in names

    def test_bubble2_segments(self, graph: GfaGraph) -> None:
        segs = query_segments_at_locus(graph, "chr1", 60, 67)
        names = {s.name for s in segs}
        assert names == {"s4_target", "s4_offt"}


# ---------------------------------------------------------------------------
# TestGetSamplesTraversingSegment
# ---------------------------------------------------------------------------


class TestGetSamplesTraversingSegment:
    def test_s2_target_traversed_by_targets(self, graph: GfaGraph) -> None:
        assert get_samples_traversing_segment(graph, "s2_target") == {"T1", "T2"}

    def test_s2_offt_traversed_by_offtargets(self, graph: GfaGraph) -> None:
        assert get_samples_traversing_segment(graph, "s2_offt") == {"O1", "O2", "O3"}

    def test_s4_target_t1_only(self, graph: GfaGraph) -> None:
        assert get_samples_traversing_segment(graph, "s4_target") == {"T1"}

    def test_backbone_s1_all_samples(self, graph: GfaGraph) -> None:
        assert get_samples_traversing_segment(graph, "s1") == {"T1", "T2", "O1", "O2", "O3"}

    def test_unknown_segment_returns_empty(self, graph: GfaGraph) -> None:
        assert get_samples_traversing_segment(graph, "does_not_exist") == set()


# ---------------------------------------------------------------------------
# TestGetSamplesPresentAtLocus
# ---------------------------------------------------------------------------


class TestGetSamplesPresentAtLocus:
    def test_bubble1_all_present(self, graph: GfaGraph) -> None:
        # All five samples have walks covering chr1:8-18 (via chr1:0-26 walks)
        present = get_samples_present_at_locus(graph, "chr1", 8, 18)
        assert present == {"T1", "T2", "O1", "O2", "O3"}

    def test_bubble2_t2_missing(self, graph: GfaGraph) -> None:
        # T2 has no walk covering chr1:60-67
        present = get_samples_present_at_locus(graph, "chr1", 60, 67)
        assert present == {"T1", "O1", "O2", "O3"}
        assert "T2" not in present

    def test_gap_region_empty(self, graph: GfaGraph) -> None:
        # No walks or coordinate segments at chr1:200-300
        present = get_samples_present_at_locus(graph, "chr1", 200, 300)
        assert present == set()

    def test_unknown_contig_empty(self, graph: GfaGraph) -> None:
        present = get_samples_present_at_locus(graph, "chrX", 0, 100)
        assert present == set()


# ---------------------------------------------------------------------------
# TestExtractCohortSegmentCounts
# ---------------------------------------------------------------------------


class TestExtractCohortSegmentCounts:
    def test_s2_target_strict_complete(
        self, graph: GfaGraph, targets: list[str], offtargets: list[str]
    ) -> None:
        """s2_target: all targets traverse, all off-targets present but absent → strict_complete."""
        seg = graph.segments["s2_target"]
        ts, tt, os_, ot, tm, om = extract_cohort_segment_counts(
            graph, "s2_target",
            seg.ref_contig, seg.ref_start, seg.ref_end,
            targets, offtargets,
        )
        assert ts == 2   # T1+T2 traverse
        assert tt == 2
        assert os_ == 0  # O1+O2+O3 present but traverse s2_offt
        assert ot == 3
        assert tm == 0   # no target missing
        assert om == 0   # no off-target missing

    def test_s4_target_strict_target_missing(
        self, graph: GfaGraph, targets: list[str], offtargets: list[str]
    ) -> None:
        """s4_target: T1 traverses, T2 missing, O* present but absent."""
        seg = graph.segments["s4_target"]
        ts, tt, os_, ot, tm, om = extract_cohort_segment_counts(
            graph, "s4_target",
            seg.ref_contig, seg.ref_start, seg.ref_end,
            targets, offtargets,
        )
        assert ts == 1   # T1 only
        assert tt == 2
        assert os_ == 0  # O* traverse s4_offt
        assert ot == 3
        assert tm == 1   # T2 missing (no walk at 60-67)
        assert om == 0

    def test_s2_offt_contradicted(
        self, graph: GfaGraph, targets: list[str], offtargets: list[str]
    ) -> None:
        """s2_offt: off-targets traverse it → should be contradicted."""
        seg = graph.segments["s2_offt"]
        ts, tt, os_, ot, tm, om = extract_cohort_segment_counts(
            graph, "s2_offt",
            seg.ref_contig, seg.ref_start, seg.ref_end,
            targets, offtargets,
        )
        assert ts == 0
        assert os_ == 3   # all off-targets traverse s2_offt

    def test_backbone_s1_all_traverse(
        self, graph: GfaGraph, targets: list[str], offtargets: list[str]
    ) -> None:
        seg = graph.segments["s1"]
        ts, tt, os_, ot, tm, om = extract_cohort_segment_counts(
            graph, "s1",
            seg.ref_contig, seg.ref_start, seg.ref_end,
            targets, offtargets,
        )
        assert ts == 2
        assert os_ == 3

    def test_ghost_sample_counted_as_missing(
        self, graph: GfaGraph, offtargets: list[str]
    ) -> None:
        """A target sample name not present in the GFA at all counts as missing."""
        seg = graph.segments["s2_target"]
        ts, tt, os_, ot, tm, om = extract_cohort_segment_counts(
            graph, "s2_target",
            seg.ref_contig, seg.ref_start, seg.ref_end,
            ["T1", "GHOST"],   # GHOST not in GFA
            offtargets,
        )
        assert tt == 2
        assert tm == 1   # GHOST → missing

    def test_no_coord_info_no_missing_detection(
        self, graph: GfaGraph, targets: list[str], offtargets: list[str]
    ) -> None:
        """When coord info is omitted, missing detection is disabled."""
        ts, tt, os_, ot, tm, om = extract_cohort_segment_counts(
            graph, "s4_target",
            None, None, None,   # no coordinate info
            targets, offtargets,
        )
        assert ts == 1
        assert tm == 0   # T2 not detected as missing without coords


# ---------------------------------------------------------------------------
# TestParseErrors
# ---------------------------------------------------------------------------


class TestParseErrors:
    def test_gfa2_raises(self, tmp_path: Path) -> None:
        content = "H\tVN:Z:2.0\n"
        p = tmp_path / "gfa2.gfa"
        p.write_text(content)
        with pytest.raises(ValueError, match="GFA version 2 is not supported"):
            parse_gfa(p)

    def test_file_not_found(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            parse_gfa(tmp_path / "missing.gfa")

    def test_s_line_too_few_fields(self, tmp_path: Path) -> None:
        content = "H\tVN:Z:1.0\nS\tname\n"
        p = tmp_path / "bad.gfa"
        p.write_text(content)
        with pytest.raises(ValueError, match="S-line"):
            parse_gfa(p)

    def test_w_line_too_few_fields(self, tmp_path: Path) -> None:
        content = "H\tVN:Z:1.1\nW\tS1\t1\tchr1\n"
        p = tmp_path / "bad_w.gfa"
        p.write_text(content)
        with pytest.raises(ValueError, match="W-line"):
            parse_gfa(p)

    def test_comment_lines_skipped(self, tmp_path: Path) -> None:
        content = "H\tVN:Z:1.0\n# this is a comment\nS\tsegX\t*\tLN:i:5\n"
        p = tmp_path / "comments.gfa"
        p.write_text(content)
        g = parse_gfa(p)
        assert "segX" in g.segments

    def test_empty_file_parses_cleanly(self, tmp_path: Path) -> None:
        p = tmp_path / "empty.gfa"
        p.write_text("")
        g = parse_gfa(p)
        assert len(g.segments) == 0
        assert len(g.walks) == 0
