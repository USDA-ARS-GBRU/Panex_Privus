"""Unit tests for hierarchical base assignment + run-length smoothing."""

from __future__ import annotations

import pytest

from privy.core.intervals import (
    class_proportions,
    hierarchical_base_assignment,
    run_length_smooth,
)


class TestHierarchicalBaseAssignment:
    def test_priority_resolves_overlap(self):
        feats = [(0, 50, "A"), (40, 80, "B")]
        counts = hierarchical_base_assignment(0, 100, feats, ["A", "B"])
        assert counts == {"A": 50, "B": 30, "missing": 20}
        assert sum(counts.values()) == 100

    def test_priority_order_matters(self):
        feats = [(0, 50, "A"), (40, 80, "B")]
        counts = hierarchical_base_assignment(0, 100, feats, ["B", "A"])
        assert counts == {"B": 40, "A": 40, "missing": 20}

    def test_no_features_all_missing(self):
        counts = hierarchical_base_assignment(0, 100, [], ["A"])
        assert counts == {"A": 0, "missing": 100}

    def test_unknown_class_ignored(self):
        feats = [(0, 30, "A"), (0, 100, "ignore_me")]
        counts = hierarchical_base_assignment(0, 100, feats, ["A"])
        assert counts == {"A": 30, "missing": 70}

    def test_features_clipped_to_window(self):
        feats = [(-50, 200, "A")]
        counts = hierarchical_base_assignment(10, 60, feats, ["A"])
        assert counts == {"A": 50, "missing": 0}


class TestClassProportions:
    def test_proportions_sum_to_one(self):
        props = class_proportions({"A": 50, "B": 30, "missing": 20})
        assert props == {"A": 0.5, "B": 0.3, "missing": 0.2}
        assert sum(props.values()) == pytest.approx(1.0)

    def test_empty_is_zero(self):
        assert class_proportions({"A": 0, "missing": 0}) == {"A": 0.0, "missing": 0.0}


class TestRunLengthSmooth:
    def test_short_interior_run_adopts_previous(self):
        assert run_length_smooth(["A", "A", "B", "A", "A"], 2) == ["A", "A", "A", "A", "A"]

    def test_leading_short_run_adopts_next(self):
        assert run_length_smooth(["B", "A", "A", "A"], 2) == ["A", "A", "A", "A"]

    def test_min_run_one_is_noop(self):
        assert run_length_smooth(["A", "B", "A"], 1) == ["A", "B", "A"]

    def test_all_short_unchanged(self):
        assert run_length_smooth(["A", "B", "C"], 2) == ["A", "B", "C"]

    def test_keeps_long_runs(self):
        assert run_length_smooth(["A", "A", "A", "B", "B", "B"], 2) == [
            "A", "A", "A", "B", "B", "B",
        ]
