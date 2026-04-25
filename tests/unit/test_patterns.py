"""Unit tests for privy.core.patterns — the core classification logic.

These tests cover all six strictness classes and boundary conditions.
They are the most critical unit tests in the package: the classification
logic is the conceptual center of Panex Privus.
"""

import pytest

from privy.core.patterns import (
    StrictnessClass,
    build_allele_pattern,
    classify_strictness,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def classify(
    ts: int, tt: int,
    os: int, ot: int,
    tm: int = 0, om: int = 0,
    min_ts: float = 1.0,
    max_os: float = 0.0,
    relax_tm: float | None = None,
    relax_om: float | None = None,
) -> tuple[StrictnessClass, bool, str]:
    return classify_strictness(
        target_support_n=ts, target_total_n=tt,
        offtarget_support_n=os, offtarget_total_n=ot,
        target_missing_n=tm, offtarget_missing_n=om,
        min_target_support=min_ts,
        max_offtarget_support=max_os,
        relaxed_target_missing=relax_tm,
        relaxed_offtarget_missing=relax_om,
    )


# ---------------------------------------------------------------------------
# strict_complete
# ---------------------------------------------------------------------------

class TestStrictComplete:
    def test_all_targets_support_no_missing(self) -> None:
        klass, passes, reason = classify(ts=3, tt=3, os=0, ot=3)
        assert klass == StrictnessClass.STRICT_COMPLETE
        assert passes is True

    def test_single_target_single_offtarget(self) -> None:
        klass, passes, _ = classify(ts=1, tt=1, os=0, ot=1)
        assert klass == StrictnessClass.STRICT_COMPLETE
        assert passes is True

    def test_partial_target_support_at_threshold_100(self) -> None:
        # 2/3 targets supporting with min=1.0 → should NOT be strict_complete
        klass, passes, _ = classify(ts=2, tt=3, os=0, ot=3, min_ts=1.0)
        assert klass != StrictnessClass.STRICT_COMPLETE
        assert passes is False

    def test_zero_offtarget_support_required(self) -> None:
        klass, _, _ = classify(ts=3, tt=3, os=0, ot=5)
        assert klass == StrictnessClass.STRICT_COMPLETE


# ---------------------------------------------------------------------------
# contradicted
# ---------------------------------------------------------------------------

class TestContradicted:
    def test_any_offtarget_support_contradicts(self) -> None:
        klass, passes, reason = classify(ts=3, tt=3, os=1, ot=3)
        assert klass == StrictnessClass.CONTRADICTED
        assert passes is False
        assert "off-target support" in reason

    def test_all_offtarget_support_contradicts(self) -> None:
        klass, passes, _ = classify(ts=3, tt=3, os=3, ot=3)
        assert klass == StrictnessClass.CONTRADICTED
        assert passes is False

    def test_relaxed_max_offtarget_allows_some_support(self) -> None:
        # Allow up to 10% off-target support
        klass, passes, _ = classify(ts=3, tt=3, os=0, ot=10, max_os=0.1)
        assert klass == StrictnessClass.STRICT_COMPLETE
        assert passes is True

    def test_exceeding_relaxed_max_contradicts(self) -> None:
        klass, passes, _ = classify(ts=3, tt=3, os=2, ot=10, max_os=0.1)
        assert klass == StrictnessClass.CONTRADICTED
        assert passes is False


# ---------------------------------------------------------------------------
# strict_target_missing
# ---------------------------------------------------------------------------

class TestStrictTargetMissing:
    def test_one_target_missing(self) -> None:
        # 2 targets support, 1 missing, 0 offtarget missing
        klass, passes, _ = classify(ts=2, tt=3, os=0, ot=3, tm=1, om=0, min_ts=1.0)
        assert klass == StrictnessClass.STRICT_TARGET_MISSING
        assert passes is True

    def test_all_targets_missing(self) -> None:
        # All targets missing → cannot confirm target support
        klass, passes, _ = classify(ts=0, tt=3, os=0, ot=3, tm=3, om=0)
        assert klass == StrictnessClass.STRICT_TARGET_MISSING
        assert passes is False


# ---------------------------------------------------------------------------
# strict_offtarget_missing
# ---------------------------------------------------------------------------

class TestStrictOfftargetMissing:
    def test_one_offtarget_missing(self) -> None:
        klass, passes, _ = classify(ts=3, tt=3, os=0, ot=3, tm=0, om=1)
        assert klass == StrictnessClass.STRICT_OFFTARGET_MISSING
        assert passes is True


# ---------------------------------------------------------------------------
# strict_both_missing
# ---------------------------------------------------------------------------

class TestStrictBothMissing:
    def test_both_groups_have_missing(self) -> None:
        klass, passes, _ = classify(ts=2, tt=3, os=0, ot=3, tm=1, om=1, min_ts=1.0)
        assert klass == StrictnessClass.STRICT_BOTH_MISSING
        assert passes is True

    def test_all_missing_in_both(self) -> None:
        klass, passes, _ = classify(ts=0, tt=3, os=0, ot=3, tm=3, om=3)
        assert klass == StrictnessClass.STRICT_BOTH_MISSING
        assert passes is False


# ---------------------------------------------------------------------------
# relaxed_threshold
# ---------------------------------------------------------------------------

class TestRelaxedThreshold:
    def test_target_support_below_threshold(self) -> None:
        # 2/3 supporting with min=1.0 → relaxed
        klass, passes, reason = classify(ts=2, tt=3, os=0, ot=3, min_ts=1.0)
        assert klass == StrictnessClass.RELAXED_THRESHOLD
        assert passes is False

    def test_passes_at_lower_threshold(self) -> None:
        klass, passes, _ = classify(ts=2, tt=3, os=0, ot=3, min_ts=0.5)
        assert passes is True

    def test_relaxed_target_missing_downgrades(self) -> None:
        # 1/3 targets missing (33%) with relaxed_target_missing=0.2 → downgrade
        klass, passes, _ = classify(
            ts=2, tt=3, os=0, ot=3, tm=1, om=0,
            min_ts=1.0,
            relax_tm=0.2,
        )
        assert klass == StrictnessClass.RELAXED_THRESHOLD
        assert passes is True

    def test_relaxed_offtarget_missing_downgrades(self) -> None:
        klass, passes, _ = classify(
            ts=3, tt=3, os=0, ot=4, tm=0, om=2,
            min_ts=1.0,
            relax_om=0.4,
        )
        assert klass == StrictnessClass.RELAXED_THRESHOLD


# ---------------------------------------------------------------------------
# AllelePattern computed properties
# ---------------------------------------------------------------------------

class TestAllelePatternProperties:
    def test_target_support_fraction(self) -> None:
        p = build_allele_pattern(
            allele_key="chr1:100:A:T",
            target_support_n=3, target_total_n=3,
            offtarget_support_n=0, offtarget_total_n=3,
            target_missing_n=0, offtarget_missing_n=0,
        )
        assert p.target_support_fraction == 1.0
        assert p.offtarget_support_fraction == 0.0

    def test_none_fraction_when_all_missing(self) -> None:
        p = build_allele_pattern(
            allele_key="chr1:100:A:T",
            target_support_n=0, target_total_n=3,
            offtarget_support_n=0, offtarget_total_n=3,
            target_missing_n=3, offtarget_missing_n=0,
        )
        assert p.target_support_fraction is None

    def test_missing_fraction(self) -> None:
        p = build_allele_pattern(
            allele_key="key",
            target_support_n=2, target_total_n=4,
            offtarget_support_n=0, offtarget_total_n=4,
            target_missing_n=2, offtarget_missing_n=1,
        )
        assert p.target_missing_fraction == pytest.approx(0.5)
        assert p.offtarget_missing_fraction == pytest.approx(0.25)

    def test_build_allele_pattern_returns_correct_class(self) -> None:
        p = build_allele_pattern(
            allele_key="chr1:1000:G:C",
            target_support_n=5, target_total_n=5,
            offtarget_support_n=0, offtarget_total_n=5,
            target_missing_n=0, offtarget_missing_n=0,
        )
        assert p.strictness_class == StrictnessClass.STRICT_COMPLETE
        assert p.pattern_pass is True
        assert p.allele_key == "chr1:1000:G:C"


# ---------------------------------------------------------------------------
# StrictnessClass enum values
# ---------------------------------------------------------------------------

class TestStrictnessClassEnum:
    def test_all_values_exist(self) -> None:
        expected = {
            "strict_complete",
            "strict_target_missing",
            "strict_offtarget_missing",
            "strict_both_missing",
            "relaxed_threshold",
            "contradicted",
        }
        actual = {c.value for c in StrictnessClass}
        assert expected == actual

    def test_string_comparison(self) -> None:
        assert StrictnessClass.STRICT_COMPLETE == "strict_complete"
