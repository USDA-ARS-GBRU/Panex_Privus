"""Unit tests for privy.core.scoring."""

import pytest

from privy.core.patterns import StrictnessClass, build_allele_pattern
from privy.core.scoring import (
    ScoredHit,
    compute_discovery_score,
    compute_final_score,
    compute_penalty_score,
    compute_support_score,
    rank_scored_hits,
)


def strict_complete_pattern(n: int = 3) -> object:
    return build_allele_pattern(
        allele_key="chr1:100:A:T",
        target_support_n=n, target_total_n=n,
        offtarget_support_n=0, offtarget_total_n=n,
        target_missing_n=0, offtarget_missing_n=0,
    )


def contradicted_pattern() -> object:
    return build_allele_pattern(
        allele_key="chr1:200:G:C",
        target_support_n=3, target_total_n=3,
        offtarget_support_n=1, offtarget_total_n=3,
        target_missing_n=0, offtarget_missing_n=0,
    )


class TestComputeDiscoveryScore:
    def test_strict_complete_scores_higher_than_partial(self) -> None:
        p_complete = strict_complete_pattern()
        p_partial = build_allele_pattern(
            allele_key="key",
            target_support_n=2, target_total_n=3,
            offtarget_support_n=0, offtarget_total_n=3,
            target_missing_n=1, offtarget_missing_n=0,
            min_target_support=0.5,
        )
        assert compute_discovery_score(p_complete) > compute_discovery_score(p_partial)

    def test_zero_target_support_scores_low(self) -> None:
        p = build_allele_pattern(
            allele_key="key",
            target_support_n=0, target_total_n=3,
            offtarget_support_n=0, offtarget_total_n=3,
            target_missing_n=3, offtarget_missing_n=0,
        )
        score = compute_discovery_score(p)
        assert score >= 0.0

    def test_discovery_weight_scales_score(self) -> None:
        p = strict_complete_pattern()
        s1 = compute_discovery_score(p, discovery_weight=1.0)
        s2 = compute_discovery_score(p, discovery_weight=2.0)
        assert s2 == pytest.approx(s1 * 2.0)

    def test_qual_increases_score(self) -> None:
        p = strict_complete_pattern()
        no_qual = compute_discovery_score(p, variant_qual=None)
        with_qual = compute_discovery_score(p, variant_qual=60.0)
        assert with_qual > no_qual


class TestComputeSupportScore:
    def test_empty_evidence_returns_zero(self) -> None:
        assert compute_support_score([]) == 0.0

    def test_mean_of_values(self) -> None:
        score = compute_support_score([0.8, 0.6, 1.0])
        assert score == pytest.approx((0.8 + 0.6 + 1.0) / 3.0 * 1.0)

    def test_capped_at_weight(self) -> None:
        score = compute_support_score([1.0, 1.0, 1.0], support_weight=0.7)
        assert score == pytest.approx(0.7)


class TestComputePenaltyScore:
    def test_contradicted_max_penalty(self) -> None:
        p = contradicted_pattern()
        score = compute_penalty_score(p, penalty_weight=1.0)
        assert score == pytest.approx(1.0)

    def test_no_missingness_zero_penalty(self) -> None:
        p = strict_complete_pattern()
        score = compute_penalty_score(p)
        assert score == pytest.approx(0.0)

    def test_target_missingness_penalised(self) -> None:
        p = build_allele_pattern(
            allele_key="key",
            target_support_n=2, target_total_n=4,
            offtarget_support_n=0, offtarget_total_n=4,
            target_missing_n=2, offtarget_missing_n=0,
            min_target_support=0.5,
        )
        score = compute_penalty_score(p)
        assert score > 0.0


class TestComputeFinalScore:
    def test_additive_formula(self) -> None:
        assert compute_final_score(1.0, 0.5, 0.2) == pytest.approx(1.3)

    def test_negative_score_allowed(self) -> None:
        result = compute_final_score(0.1, 0.0, 0.8)
        assert result < 0.0

    def test_zero_inputs(self) -> None:
        assert compute_final_score(0.0, 0.0, 0.0) == pytest.approx(0.0)


class TestRankScoredHits:
    def _make_hit(self, locus_id: str, final_score: float) -> ScoredHit:
        return ScoredHit(
            locus_id=locus_id,
            discovery_score=final_score,
            support_score=0.0,
            penalty_score=0.0,
            final_score=final_score,
            rank=0,
            strictness_class="strict_complete",
            summary_label="test",
        )

    def test_rank_order(self) -> None:
        hits = [
            self._make_hit("L3", 0.5),
            self._make_hit("L1", 1.0),
            self._make_hit("L2", 0.8),
        ]
        ranked = rank_scored_hits(hits)
        assert ranked[0].locus_id == "L1"
        assert ranked[0].rank == 1
        assert ranked[1].locus_id == "L2"
        assert ranked[1].rank == 2
        assert ranked[2].locus_id == "L3"
        assert ranked[2].rank == 3

    def test_empty_input(self) -> None:
        assert rank_scored_hits([]) == []

    def test_single_hit_rank_one(self) -> None:
        hit = self._make_hit("L1", 0.9)
        ranked = rank_scored_hits([hit])
        assert ranked[0].rank == 1
