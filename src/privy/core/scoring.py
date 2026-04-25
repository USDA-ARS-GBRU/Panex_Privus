"""ScoredHit and scoring functions for Panex Privus.

Scoring is transparent and additive::

    final_score = discovery_score + support_score - penalty_score

All components are in the range [0, ∞).  Weights are stored in
:class:`~privy.core.config.ScoringConfig` and written to ``run.json``
so every run is reproducible.

Design principle:
    Scores must be explainable.  A bioinformatics tool that produces a
    black-box rank is not useful for publication.  Every final_score must
    be decomposable into its three components.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from privy.core.patterns import AllelePattern, StrictnessClass


@dataclass
class ScoredHit:
    """Final ranked output for a locus, combining discovery and support evidence.

    Attributes:
        locus_id: Locus identifier.
        discovery_score: Score derived from the VCF cohort pattern.
        support_score: Score from secondary evidence (BAM/GFA/XMFA).
            Zero if no support evidence was collected.
        penalty_score: Deduction for missingness, contradiction, or
            ambiguity.
        final_score: ``discovery_score + support_score - penalty_score``.
        rank: Rank by ``final_score`` (1 = highest).
        strictness_class: Strictness class string from the source
            :class:`~privy.core.patterns.AllelePattern`.
        summary_label: Human-readable one-line description.
        metadata: Optional additional scoring components for transparency.
    """

    locus_id: str
    discovery_score: float
    support_score: float
    penalty_score: float
    final_score: float
    rank: int
    strictness_class: str
    summary_label: str
    metadata: dict[str, float] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Score computation functions
# ---------------------------------------------------------------------------

def compute_discovery_score(
    pattern: AllelePattern,
    variant_qual: float | None = None,
    discovery_weight: float = 1.0,
) -> float:
    """Compute a discovery score from an :class:`~privy.core.patterns.AllelePattern`.

    Components (all normalised to [0, 1] before weighting):
        - target support fraction
        - off-target exclusion fraction (1 - off_target_support)
        - private-allele specificity bonus for ``strict_complete``
        - optional variant quality contribution (capped and scaled)

    Args:
        pattern: Classified allele pattern from VCF.
        variant_qual: VCF QUAL value (optional; capped at 60 for normalisation).
        discovery_weight: Multiplier from :class:`~privy.core.config.ScoringConfig`.

    Returns:
        Non-negative discovery score.
    """
    target_frac = pattern.target_support_fraction or 0.0
    offtarget_frac = pattern.offtarget_support_fraction or 0.0
    offtarget_exclusion = 1.0 - offtarget_frac

    specificity_bonus = 0.2 if pattern.strictness_class == StrictnessClass.STRICT_COMPLETE else 0.0

    qual_component = 0.0
    if variant_qual is not None:
        # Scale quality: cap at 60, contribute up to 0.1
        qual_component = min(variant_qual / 60.0, 1.0) * 0.1

    raw = (target_frac + offtarget_exclusion) / 2.0 + specificity_bonus + qual_component
    return round(min(raw, 2.0) * discovery_weight, 6)


def compute_support_score(
    evidence_values: list[float],
    support_weight: float = 1.0,
) -> float:
    """Compute a support score from a list of normalised secondary evidence values.

    Each evidence value should be in [0, 1].  The support score is their mean,
    scaled by *support_weight*.

    Args:
        evidence_values: Normalised support metrics from BAM/GFA/XMFA.
        support_weight: Multiplier from :class:`~privy.core.config.ScoringConfig`.

    Returns:
        Non-negative support score (0.0 if no evidence values provided).
    """
    if not evidence_values:
        return 0.0
    mean_val = sum(evidence_values) / len(evidence_values)
    return round(min(mean_val, 1.0) * support_weight, 6)


def compute_penalty_score(
    pattern: AllelePattern,
    penalty_weight: float = 1.0,
) -> float:
    """Compute a penalty score for missingness and contradiction.

    Higher penalty = less reliable finding.

    Penalty components:
        - ``contradicted`` class: maximum penalty (1.0 * weight)
        - target missingness fraction: up to 0.4 * weight
        - off-target missingness fraction: up to 0.3 * weight

    Args:
        pattern: Classified allele pattern.
        penalty_weight: Multiplier from :class:`~privy.core.config.ScoringConfig`.

    Returns:
        Non-negative penalty score.
    """
    if pattern.strictness_class == StrictnessClass.CONTRADICTED:
        return round(1.0 * penalty_weight, 6)

    target_miss_pen = pattern.target_missing_fraction * 0.4
    offtarget_miss_pen = pattern.offtarget_missing_fraction * 0.3
    raw = target_miss_pen + offtarget_miss_pen
    return round(min(raw, 1.0) * penalty_weight, 6)


def compute_final_score(
    discovery_score: float,
    support_score: float,
    penalty_score: float,
) -> float:
    """Return the transparent additive final score.

    ``final_score = discovery_score + support_score - penalty_score``

    The result is not floored at zero — a negative final_score is a valid
    signal that the locus is contradicted or heavily penalised.
    """
    return round(discovery_score + support_score - penalty_score, 6)


def rank_scored_hits(hits: list[ScoredHit]) -> list[ScoredHit]:
    """Assign ranks to a list of :class:`ScoredHit` objects by ``final_score`` (descending).

    Modifies the ``rank`` attribute of each hit in-place and returns the
    sorted list.
    """
    hits_sorted = sorted(hits, key=lambda h: h.final_score, reverse=True)
    for i, hit in enumerate(hits_sorted, start=1):
        hit.rank = i
    return hits_sorted
