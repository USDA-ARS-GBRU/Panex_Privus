"""AllelePattern and StrictnessClass — cohort-level VCF allele pattern logic.

This module contains the central biological logic of Panex Privus:
the private-allele discovery rule and the strictness classification system
that separates biological support from technical missingness.

The :class:`StrictnessClass` enum and :func:`classify_strictness` function
are the decision kernel.  Every VCF record passes through here before
any downstream scoring or region-building.

Design principle:
    Missingness must not be silently folded into pass/fail logic.  It must
    surface as a named :class:`StrictnessClass` so that biological signal and
    technical incompleteness remain separable in all outputs.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class StrictnessClass(str, Enum):
    """Strictness classes separate biological support from technical missingness.

    These values appear verbatim in the ``strictness_class`` column of
    ``hits.tsv`` and related outputs.  They must never be collapsed into a
    generic pass/fail.
    """

    STRICT_COMPLETE = "strict_complete"
    """All targets support the private allele; all off-targets are confidently
    absent; no missing calls among any required samples."""

    STRICT_TARGET_MISSING = "strict_target_missing"
    """Off-target exclusion holds, but one or more target samples are missing
    or uninformative.  The biological signal may be real but is incomplete."""

    STRICT_OFFTARGET_MISSING = "strict_offtarget_missing"
    """Target support holds, but one or more off-target samples are missing
    or uninformative.  Off-target exclusion cannot be fully confirmed."""

    STRICT_BOTH_MISSING = "strict_both_missing"
    """The pattern is otherwise consistent with target-private status, but
    missingness exists in both the target and off-target cohorts."""

    RELAXED_THRESHOLD = "relaxed_threshold"
    """Passes user-defined support thresholds but does not meet the criteria
    for any strict_* class (e.g., support fraction below strict threshold, or
    missingness above relaxed tolerance)."""

    CONTRADICTED = "contradicted"
    """The private-allele model fails: at least one off-target sample carries
    the allele above the allowed threshold, or target consistency fails."""


@dataclass
class AllelePattern:
    """Cohort-level summary of a candidate private allele from VCF evidence.

    This is the primary output of the private-allele discovery logic.  One
    :class:`AllelePattern` is produced per candidate alternate allele per
    VCF record.

    Attributes:
        allele_key: String key identifying the allele (``"contig:pos:ref:alt"``).
        target_support_n: Target samples carrying this allele (GT != 0/0 and not missing).
        target_total_n: Total target samples in the cohort.
        offtarget_support_n: Off-target samples carrying this allele.
        offtarget_total_n: Total off-target samples in the cohort.
        target_missing_n: Target samples with missing/uninformative genotype.
        offtarget_missing_n: Off-target samples with missing/uninformative genotype.
        strictness_class: Classification of this pattern (see :class:`StrictnessClass`).
        pattern_pass: Whether this allele passes discovery criteria.
        pattern_reason: Human-readable explanation of the classification decision.
    """

    allele_key: str
    target_support_n: int
    target_total_n: int
    offtarget_support_n: int
    offtarget_total_n: int
    target_missing_n: int
    offtarget_missing_n: int
    strictness_class: StrictnessClass
    pattern_pass: bool
    pattern_reason: str

    # -------------------------------------------------------------- computed

    @property
    def target_support_fraction(self) -> float | None:
        """Fraction of *called* (non-missing) target samples supporting the allele."""
        called = self.target_total_n - self.target_missing_n
        if called == 0:
            return None
        return self.target_support_n / called

    @property
    def offtarget_support_fraction(self) -> float | None:
        """Fraction of *called* (non-missing) off-target samples supporting the allele."""
        called = self.offtarget_total_n - self.offtarget_missing_n
        if called == 0:
            return None
        return self.offtarget_support_n / called

    @property
    def target_missing_fraction(self) -> float:
        """Fraction of target samples with missing genotype."""
        if self.target_total_n == 0:
            return 0.0
        return self.target_missing_n / self.target_total_n

    @property
    def offtarget_missing_fraction(self) -> float:
        """Fraction of off-target samples with missing genotype."""
        if self.offtarget_total_n == 0:
            return 0.0
        return self.offtarget_missing_n / self.offtarget_total_n


# ---------------------------------------------------------------------------
# Classification logic
# ---------------------------------------------------------------------------

def classify_strictness(
    target_support_n: int,
    target_total_n: int,
    offtarget_support_n: int,
    offtarget_total_n: int,
    target_missing_n: int,
    offtarget_missing_n: int,
    min_target_support: float = 1.0,
    max_offtarget_support: float = 0.0,
    relaxed_target_missing: float | None = None,
    relaxed_offtarget_missing: float | None = None,
) -> tuple[StrictnessClass, bool, str]:
    """Classify the strictness of a candidate private-allele pattern.

    This is the core decision function.  It is called once per candidate
    allele per VCF record.

    Args:
        target_support_n: Number of target samples carrying the allele.
        target_total_n: Total target cohort size.
        offtarget_support_n: Number of off-target samples carrying the allele.
        offtarget_total_n: Total off-target cohort size.
        target_missing_n: Target samples with missing/uninformative genotype.
        offtarget_missing_n: Off-target samples with missing/uninformative genotype.
        min_target_support: Required fraction of called target samples that
            must carry the allele (default 1.0 = all called targets).
        max_offtarget_support: Maximum fraction of called off-target samples
            allowed to carry the allele (default 0.0 = none allowed).
        relaxed_target_missing: If set, target missingness above this fraction
            downgrades to ``relaxed_threshold``.
        relaxed_offtarget_missing: If set, off-target missingness above this
            fraction downgrades to ``relaxed_threshold``.

    Returns:
        A tuple of ``(StrictnessClass, pattern_pass, reason_string)``.
    """
    target_called = target_total_n - target_missing_n
    offtarget_called = offtarget_total_n - offtarget_missing_n

    target_support_frac: float | None = (
        target_support_n / target_called if target_called > 0 else None
    )
    offtarget_support_frac: float | None = (
        offtarget_support_n / offtarget_called if offtarget_called > 0 else None
    )

    has_target_missing = target_missing_n > 0
    has_offtarget_missing = offtarget_missing_n > 0

    # ── Contradiction: off-target carries the allele above threshold ──────────
    if offtarget_support_frac is not None and offtarget_support_frac > max_offtarget_support:
        return (
            StrictnessClass.CONTRADICTED,
            False,
            (
                f"off-target support {offtarget_support_frac:.4f} "
                f"> threshold {max_offtarget_support}"
            ),
        )

    # ── All targets missing ───────────────────────────────────────────────────
    if target_support_frac is None:
        klass = (
            StrictnessClass.STRICT_BOTH_MISSING
            if has_offtarget_missing
            else StrictnessClass.STRICT_TARGET_MISSING
        )
        return (klass, False, "all target samples are missing/uninformative")

    # ── Target support below threshold ───────────────────────────────────────
    if target_support_frac < min_target_support:
        # Could still pass as relaxed if no off-target support
        if offtarget_support_frac is None or offtarget_support_frac <= max_offtarget_support:
            return (
                StrictnessClass.RELAXED_THRESHOLD,
                False,
                (
                    f"target support {target_support_frac:.4f} "
                    f"< threshold {min_target_support}"
                ),
            )
        return (
            StrictnessClass.CONTRADICTED,
            False,
            (
                f"target support {target_support_frac:.4f} below threshold "
                f"and off-target support {offtarget_support_frac:.4f} > {max_offtarget_support}"
            ),
        )

    # ── Target support passes threshold: classify by missingness ─────────────
    if not has_target_missing and not has_offtarget_missing:
        return (
            StrictnessClass.STRICT_COMPLETE,
            True,
            "all targets support; all off-targets confidently absent; no missing data",
        )

    if has_target_missing and has_offtarget_missing:
        klass = StrictnessClass.STRICT_BOTH_MISSING
    elif has_target_missing:
        klass = StrictnessClass.STRICT_TARGET_MISSING
    else:
        klass = StrictnessClass.STRICT_OFFTARGET_MISSING

    # Check whether missingness exceeds relaxed tolerances
    target_miss_frac = (
        target_missing_n / target_total_n if target_total_n > 0 else 0.0
    )
    offtarget_miss_frac = (
        offtarget_missing_n / offtarget_total_n if offtarget_total_n > 0 else 0.0
    )

    if relaxed_target_missing is not None and target_miss_frac > relaxed_target_missing:
        return (
            StrictnessClass.RELAXED_THRESHOLD,
            True,
            (
                f"target missingness {target_miss_frac:.4f} "
                f"> relaxed threshold {relaxed_target_missing}"
            ),
        )

    if relaxed_offtarget_missing is not None and offtarget_miss_frac > relaxed_offtarget_missing:
        return (
            StrictnessClass.RELAXED_THRESHOLD,
            True,
            (
                f"off-target missingness {offtarget_miss_frac:.4f} "
                f"> relaxed threshold {relaxed_offtarget_missing}"
            ),
        )

    return (
        klass,
        True,
        (
            f"passes with missingness: "
            f"target_missing={target_missing_n}, "
            f"offtarget_missing={offtarget_missing_n}"
        ),
    )


def build_allele_pattern(
    allele_key: str,
    target_support_n: int,
    target_total_n: int,
    offtarget_support_n: int,
    offtarget_total_n: int,
    target_missing_n: int,
    offtarget_missing_n: int,
    min_target_support: float = 1.0,
    max_offtarget_support: float = 0.0,
    relaxed_target_missing: float | None = None,
    relaxed_offtarget_missing: float | None = None,
) -> AllelePattern:
    """Construct a fully classified :class:`AllelePattern` from raw counts.

    This is the primary factory function used by the VCF scan backend.
    """
    strictness_class, pattern_pass, pattern_reason = classify_strictness(
        target_support_n=target_support_n,
        target_total_n=target_total_n,
        offtarget_support_n=offtarget_support_n,
        offtarget_total_n=offtarget_total_n,
        target_missing_n=target_missing_n,
        offtarget_missing_n=offtarget_missing_n,
        min_target_support=min_target_support,
        max_offtarget_support=max_offtarget_support,
        relaxed_target_missing=relaxed_target_missing,
        relaxed_offtarget_missing=relaxed_offtarget_missing,
    )
    return AllelePattern(
        allele_key=allele_key,
        target_support_n=target_support_n,
        target_total_n=target_total_n,
        offtarget_support_n=offtarget_support_n,
        offtarget_total_n=offtarget_total_n,
        target_missing_n=target_missing_n,
        offtarget_missing_n=offtarget_missing_n,
        strictness_class=strictness_class,
        pattern_pass=pattern_pass,
        pattern_reason=pattern_reason,
    )
