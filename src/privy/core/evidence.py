"""EvidenceRecord, ComparisonRecord, and associated enumerations.

These are the format-agnostic normalisation types at the heart of Panex
Privus.  Every data source (VCF, BAM, GFA, XMFA) maps its outputs into
:class:`EvidenceRecord` instances so that all downstream logic is
format-agnostic.

:class:`ComparisonRecord` captures the result of a cross-source comparison,
as produced by ``privy compare``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any


class EvidenceClass(str, Enum):
    """Classification of a piece of evidence relative to target-private status.

    Every :class:`EvidenceRecord` carries exactly one of these classes.
    """

    SUPPORT = "support"
    """Evidence positively supports target-private status at this locus."""

    ABSENCE = "absence"
    """Evidence confirms absence of the allele/signal in off-target samples,
    supporting target-private status."""

    AMBIGUOUS = "ambiguous"
    """Evidence exists but cannot be classified as supporting or contradicting
    (e.g., low coverage, multi-allelic complexity, or alignment uncertainty)."""

    CONTRADICTION = "contradiction"
    """Evidence contradicts target-private status (e.g., off-target BAM reads
    carrying the private allele)."""

    UNINFORMATIVE = "uninformative"
    """Evidence is present but provides no useful signal for or against
    target-private status at this locus (e.g., a graph node with no path
    information)."""


class SourceType(str, Enum):
    """Data source that produced an :class:`EvidenceRecord`."""

    VCF = "vcf"
    BAM = "bam"
    GFA = "gfa"
    XMFA = "xmfa"


class MatchClass(str, Enum):
    """Cross-source comparison outcome for a locus.

    Used in :class:`ComparisonRecord` and in ``compare.tsv``.
    """

    SUPPORTED = "supported"
    """Two or more evidence sources agree that the locus is target-private."""

    PARTIALLY_SUPPORTED = "partially_supported"
    """Sources partially agree — overlap present but evidence compatibility
    is incomplete or weak."""

    CONTRADICTED = "contradicted"
    """A secondary source contradicts the primary source's finding at this
    locus (e.g., VCF says private; BAM shows off-target reads)."""

    SOURCE_SPECIFIC = "source_specific"
    """The locus exists in only one evidence source; no counterpart in the
    comparison source."""

    UNINFORMATIVE = "uninformative"
    """The comparison source is present but provides no usable signal
    (e.g., insufficient depth, alignment gaps, ambiguous graph structure)."""

    MISSING_DATA = "missing_data"
    """One or both sources have no data at this locus."""


@dataclass
class EvidenceRecord:
    """A normalised piece of evidence for a locus from any data source.

    All heterogeneous evidence flows into this form before reaching the
    comparison or scoring layers.

    Attributes:
        locus_id: Identifier of the :class:`~privy.core.locus.Locus` this
            evidence pertains to.
        source_type: Which data source produced this record.
        evidence_class: Classification of the evidence relative to
            target-private status.
        metric_name: Name of the metric being reported (e.g., ``"depth"``,
            ``"allele_fraction"``, ``"path_count"``).
        metric_value: Numeric value for the metric.
        sample_id: Optional sample identifier (for per-sample evidence).
        group_id: Optional group identifier (for cohort-aggregated evidence).
        qualifiers: Optional auxiliary key-value data (read counts, MAPQ, etc.).
        provenance: Human-readable description of how this record was produced.
    """

    locus_id: str
    source_type: SourceType
    evidence_class: EvidenceClass
    metric_name: str
    metric_value: float
    sample_id: str | None = None
    group_id: str | None = None
    qualifiers: dict[str, Any] = field(default_factory=dict)
    provenance: str = ""

    def is_supporting(self) -> bool:
        """Return True if this evidence supports target-private status."""
        return self.evidence_class in (EvidenceClass.SUPPORT, EvidenceClass.ABSENCE)

    def is_contradicting(self) -> bool:
        """Return True if this evidence contradicts target-private status."""
        return self.evidence_class == EvidenceClass.CONTRADICTION


@dataclass
class ComparisonRecord:
    """The outcome of a cross-source comparison for a locus.

    Produced by ``privy compare`` when two evidence sources are evaluated
    for the same locus.

    Attributes:
        locus_id: Identifier of the locus being compared.
        source_a: Primary evidence source.
        source_b: Secondary evidence source.
        match_class: Overall comparison outcome (see :class:`MatchClass`).
        coordinate_overlap: Reciprocal overlap fraction [0, 1].
        state_compatibility: Whether allele/state calls are compatible.
        support_summary: Human-readable summary of supporting evidence.
        contradiction_summary: Human-readable summary of contradictory evidence.
        comparison_score: Numeric concordance score in [0, 1].
        metadata: Optional auxiliary fields.
    """

    locus_id: str
    source_a: SourceType
    source_b: SourceType
    match_class: MatchClass
    coordinate_overlap: float
    state_compatibility: bool
    support_summary: str
    contradiction_summary: str
    comparison_score: float
    metadata: dict[str, Any] = field(default_factory=dict)
