"""Core domain model for Panex Privus.

This subpackage contains the logic-centered domain objects that are
format-agnostic.  All evidence sources (VCF, BAM, GFA, XMFA) ultimately
feed into these objects.

Key types:
    - :class:`~privy.core.cohort.CohortDefinition`
    - :class:`~privy.core.locus.Locus`
    - :class:`~privy.core.patterns.AllelePattern`
    - :class:`~privy.core.patterns.StrictnessClass`
    - :class:`~privy.core.evidence.EvidenceRecord`
    - :class:`~privy.core.evidence.ComparisonRecord`
    - :class:`~privy.core.scoring.ScoredHit`
    - :class:`~privy.core.config.PrivyConfig`
"""
