"""Match-class assignment for privy compare.

Classifies loci into one of six :class:`~privy.core.evidence.MatchClass`
values based on interval overlap and evidence compatibility.

TODO (Phase 5): implement :func:`classify_comparison`.
"""

from __future__ import annotations

from privy.core.evidence import EvidenceRecord, MatchClass
from privy.core.locus import Locus


def classify_comparison(
    query_locus: Locus,
    query_evidence: list[EvidenceRecord],
    target_locus: Locus | None,
    target_evidence: list[EvidenceRecord],
    coordinate_overlap: float,
    breakpoint_tolerance_bp: int = 200,
    require_state_compatibility: bool = False,
) -> MatchClass:
    """Assign a :class:`~privy.core.evidence.MatchClass` to a locus comparison.

    Decision logic (in priority order):
        1. ``missing_data``      — no target locus or evidence at all.
        2. ``contradicted``      — overlap present but state incompatible.
        3. ``supported``         — high overlap + compatible states.
        4. ``partially_supported`` — moderate overlap or weak compatibility.
        5. ``uninformative``     — overlap but no usable state signal.
        6. ``source_specific``   — no overlap found.

    TODO (Phase 5): implement.
    """
    raise NotImplementedError("classify_comparison is not yet implemented.")
