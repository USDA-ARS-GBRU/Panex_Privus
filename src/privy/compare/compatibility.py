"""Allele and state compatibility checks for privy compare.

TODO (Phase 5): implement all compatibility dimensions.
"""

from __future__ import annotations

from privy.core.evidence import EvidenceRecord


def are_states_compatible(a: EvidenceRecord, b: EvidenceRecord) -> bool:
    """Return True if two evidence records have compatible states.

    Compatibility dimensions:
        - allele identity (same ALT if both from VCF)
        - support direction (both supporting OR both absent)
        - evidence class agreement

    TODO (Phase 5): implement.
    """
    raise NotImplementedError("are_states_compatible is not yet implemented.")
