"""GFA graph-context support layer.

Annotates candidate loci with path membership, junction proximity, and
local graph complexity from a GFA graph file.

TODO (Phase 4): implement all functions.
"""

from __future__ import annotations

from pathlib import Path

from privy.core.config import GfaConfig
from privy.core.evidence import EvidenceRecord
from privy.core.locus import Locus


def annotate_loci_with_gfa(
    loci: list[Locus],
    gfa_path: Path,
    cfg: GfaConfig,
) -> list[EvidenceRecord]:
    """Query a GFA graph around candidate loci and return context evidence.

    TODO (Phase 4): implement path membership and junction queries.
    """
    raise NotImplementedError("annotate_loci_with_gfa is not yet implemented.")
