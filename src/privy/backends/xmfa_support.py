"""XMFA alignment-corroboration support layer.

Provides alignment-based corroboration for candidate loci from XMFA
whole-genome alignment files.

TODO (Phase 5): implement all functions.
"""

from __future__ import annotations

from pathlib import Path

from privy.core.config import XmfaConfig
from privy.core.evidence import EvidenceRecord
from privy.core.locus import Locus


def annotate_loci_with_xmfa(
    loci: list[Locus],
    xmfa_path: Path,
    cfg: XmfaConfig,
) -> list[EvidenceRecord]:
    """Find XMFA alignment blocks overlapping candidate loci and return evidence.

    TODO (Phase 5): implement gap-aware corroboration logic.
    """
    raise NotImplementedError("annotate_loci_with_xmfa is not yet implemented.")
