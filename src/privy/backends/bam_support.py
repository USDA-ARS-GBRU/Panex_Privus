"""BAM support layer — read-level evidence at candidate loci.

BAM is not a discovery caller.  It provides corroborating or contradicting
evidence at loci already identified by the VCF backend.

TODO (Phase 3): implement all functions.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from privy.core.cohort import CohortDefinition
from privy.core.config import BamConfig
from privy.core.evidence import EvidenceRecord
from privy.core.locus import Locus


def annotate_loci_with_bam(
    loci: list[Locus],
    bam_paths: list[Path],
    cohort: CohortDefinition,
    cfg: BamConfig,
) -> list[EvidenceRecord]:
    """Query BAM files at candidate loci and return evidence records.

    TODO (Phase 3): implement depth, allele-count, and allele-fraction queries.
    """
    raise NotImplementedError("annotate_loci_with_bam is not yet implemented.")
