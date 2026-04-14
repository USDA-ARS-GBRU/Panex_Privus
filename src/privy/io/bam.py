"""BAM depth and allele-level queries for Panex Privus.

BAM is a support layer — it provides read-level evidence at known candidate
loci, not a de-novo variant caller.

TODO (Phase 3):
    - Implement :func:`query_depth_at_locus` using pysam.AlignmentFile.count_coverage.
    - Implement :func:`query_allele_counts_at_locus` for ref/alt read counting.
    - Implement :func:`query_allele_fraction` wrapping count-based logic.
    - Add soft-clip / split-read summarisers near SV-like intervals.
    - Support manifest-driven multi-BAM iteration.
    - Validate BAM index (.bai or .csi) presence.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional


def validate_bam_index(bam_path: Path) -> None:
    """Raise FileNotFoundError if the BAM index is absent.

    TODO (Phase 3): implement.
    """
    raise NotImplementedError("validate_bam_index is not yet implemented.")


def query_depth_at_locus(
    bam_path: Path,
    contig: str,
    start: int,
    end: int,
    min_mapq: int = 20,
) -> list[int]:
    """Return per-position depth over [start, end) in a BAM file.

    TODO (Phase 3): implement using pysam.AlignmentFile.count_coverage.
    """
    raise NotImplementedError("query_depth_at_locus is not yet implemented.")


def query_allele_counts_at_locus(
    bam_path: Path,
    contig: str,
    pos: int,
    ref_allele: str,
    alt_allele: str,
    min_mapq: int = 20,
    min_baseq: int = 20,
) -> tuple[int, int, int]:
    """Return (ref_count, alt_count, other_count) at a SNP position.

    TODO (Phase 3): implement using pysam pileup.
    """
    raise NotImplementedError("query_allele_counts_at_locus is not yet implemented.")


def load_bam_manifest(manifest_path: Path) -> list[dict[str, str]]:
    """Load a TSV BAM manifest mapping file paths to sample names/groups.

    Expected columns: ``bam_path``, ``sample_id``, optionally ``group``.

    TODO (Phase 3): implement.
    """
    raise NotImplementedError("load_bam_manifest is not yet implemented.")
