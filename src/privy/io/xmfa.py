"""XMFA alignment-corroboration reader for Panex Privus.

XMFA (eXtended Multi-FastA Alignment) blocks are used as optional
alignment-level corroboration for candidate loci identified by VCF.

TODO (Phase 5):
    - Implement XMFA block parser (header ``=`` separators, strand-aware coords).
    - Implement reference-anchored coordinate extraction.
    - Implement gap-aware alignment state evaluation at candidate positions.
    - Support interval aggregation across discriminant columns.
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path


@dataclass
class XmfaBlock:
    """A single XMFA alignment block (between two ``=`` separators).

    TODO (Phase 5): flesh out with strand, genome ID, and sequence fields.
    """

    genome_id: str
    contig: str
    start: int  # 0-based
    end: int    # 0-based, exclusive
    strand: str  # '+' or '-'
    sequence: str


def iter_xmfa_blocks(xmfa_path: Path) -> Iterator[list[XmfaBlock]]:
    """Yield alignment blocks from an XMFA file.

    Each yielded value is a list of :class:`XmfaBlock` objects, one per
    genome in the alignment block.

    TODO (Phase 5): implement.
    """
    raise NotImplementedError("iter_xmfa_blocks is not yet implemented.")


def find_blocks_overlapping(
    xmfa_path: Path,
    contig: str,
    start: int,
    end: int,
    reference_genome_id: str | None = None,
) -> list[list[XmfaBlock]]:
    """Return alignment blocks where the reference genome overlaps [start, end).

    TODO (Phase 5): implement with appropriate indexing strategy.
    """
    raise NotImplementedError("find_blocks_overlapping is not yet implemented.")
