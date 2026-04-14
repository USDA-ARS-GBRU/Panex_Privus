"""BED region file reader/writer for Panex Privus.

BED files are used for restricting scans to genomic windows and for
exporting region intervals.

TODO (Phase 2):
    - Implement :func:`read_bed` for 3-column and 6-column BED files.
    - Implement :func:`write_bed` from a list of :class:`~privy.core.locus.Locus`.
    - Support gzipped BED input.
"""

from __future__ import annotations

from pathlib import Path

from privy.core.locus import Locus


def read_bed(bed_path: Path) -> list[tuple[str, int, int]]:
    """Read a BED file and return (contig, start, end) tuples.

    Coordinates are 0-based, half-open — same as BED spec and
    :class:`~privy.core.locus.Locus`.

    TODO (Phase 2): implement.
    """
    raise NotImplementedError("read_bed is not yet implemented.")


def write_bed(loci: list[Locus], bed_path: Path) -> None:
    """Write a list of Loci to a 3-column BED file.

    TODO (Phase 2): implement.
    """
    raise NotImplementedError("write_bed is not yet implemented.")
