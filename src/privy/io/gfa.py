"""GFA graph-context reader for Panex Privus.

GFA is a structural context layer — it annotates whether candidate loci
fall near graph branches, bubbles, or junctions, and reports path membership.

TODO (Phase 4):
    - Implement GFA line parser (S/L/P/W record types).
    - Implement path membership query for a genomic interval.
    - Implement local graph complexity (junction count, branch density) near a locus.
    - Support GFA1 and GFA2; raise NotImplementedError for unsupported versions.
    - Index graph positions against a linear reference coordinate system.
"""

from __future__ import annotations

from pathlib import Path


def parse_gfa_header(gfa_path: Path) -> dict[str, str]:
    """Return metadata from the GFA H-line header.

    TODO (Phase 4): implement.
    """
    raise NotImplementedError("parse_gfa_header is not yet implemented.")


def query_path_membership(
    gfa_path: Path,
    contig: str,
    start: int,
    end: int,
) -> list[str]:
    """Return path names that overlap the given linear coordinate window.

    TODO (Phase 4): implement.
    """
    raise NotImplementedError("query_path_membership is not yet implemented.")


def query_junction_count(
    gfa_path: Path,
    contig: str,
    start: int,
    end: int,
    window_bp: int = 1000,
) -> int:
    """Return the number of graph junctions within *window_bp* of [start, end).

    TODO (Phase 4): implement.
    """
    raise NotImplementedError("query_junction_count is not yet implemented.")
