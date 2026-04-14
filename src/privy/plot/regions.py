"""Region summary and variant density plot generators.

TODO (Phase 2): implement region-level plot types.
"""

from __future__ import annotations

from pathlib import Path

from privy.core.locus import Locus


def plot_region_summary(
    region: Locus,
    constituent_loci: list[Locus],
    outdir: Path,
    width: float = 12.0,
    height: float = 4.0,
    dpi: int = 150,
    output_format: str = "png",
) -> Path:
    """Plot a region summary with variant density and strictness distribution.

    TODO (Phase 2): implement.
    """
    raise NotImplementedError("plot_region_summary is not yet implemented.")
