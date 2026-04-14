"""Per-locus evidence panel and genotype heatmap plots.

TODO (Phase 2): implement :func:`run_plot` and individual plot types.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

from privy.core.config import PrivyConfig

log = logging.getLogger("privy.plot.loci")


def run_plot(
    hits: Path,
    regions: Optional[Path],
    evidence: Optional[Path],
    vcf: Optional[Path],
    bam: Optional[list[Path]],
    bam_manifest: Optional[Path],
    gfa: Optional[Path],
    xmfa: Optional[Path],
    cfg: PrivyConfig,
    locus_id: Optional[str],
    region_id: Optional[str],
    top_n: Optional[int],
    contig: Optional[str],
    region: Optional[str],
    plot_type: str,
    width: float,
    height: float,
    dpi: int,
    output_format: str,
    show_labels: bool,
    outdir: Path,
) -> None:
    """Generate focused plots for loci, regions, or summaries.

    TODO (Phase 2): implement locus_panel, region_summary, genotype_heatmap,
    strictness_bar, support_bar, and depth_panel plot types.
    """
    log.warning("privy plot is not yet implemented.  Coming in Phase 2.")
    raise NotImplementedError(
        "privy plot is not yet implemented.  "
        "The scaffold is in place; implementation arrives in Phase 2."
    )
