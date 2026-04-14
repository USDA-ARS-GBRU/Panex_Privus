"""Compare engine — orchestrates a full ``privy compare`` run.

TODO (Phase 5): implement :func:`run_compare`.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

from privy.core.config import PrivyConfig

log = logging.getLogger("privy.compare.engine")


def run_compare(
    hits: Optional[Path],
    regions: Optional[Path],
    vcf: Optional[Path],
    bam: Optional[list[Path]],
    bam_manifest: Optional[Path],
    gfa: Optional[Path],
    xmfa: Optional[Path],
    source_a: Optional[Path],
    source_b: Optional[Path],
    cfg: PrivyConfig,
    mode: str,
    outdir: Path,
    write_compare_tsv: bool = True,
    write_summary_tsv: bool = True,
    write_json: bool = True,
) -> None:
    """Orchestrate a cross-evidence comparison run.

    TODO (Phase 5): implement vcf_vs_bam, vcf_vs_gfa, vcf_vs_xmfa,
    scan_vs_scan, and multi_evidence comparison modes.
    """
    log.warning("privy compare is not yet implemented.  Coming in Phase 5.")
    raise NotImplementedError(
        "privy compare is not yet implemented.  "
        "The scaffold is in place; full implementation arrives in Phase 5.\n\n"
        "Planned output files:\n"
        "  compare.tsv         — per-locus comparison outcomes\n"
        "  compare_summary.tsv — match class distribution\n"
        "  compare.json        — run metadata\n"
    )
