"""TSV summary and ranked-hit table generators for privy report.

TODO (Phase 2): implement :func:`run_report` and supporting functions.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

from privy.core.config import PrivyConfig

log = logging.getLogger("privy.report.summary")


def run_report(
    hits: Path,
    regions: Optional[Path],
    evidence: Optional[Path],
    compare: Optional[Path],
    qc: Optional[Path],
    run_json: Optional[Path],
    cfg: PrivyConfig,
    fmt: str,
    top_n: int,
    include_qc: bool,
    include_strictness: bool,
    include_compare: bool,
    include_regions: bool,
    title: str,
    outdir: Path,
) -> None:
    """Generate ranked summaries and a human-readable report.

    TODO (Phase 2): implement ranked_hits.tsv, strictness_summary.tsv,
    support_summary.tsv, contradiction_summary.tsv, and Markdown report.
    """
    log.warning("privy report is not yet implemented.  Coming in Phase 2.")
    raise NotImplementedError(
        "privy report is not yet implemented.  "
        "The scaffold is in place; implementation arrives in Phase 2.\n\n"
        "Planned output files:\n"
        "  summary.tsv             — run-level summary\n"
        "  ranked_hits.tsv         — loci ranked by final_score\n"
        "  strictness_summary.tsv  — strictness class distribution\n"
        "  support_summary.tsv     — per-source support summary\n"
        "  contradiction_summary.tsv\n"
        "  report.md               — Markdown report\n"
        "  report.html             — HTML report (optional)\n"
    )
