"""Compatibility facade for the implemented ``privy compare`` backend.

The public CLI uses :mod:`privy.backends.compare` directly.  This module keeps
the older orchestration import path working for callers that still reach for
``privy.compare.engine.run_compare``.
"""

from __future__ import annotations

import logging
from pathlib import Path

from privy.backends.compare import run_compare as run_hits_compare
from privy.core.config import PrivyConfig

log = logging.getLogger("privy.compare.engine")


def run_compare(
    hits: Path | None,
    regions: Path | None,
    vcf: Path | None,
    bam: list[Path] | None,
    bam_manifest: Path | None,
    gfa: Path | None,
    xmfa: Path | None,
    source_a: Path | None,
    source_b: Path | None,
    cfg: PrivyConfig,
    mode: str,
    outdir: Path,
    write_compare_tsv: bool = True,
    write_summary_tsv: bool = True,
    write_json: bool = True,
) -> None:
    """Compare two scan ``hits.tsv`` files through the current backend.

    The legacy signature accepted many evidence-source paths.  The implemented
    compare backend currently supports scan-vs-scan comparison, so this facade
    resolves the two hits tables from ``source_a``/``source_b`` first, then from
    ``hits`` plus ``gfa``/``vcf`` when those paths point at TSV files.
    """
    hits_a, hits_b = _resolve_hit_tables(hits, vcf, gfa, source_a, source_b)
    if hits_a is None or hits_b is None:
        raise ValueError(
            "privy.compare.engine.run_compare requires two hits.tsv files. "
            "Pass them as source_a/source_b, or use privy.backends.compare.run_compare "
            "with hits_a and hits_b."
        )

    if mode not in {"scan_vs_scan", "vcf_vs_gfa", "auto"}:
        log.warning("Compare mode %s is handled as scan_vs_scan.", mode)

    run_hits_compare(
        hits_a=hits_a,
        hits_b=hits_b,
        outdir=outdir,
        cfg=cfg,
        write_compare_tsv=write_compare_tsv,
        write_summary_tsv=write_summary_tsv,
        write_json=write_json,
    )


def _resolve_hit_tables(
    hits: Path | None,
    vcf: Path | None,
    gfa: Path | None,
    source_a: Path | None,
    source_b: Path | None,
) -> tuple[Path | None, Path | None]:
    """Resolve legacy compare arguments to the two implemented hits tables."""
    first = source_a or hits
    second = source_b
    if second is None and gfa is not None and _looks_like_hits_tsv(gfa):
        second = gfa
    if second is None and vcf is not None and _looks_like_hits_tsv(vcf):
        second = vcf
    return first, second


def _looks_like_hits_tsv(path: Path) -> bool:
    return path.name == "hits.tsv" or path.suffix.lower() == ".tsv"
