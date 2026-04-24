"""Panex Privus — a comparative genomics toolkit for discovering target-private genomic signal.

CLI name: ``privy``

Primary use case:
    Given a target cohort and an off-target cohort, identify alleles and genomic
    regions present within the target group and absent from the off-target group,
    with explicit missingness reporting and optional cross-evidence support from
    BAM, GFA, and XMFA.

Commands:
    privy scan      — Discover target-private alleles and regions (VCF-first)
    privy compare   — Compare loci or regions across evidence sources
    privy report    — Generate ranked summaries and human-readable reports
    privy plot      — Create focused locus and region visualizations
    privy annotate  — Intersect private loci with GFF3 gene annotations
"""

from importlib.metadata import PackageNotFoundError, version
from pathlib import Path


def _version_from_checkout() -> str | None:
    """Return the local pyproject version when running from a source checkout."""
    pyproject = Path(__file__).resolve().parents[2] / "pyproject.toml"
    if not pyproject.exists():
        return None

    in_project = False
    for raw_line in pyproject.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if line == "[project]":
            in_project = True
            continue
        if in_project and line.startswith("["):
            return None
        if in_project and line.startswith("version"):
            _, _, value = line.partition("=")
            return value.strip().strip('"')
    return None


try:
    __version__: str = _version_from_checkout() or version("panex-privus")
except PackageNotFoundError:
    __version__ = "unknown"

__all__ = ["__version__"]
