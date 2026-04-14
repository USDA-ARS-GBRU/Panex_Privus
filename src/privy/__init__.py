"""Panex Privus — a comparative genomics toolkit for discovering target-private genomic signal.

CLI name: ``privy``

Primary use case:
    Given a target cohort and an off-target cohort, identify alleles and genomic
    regions present within the target group and absent from the off-target group,
    with explicit missingness reporting and optional cross-evidence support from
    BAM, GFA, and XMFA.

Commands:
    privy scan     — Discover target-private alleles and regions (VCF-first)
    privy compare  — Compare loci or regions across evidence sources
    privy report   — Generate ranked summaries and human-readable reports
    privy plot     — Create focused locus and region visualizations
"""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__: str = version("panex-privus")
except PackageNotFoundError:
    __version__ = "unknown"

__all__ = ["__version__"]
