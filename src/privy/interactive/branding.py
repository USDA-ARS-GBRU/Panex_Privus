"""Shared provenance branding for self-contained interactive dashboards."""

from __future__ import annotations

PANEX_PRIVUS_NAME = "Panex Privus"
PRIVY_NAME = "privy"
REPOSITORY_LABEL = "USDA-ARS-GBRU/Panex_Privus"
REPOSITORY_URL = "https://github.com/USDA-ARS-GBRU/Panex_Privus"

TOP_CREDIT_HTML = (
    'Generated with <a href="'
    + REPOSITORY_URL
    + '" rel="noreferrer">'
    + PRIVY_NAME
    + "</a>"
)

FOOTER_HTML = (
    "This self-contained HTML report was generated with "
    '<a href="'
    + REPOSITORY_URL
    + '" rel="noreferrer">'
    + PANEX_PRIVUS_NAME
    + " ("
    + PRIVY_NAME
    + ")</a>, open-source software for target-private comparative genomics. "
    'Repository: <a href="'
    + REPOSITORY_URL
    + '" rel="noreferrer">'
    + REPOSITORY_LABEL
    + "</a>."
)


def default_dashboard_title(mode: str) -> str:
    """Return the standard public-facing title for an interactive dashboard."""
    return f"{PANEX_PRIVUS_NAME} Interactive {mode} Report"
