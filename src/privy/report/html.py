"""HTML report generator for privy report.

TODO (Phase 2): implement :func:`render_html_report` from Markdown source.
"""

from __future__ import annotations

from pathlib import Path


def render_html_report(markdown_path: Path, outdir: Path) -> Path:
    """Convert a Markdown report to HTML.

    Args:
        markdown_path: Path to ``report.md``.
        outdir: Output directory.

    Returns:
        Path to written ``report.html``.

    TODO (Phase 2): implement using markdown or mistune library.
    """
    raise NotImplementedError("render_html_report is not yet implemented.")
