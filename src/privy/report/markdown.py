"""Markdown report generator for privy report.

TODO (Phase 2): implement :func:`render_markdown_report`.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any


def render_markdown_report(
    sections: dict[str, Any],
    title: str,
    outdir: Path,
) -> Path:
    """Render a Markdown report from assembled sections.

    Args:
        sections: Dict mapping section name to content (tables, strings, etc.).
        title: Report title.
        outdir: Output directory.

    Returns:
        Path to the written ``report.md``.

    TODO (Phase 2): implement template-based Markdown rendering.
    """
    raise NotImplementedError("render_markdown_report is not yet implemented.")
