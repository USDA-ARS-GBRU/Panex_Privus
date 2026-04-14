"""Summary diagnostic plot generators (strictness bars, support bars).

TODO (Phase 2): implement strictness_bar and support_bar plot types.
"""

from __future__ import annotations

from pathlib import Path


def plot_strictness_bar(
    strictness_counts: dict[str, int],
    outdir: Path,
    width: float = 8.0,
    height: float = 4.0,
    dpi: int = 150,
    output_format: str = "png",
) -> Path:
    """Bar chart of strictness class distribution.

    TODO (Phase 2): implement using :data:`~privy.plot.themes.STRICTNESS_COLOURS`.
    """
    raise NotImplementedError("plot_strictness_bar is not yet implemented.")


def plot_support_bar(
    support_counts: dict[str, int],
    outdir: Path,
    width: float = 8.0,
    height: float = 4.0,
    dpi: int = 150,
    output_format: str = "png",
) -> Path:
    """Stacked bar chart of evidence class distribution by source.

    TODO (Phase 2): implement.
    """
    raise NotImplementedError("plot_support_bar is not yet implemented.")
