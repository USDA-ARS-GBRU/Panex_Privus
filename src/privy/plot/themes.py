"""Matplotlib theme and colour palette for Panex Privus plots.

Apply with :func:`apply_privy_theme` at the start of any plot function.
"""

from __future__ import annotations

# ---------------------------------------------------------------------------
# Colour palettes
# ---------------------------------------------------------------------------

STRICTNESS_COLOURS: dict[str, str] = {
    "strict_complete": "#2ecc71",
    "strict_target_missing": "#f39c12",
    "strict_offtarget_missing": "#e67e22",
    "strict_both_missing": "#e74c3c",
    "relaxed_threshold": "#9b59b6",
    "contradicted": "#2c3e50",
}

EVIDENCE_COLOURS: dict[str, str] = {
    "support": "#2ecc71",
    "absence": "#3498db",
    "ambiguous": "#f39c12",
    "contradiction": "#e74c3c",
    "uninformative": "#95a5a6",
}

MATCH_COLOURS: dict[str, str] = {
    "supported": "#2ecc71",
    "partially_supported": "#f39c12",
    "contradicted": "#e74c3c",
    "source_specific": "#3498db",
    "uninformative": "#95a5a6",
    "missing_data": "#bdc3c7",
}

PANGENOME_GROUP_COLOURS: dict[str, str] = {
    "full": "#933b41",
    "target": "#0868ac",
    "off_target": "#542788",
}

PANGENOME_CATEGORY_COLOURS: dict[str, str] = {
    "core": "#0868ac",
    "accessory": "#b35806",
    "private": "#542788",
    "absent": "#bdbdbd",
}

# ---------------------------------------------------------------------------
# Canonical ordering for consistent plot axes
# ---------------------------------------------------------------------------

STRICTNESS_ORDER: list[str] = [
    "strict_complete",
    "strict_target_missing",
    "strict_offtarget_missing",
    "strict_both_missing",
    "relaxed_threshold",
    "contradicted",
]

MATCH_ORDER: list[str] = [
    "supported",
    "partially_supported",
    "contradicted",
    "source_specific",
    "uninformative",
    "missing_data",
]

PANGENOME_GROUP_ORDER: list[str] = ["full", "target", "off_target"]

PANGENOME_CATEGORY_ORDER: list[str] = ["core", "accessory", "private", "absent"]


def apply_privy_theme() -> None:
    """Apply the Panex Privus matplotlib rc theme.

    Call once at the top of any plotting function (after setting backend).
    """
    try:
        import matplotlib.pyplot as plt  # noqa: PLC0415

        plt.rcParams.update({
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "axes.grid.axis": "x",
            "grid.alpha": 0.35,
            "grid.linewidth": 0.5,
            "grid.color": "#dddddd",
            "figure.dpi": 150,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "legend.framealpha": 0.85,
            "legend.edgecolor": "none",
            "legend.fontsize": 8,
        })
    except ImportError:
        pass  # matplotlib is optional at import time; will fail at plot-time
