"""Matplotlib theme and colour palette for Panex Privus plots.

Apply with :func:`apply_privy_theme` at the start of any plot function.
"""

from __future__ import annotations

# Colour palette — aligned to strictness classes and evidence classes
STRICTNESS_COLOURS: dict[str, str] = {
    "strict_complete": "#2ecc71",        # green
    "strict_target_missing": "#f39c12",  # amber
    "strict_offtarget_missing": "#e67e22",  # orange
    "strict_both_missing": "#e74c3c",    # red
    "relaxed_threshold": "#9b59b6",      # purple
    "contradicted": "#2c3e50",           # dark
}

EVIDENCE_COLOURS: dict[str, str] = {
    "support": "#2ecc71",
    "absence": "#3498db",
    "ambiguous": "#f39c12",
    "contradiction": "#e74c3c",
    "uninformative": "#95a5a6",
}


def apply_privy_theme() -> None:
    """Apply the Panex Privus matplotlib rc theme.

    Call once at the top of any plotting function.

    TODO (Phase 2): expand with publication-ready rcParams.
    """
    try:
        import matplotlib.pyplot as plt  # noqa: PLC0415

        plt.rcParams.update(
            {
                "font.family": "DejaVu Sans",
                "font.size": 10,
                "axes.spines.top": False,
                "axes.spines.right": False,
                "figure.dpi": 150,
                "savefig.dpi": 300,
                "savefig.bbox": "tight",
            }
        )
    except ImportError:
        pass  # matplotlib is optional at import time; will fail at plot-time
