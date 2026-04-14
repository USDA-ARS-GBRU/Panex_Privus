"""Miscellaneous utilities for Panex Privus."""

from __future__ import annotations

import datetime
from pathlib import Path


def now_iso() -> str:
    """Return the current UTC time as an ISO 8601 string."""
    return datetime.datetime.now(tz=datetime.timezone.utc).isoformat()


def ensure_dir(path: Path) -> Path:
    """Create *path* and all parents if they do not exist.  Return *path*."""
    path.mkdir(parents=True, exist_ok=True)
    return path


def format_bp(n: int) -> str:
    """Format a base-pair count with SI-like suffix for human readability.

    Examples::

        format_bp(500)      → "500 bp"
        format_bp(3_000)    → "3.0 kb"
        format_bp(5_000_000) → "5.0 Mb"
    """
    if n < 1_000:
        return f"{n} bp"
    if n < 1_000_000:
        return f"{n / 1_000:.1f} kb"
    return f"{n / 1_000_000:.1f} Mb"
