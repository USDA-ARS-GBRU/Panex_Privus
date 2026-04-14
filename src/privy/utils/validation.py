"""Input validation helpers for Panex Privus.

TODO (Phase 2): implement file-existence, format, and cross-input checks.
"""

from __future__ import annotations

from pathlib import Path


def require_file(path: Path, label: str = "File") -> None:
    """Raise FileNotFoundError if *path* does not exist.

    Args:
        path: Path to check.
        label: Human-readable label for error messages.
    """
    if not path.exists():
        raise FileNotFoundError(f"{label} not found: {path}")


def require_suffix(path: Path, *suffixes: str) -> None:
    """Raise ValueError if *path* does not end with one of *suffixes*.

    Args:
        path: Path to check.
        suffixes: Acceptable suffixes (e.g., ``".vcf.gz"``, ``".bcf"``).
    """
    name = path.name
    if not any(name.endswith(s) for s in suffixes):
        raise ValueError(
            f"Expected one of {suffixes} for {path}, got suffix {path.suffix!r}."
        )
