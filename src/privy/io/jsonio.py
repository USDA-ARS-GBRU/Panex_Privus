"""JSON I/O utilities for run.json and compare.json outputs.

``run.json`` is the primary reproducibility artefact: it records the full
resolved configuration, cohort definition, software version, and summary
statistics for every run.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def write_run_json(path: Path, data: dict[str, Any]) -> None:
    """Write run metadata to a JSON file with consistent 2-space indentation.

    Args:
        path: Destination path (will be created or overwritten).
        data: Serialisable dict.  Non-JSON-native types (``Path``,
            ``datetime``, etc.) are converted via ``default=str``.
    """
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2, default=str)
        fh.write("\n")


def read_run_json(path: Path) -> dict[str, Any]:
    """Read and return run metadata from a JSON file.

    Args:
        path: Path to ``run.json``.

    Returns:
        Parsed dict.

    Raises:
        FileNotFoundError: If *path* does not exist.
        json.JSONDecodeError: If the file is not valid JSON.
    """
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)  # type: ignore[no-any-return]


def write_compare_json(path: Path, data: dict[str, Any]) -> None:
    """Write compare metadata to a JSON file."""
    write_run_json(path, data)
