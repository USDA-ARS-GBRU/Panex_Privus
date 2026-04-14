"""Global run-state shared across CLI subcommands.

The ``_state`` singleton is populated by the root ``app.callback()`` in
:mod:`privy.cli.main` before any subcommand executes.  Subcommands read from
``get_state()`` after the callback has run.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


@dataclass
class _GlobalState:
    """Mutable container for options resolved before any subcommand runs."""

    config_path: Optional[Path] = None
    project_name: Optional[str] = None
    outdir: Path = field(default_factory=lambda: Path("."))
    threads: int = 1
    log_level: str = "info"
    quiet: bool = False


# Module-level singleton — mutated by the CLI callback, read by subcommands.
_state: _GlobalState = _GlobalState()


def get_state() -> _GlobalState:
    """Return the shared global run-state singleton."""
    return _state
