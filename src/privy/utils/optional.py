"""Optional-dependency helpers for Privy's tiered dependency policy.

Privy keeps a strict separation so the core always works with a light install:

* **Tier 0 — pure-Python core.**  stdlib plus the small set of hard dependencies
  declared in ``pyproject.toml``.  Every core workflow runs with only these.
* **Tier 1 — optional Python extras** (``pip install 'panex-privus[full]'``):
  scipy, scikit-learn, gfapy, pyranges, scikit-allel.  Imported lazily via
  :func:`require`; their absence degrades a feature gracefully rather than
  breaking import of the package.
* **Tier 2 — external command-line tools** (odgi, vg, minimap2, OrthoFinder, …):
  never required.  Privy ingests their *output* or, where convenient, shells out
  only when the tool is detected via :func:`require_tool`.

This module centralises the "is it here? / give me a good error if not" logic so
every call site reports a consistent, actionable message.
"""

from __future__ import annotations

import importlib
import importlib.util
import shutil
from types import ModuleType

# Friendly install hints for the Tier-1 extras we know about.  Keyed by the
# *import* name (which can differ from the PyPI name, e.g. sklearn / allel).
_EXTRA_HINTS: dict[str, str] = {
    "scipy": "pip install 'panex-privus[full]'   (or: pip install scipy)",
    "sklearn": "pip install 'panex-privus[full]'   (or: pip install scikit-learn)",
    "gfapy": "pip install 'panex-privus[full]'   (or: pip install gfapy)",
    "pyranges": "pip install 'panex-privus[full]'   (or: pip install pyranges)",
    "allel": "pip install 'panex-privus[full]'   (or: pip install scikit-allel)",
}


class MissingDependencyError(RuntimeError):
    """Raised when a Tier-1 optional Python package is needed but not installed."""


class MissingToolError(RuntimeError):
    """Raised when a Tier-2 external command-line tool is needed but not found."""


# ---------------------------------------------------------------------------
# Tier-1: optional Python packages
# ---------------------------------------------------------------------------


def is_available(module_name: str) -> bool:
    """Return True if *module_name* can be imported, without importing it."""
    try:
        return importlib.util.find_spec(module_name) is not None
    except (ImportError, ValueError):
        # ValueError: a parent package is missing / spec lookup failed.
        return False


def require(module_name: str, *, feature: str, hint: str | None = None) -> ModuleType:
    """Import and return an optional module, or raise an actionable error.

    Args:
        module_name: The import name (e.g. ``"sklearn"``, ``"scipy"``).
        feature: Human-readable name of the Privy feature needing it, used in the
            error message (e.g. ``"DAPC population structure"``).
        hint: Override the default install hint.

    Raises:
        MissingDependencyError: If the module cannot be imported.
    """
    try:
        return importlib.import_module(module_name)
    except ImportError as exc:
        install = hint or _EXTRA_HINTS.get(module_name, f"pip install {module_name}")
        raise MissingDependencyError(
            f"{feature} requires the optional package '{module_name}', which is not "
            f"installed.\nInstall it with:\n    {install}"
        ) from exc


# ---------------------------------------------------------------------------
# Tier-2: external command-line tools
# ---------------------------------------------------------------------------


def tool_available(tool_name: str) -> bool:
    """Return True if *tool_name* is found on ``PATH``."""
    return shutil.which(tool_name) is not None


def require_tool(tool_name: str, *, feature: str, hint: str | None = None) -> str:
    """Return the resolved path to an external tool, or raise an actionable error.

    Args:
        tool_name: Executable name to look up on ``PATH`` (e.g. ``"odgi"``).
        feature: Human-readable name of the Privy feature needing it.
        hint: Optional install/usage hint appended to the error.

    Raises:
        MissingToolError: If the tool is not found on ``PATH``.
    """
    path = shutil.which(tool_name)
    if path is None:
        message = (
            f"{feature} requires the external tool '{tool_name}', which was not found "
            f"on PATH. Privy can also ingest its output instead of running it."
        )
        if hint:
            message += f"\n{hint}"
        raise MissingToolError(message)
    return path
