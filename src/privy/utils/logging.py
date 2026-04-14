"""Logging configuration for Panex Privus.

All privy loggers live under the ``privy`` namespace.  Use
:func:`get_logger` within any module to obtain a namespaced logger that
inherits the root handler configured by :func:`configure_logging`.

Example::

    from privy.utils.logging import get_logger
    log = get_logger("backends.vcf_scan")
    log.info("Scanning contig %s", contig)
"""

from __future__ import annotations

import logging
import sys


def configure_logging(level: str = "info", quiet: bool = False) -> None:
    """Configure the root ``privy`` logger.

    Adds a single ``StreamHandler`` to stderr.  Idempotent: calling this
    function multiple times does not add duplicate handlers.

    Args:
        level: Logging level name: ``debug``, ``info``, ``warning``, ``error``.
        quiet: If True, override *level* with ``WARNING``.
    """
    numeric_level = getattr(logging, level.upper(), logging.INFO)
    if quiet:
        numeric_level = logging.WARNING

    root_logger = logging.getLogger("privy")
    if root_logger.handlers:
        # Reconfigure if already set up (e.g., from a second CLI invocation)
        root_logger.handlers.clear()

    handler = logging.StreamHandler(sys.stderr)
    handler.setLevel(numeric_level)
    formatter = logging.Formatter(
        fmt="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S",
    )
    handler.setFormatter(formatter)

    root_logger.setLevel(numeric_level)
    root_logger.addHandler(handler)
    root_logger.propagate = False


def get_logger(name: str) -> logging.Logger:
    """Return a logger named ``privy.<name>``.

    Args:
        name: Dotted sub-name (e.g., ``"backends.vcf_scan"``).

    Returns:
        A :class:`logging.Logger` instance.
    """
    return logging.getLogger(f"privy.{name}")
