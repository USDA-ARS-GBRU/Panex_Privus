"""Shared models for interactive dashboards.

The interactive browser displays user-facing genomic regions as 1-based,
closed intervals because that is the coordinate convention users type at the
CLI and see in VCF/GFF3 files. Backend scan logic elsewhere in Privy still uses
0-based half-open intervals internally.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

_REGION_RE = re.compile(r"^(?P<contig>[^:]+):(?P<start>[0-9,]+)-(?P<end>[0-9,]+)$")


@dataclass(frozen=True)
class FocusRegion:
    """A user-specified focus region in 1-based closed coordinates."""

    contig: str
    start: int
    end: int

    @property
    def length(self) -> int:
        """Return region length in bp."""
        return self.end - self.start + 1

    @property
    def label(self) -> str:
        """Return a compact display label."""
        return f"{self.contig}:{self.start}-{self.end}"

    @property
    def slug(self) -> str:
        """Return a filesystem-safe region slug."""
        contig = re.sub(r"[^A-Za-z0-9_.-]+", "_", self.contig).strip("_")
        return f"focus_{contig}_{self.start}_{self.end}"


def parse_focus_region(value: str) -> FocusRegion:
    """Parse ``CONTIG:START-END`` into a :class:`FocusRegion`.

    Coordinates are 1-based and inclusive. Commas in numeric coordinates are
    accepted for readability.
    """
    match = _REGION_RE.match(value.strip())
    if match is None:
        raise ValueError(
            f"Invalid focus region {value!r}; expected CONTIG:START-END, "
            "for example Gm15:1-4000000."
        )
    contig = match.group("contig").strip()
    start = int(match.group("start").replace(",", ""))
    end = int(match.group("end").replace(",", ""))
    if not contig:
        raise ValueError("Focus region contig cannot be empty.")
    if start < 1:
        raise ValueError(f"Focus region start must be >= 1: {value!r}.")
    if end < start:
        raise ValueError(f"Focus region end must be >= start: {value!r}.")
    return FocusRegion(contig=contig, start=start, end=end)
