"""Reader for ``vg deconstruct`` VCFs — graph bubble/snarl sites on a reference.

``vg deconstruct`` emits a VCF with one site per snarl (bubble) projected onto a
chosen reference path; nesting is recorded in the ``LV`` (level) and ``PS``
(parent snarl) INFO fields, and ``AT`` carries the allele traversals.  This reader
extracts those sites so Privy can join graph structural variation to its synteny
and private-allele results (e.g. keep only top-level bubbles, like vcfbub).

Pure-stdlib, gzip-aware.  Coordinates 0-based half-open internally (POS is 1-based on disk).
Wiki: https://github.com/vgteam/vg/wiki/VCF-export-with-vg-deconstruct
"""

from __future__ import annotations

import gzip
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO


@dataclass(frozen=True)
class BubbleSite:
    """One snarl/bubble site from a vg deconstruct VCF."""

    contig: str
    pos: int                      # 0-based
    ref: str
    alts: list[str]
    level: int | None = None      # LV (bubble nesting level; 0 = top-level)
    parent_snarl: str | None = None   # PS
    allele_traversals: list[str] | None = None   # AT

    @property
    def n_alleles(self) -> int:
        """Number of alleles (REF + ALTs) at the site."""
        return 1 + len(self.alts)

    @property
    def is_top_level(self) -> bool:
        """True when this is a top-level bubble (LV 0 or unknown)."""
        return self.level in (0, None)


def _open_text(path: Path) -> TextIO:
    with open(path, "rb") as probe:
        magic = probe.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, encoding="utf-8")


def read_deconstruct_sites(
    vcf_path: Path,
    *,
    top_level_only: bool = False,
) -> Iterator[BubbleSite]:
    """Stream :class:`BubbleSite` records from a vg deconstruct VCF.

    Args:
        top_level_only: keep only LV==0 bubbles (vcfbub-style top-level filter).
    """
    with _open_text(vcf_path) as handle:
        for raw in handle:
            line = raw.rstrip("\n\r")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 8:
                continue
            info = _parse_info(fields[7])
            level = int(info["LV"]) if "LV" in info and info["LV"].lstrip("-").isdigit() else None
            site = BubbleSite(
                contig=fields[0],
                pos=int(fields[1]) - 1,
                ref=fields[3],
                alts=[] if fields[4] == "." else fields[4].split(","),
                level=level,
                parent_snarl=info.get("PS"),
                allele_traversals=info["AT"].split(",") if "AT" in info else None,
            )
            if top_level_only and not site.is_top_level:
                continue
            yield site


def _parse_info(info_field: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for entry in info_field.split(";"):
        if "=" in entry:
            key, value = entry.split("=", 1)
            out[key] = value
        elif entry:
            out[entry] = ""
    return out
