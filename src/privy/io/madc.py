"""DArTag MADC importer — a lighter, amplicon-panel on-ramp to the allele model.

The Missing Allele Discovery Counts (MADC) report from the DArTag targeted-
genotyping platform is a long-format table: one row per (locus, microhaplotype
allele) carrying a per-sample read count, with an allele class (Ref / Alt /
RefMatch / AltMatch).  This reader parses that structure into per-locus allele
read-count records and offers simple dominant-allele calling.

Column names vary across panels/versions, so they are configurable; the defaults
target the standardized Breeding Insight MADC layout.  **Validate column mappings
against production MADC exports** — polyploid dosage genotype calling is out of
scope here (use polyRAD), this does straightforward dominant-allele calls.

Refs: Zhao et al. (2026) microhaplotype databases; github.com/Breeding-Insight/HapApp_utils.
"""

from __future__ import annotations

import csv
import gzip
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import TextIO

# Columns that are metadata, not per-sample read-count columns, by default.
_DEFAULT_META_COLUMNS = frozenset({
    "locus", "MarkerID", "marker", "CloneID", "clone",
    "AlleleID", "allele", "AlleleType", "AlleleClass", "class",
    "AlleleSequence", "sequence", "Sequence", "Ref", "Alt",
})


@dataclass(frozen=True)
class MadcAllele:
    """One microhaplotype allele at a locus, with per-sample read counts."""

    locus_id: str
    allele_id: str
    allele_class: str | None
    counts: dict[str, int] = field(default_factory=dict)


@dataclass(frozen=True)
class MadcLocus:
    """All alleles observed at one MADC locus."""

    locus_id: str
    samples: tuple[str, ...]
    alleles: list[MadcAllele] = field(default_factory=list)


def _open_text(path: Path) -> TextIO:
    with open(path, "rb") as probe:
        magic = probe.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, encoding="utf-8")


def read_madc(
    madc_path: Path,
    *,
    locus_col: str | None = None,
    allele_col: str | None = None,
    class_col: str | None = None,
    sample_cols: Sequence[str] | None = None,
    meta_columns: frozenset[str] = _DEFAULT_META_COLUMNS,
) -> list[MadcLocus]:
    """Parse a (optionally gzipped) DArTag MADC report into per-locus records.

    Args:
        locus_col / allele_col / class_col: column names; auto-detected from common
            aliases when None.
        sample_cols: explicit sample columns; when None, every column not in
            *meta_columns* (and not the detected id columns) is treated as a sample.

    Raises:
        ValueError: If locus/allele columns cannot be resolved.
    """
    with _open_text(madc_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("MADC file has no header row")
        fields = list(reader.fieldnames)
        locus_c = locus_col or _first_present(fields, ("MarkerID", "locus", "marker", "CloneID"))
        allele_c = allele_col or _first_present(fields, ("AlleleID", "allele"))
        class_c = class_col or _first_present(fields, ("AlleleType", "AlleleClass", "class"))
        if locus_c is None or allele_c is None:
            raise ValueError(
                "could not resolve MADC locus/allele columns; pass locus_col/allele_col"
            )
        if sample_cols is None:
            samples = [
                f for f in fields
                if f not in meta_columns and f not in {locus_c, allele_c, class_c}
            ]
        else:
            samples = list(sample_cols)

        grouped: dict[str, list[MadcAllele]] = {}
        for row in reader:
            locus_id = (row.get(locus_c) or "").strip()
            allele_id = (row.get(allele_c) or "").strip()
            if not locus_id or not allele_id:
                continue
            counts = {s: _to_int(row.get(s)) for s in samples}
            class_value = row.get(class_c) if class_c else None
            allele_class = class_value.strip() if class_value else None
            grouped.setdefault(locus_id, []).append(
                MadcAllele(
                    locus_id=locus_id,
                    allele_id=allele_id,
                    allele_class=allele_class,
                    counts=counts,
                )
            )

    return [
        MadcLocus(locus_id=lid, samples=tuple(samples), alleles=alleles)
        for lid, alleles in grouped.items()
    ]


def call_alleles(locus: MadcLocus, *, min_reads: int = 2) -> dict[str, str]:
    """Dominant-allele call per sample: the highest-read allele with ≥ *min_reads*.

    Samples with no allele meeting the threshold (or a tie) get no call.
    """
    calls: dict[str, str] = {}
    for sample in locus.samples:
        best_allele: str | None = None
        best_reads = 0
        tie = False
        for allele in locus.alleles:
            reads = allele.counts.get(sample, 0)
            if reads > best_reads:
                best_reads, best_allele, tie = reads, allele.allele_id, False
            elif reads == best_reads and reads > 0:
                tie = True
        if best_allele is not None and best_reads >= min_reads and not tie:
            calls[sample] = best_allele
    return calls


def _first_present(fields: Sequence[str], candidates: Sequence[str]) -> str | None:
    field_set = set(fields)
    for c in candidates:
        if c in field_set:
            return c
    return None


def _to_int(value: str | None) -> int:
    if value is None or value == "":
        return 0
    try:
        return int(float(value))
    except ValueError:
        return 0
