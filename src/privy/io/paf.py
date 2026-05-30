"""PAF / BEDPE reader and writer for Panex Privus.

PAF (Pairwise mApping Format) is the lingua franca for alignment-derived anchors
across the minimap2 / wfmash / odgi-untangle ecosystem.  Panex Privus uses it as
the interchange format for synteny *anchors*: a region of one sequence mapped to a
region of another.  This module is pure-Python (stdlib only) and streams records so
it scales to large pangenome alignments.

PAF specification (12 mandatory tab-delimited columns)
------------------------------------------------------
=====  =======  ====================================================
Col    Type     Field
=====  =======  ====================================================
1      str      query name
2      int      query length
3      int      query start (0-based, closed)
4      int      query end   (0-based, open)
5      char     strand, ``+`` or ``-``
6      str      target name
7      int      target length
8      int      target start (0-based, closed)
9      int      target end   (0-based, open)
10     int      residue matches (number of identical bases)
11     int      alignment block length
12     int      mapping quality (0-255; 255 = missing)
=====  =======  ====================================================

Optional SAM-style typed tags follow as ``TAG:TYPE:VALUE``.  The ones Privy cares
about are ``cg:Z:`` (CIGAR), ``cs:Z:`` (difference string), ``tp:A:`` (alignment
type), ``NM:i:`` (edit distance), and divergence (``dv:f:`` / ``de:f:``).

Coordinate convention: 0-based, half-open ``[start, end)`` — identical to PAF
itself and to :class:`~privy.core.locus.Locus`, so no conversion is needed.

BEDPE
-----
``odgi untangle`` can emit BEDPE instead of PAF.  :func:`parse_bedpe` reads the
standard 10-column BEDPE layout (``chrom1 start1 end1 chrom2 start2 end2 name
score strand1 strand2`` + optional extra columns).
"""

from __future__ import annotations

import gzip
import logging
from collections.abc import Iterable, Iterator
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, TextIO

log = logging.getLogger("privy.io.paf")

PAF_MANDATORY_COLUMNS = 12


# ---------------------------------------------------------------------------
# Domain objects
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class PafTag:
    """One optional SAM-style ``TAG:TYPE:VALUE`` field from a PAF line."""

    tag: str
    type: str   # one of i f A Z H B
    value: Any

    def __str__(self) -> str:
        return f"{self.tag}:{self.type}:{_format_tag_value(self.type, self.value)}"


@dataclass
class PafRecord:
    """One alignment record (one line) from a PAF file.

    Coordinates are 0-based half-open, matching the PAF spec and Privy's
    :class:`~privy.core.locus.Locus` domain object.
    """

    query_name: str
    query_length: int
    query_start: int
    query_end: int
    strand: str
    target_name: str
    target_length: int
    target_start: int
    target_end: int
    residue_matches: int
    alignment_block_length: int
    mapping_quality: int
    tags: dict[str, PafTag] = field(default_factory=dict)

    # -- derived convenience accessors -------------------------------------

    @property
    def is_reverse(self) -> bool:
        """True when the alignment is on the reverse strand (``-``)."""
        return self.strand == "-"

    @property
    def blast_identity(self) -> float:
        """BLAST-like identity = residue matches / alignment block length.

        Returns 0.0 when the alignment block length is zero (degenerate line).
        """
        if self.alignment_block_length <= 0:
            return 0.0
        return self.residue_matches / self.alignment_block_length

    @property
    def query_aligned_length(self) -> int:
        """Number of query bases spanned by the alignment."""
        return self.query_end - self.query_start

    @property
    def target_aligned_length(self) -> int:
        """Number of target bases spanned by the alignment."""
        return self.target_end - self.target_start

    @property
    def cigar(self) -> str | None:
        """The ``cg:Z:`` CIGAR string, or None if absent."""
        value = self.get_tag("cg")
        return None if value is None else str(value)

    @property
    def cs(self) -> str | None:
        """The ``cs:Z:`` difference string, or None if absent."""
        value = self.get_tag("cs")
        return None if value is None else str(value)

    def get_tag(self, tag: str, default: Any = None) -> Any:
        """Return the value of optional *tag*, or *default* if absent."""
        entry = self.tags.get(tag)
        return entry.value if entry is not None else default


@dataclass
class BedpeRecord:
    """One record from a BEDPE file (e.g. ``odgi untangle`` output).

    Coordinates are 0-based half-open on both mates.
    """

    chrom1: str
    start1: int
    end1: int
    chrom2: str
    start2: int
    end2: int
    name: str | None = None
    score: str | None = None
    strand1: str | None = None
    strand2: str | None = None
    extra: list[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Tag helpers
# ---------------------------------------------------------------------------


def _coerce_tag_value(type_char: str, raw: str) -> Any:
    """Coerce a raw tag value string to a Python type based on the SAM type char."""
    if type_char == "i":
        return int(raw)
    if type_char == "f":
        return float(raw)
    # A (char), Z (string), H (hex), B (array) are kept as their raw string form.
    return raw


def _format_tag_value(type_char: str, value: Any) -> str:
    """Render a coerced tag value back to its on-disk string form."""
    if type_char == "i":
        return str(int(value))
    if type_char == "f":
        return repr(float(value)) if isinstance(value, float) else str(value)
    return str(value)


def _parse_tags(fields: list[str], line_num: int) -> dict[str, PafTag]:
    """Parse trailing ``TAG:TYPE:VALUE`` fields into a tag dictionary."""
    tags: dict[str, PafTag] = {}
    for raw in fields:
        parts = raw.split(":", 2)
        if len(parts) != 3:
            raise ValueError(
                f"Malformed PAF optional field {raw!r} on line {line_num}: "
                "expected TAG:TYPE:VALUE."
            )
        tag, type_char, value = parts
        tags[tag] = PafTag(tag=tag, type=type_char, value=_coerce_tag_value(type_char, value))
    return tags


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------


def _open_text(path: Path) -> TextIO:
    """Open a plain or gzipped text file for reading."""
    with open(path, "rb") as probe:
        magic = probe.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, encoding="utf-8")


# ---------------------------------------------------------------------------
# PAF parsing / writing
# ---------------------------------------------------------------------------


def parse_paf(paf_path: Path, *, skip_malformed: bool = False) -> Iterator[PafRecord]:
    """Stream :class:`PafRecord` objects from a (optionally gzipped) PAF file.

    Args:
        paf_path: Path to a ``.paf`` or ``.paf.gz`` file.
        skip_malformed: If True, log and skip lines that cannot be parsed
            instead of raising.  Default False (strict).

    Yields:
        One :class:`PafRecord` per non-blank line.

    Raises:
        ValueError: When a line has fewer than 12 mandatory columns or a
            mandatory integer field is non-numeric (unless *skip_malformed*).
    """
    with _open_text(paf_path) as handle:
        for line_num, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n\r")
            if not line or line.startswith("#"):
                continue
            try:
                yield _parse_paf_line(line, line_num)
            except ValueError:
                if skip_malformed:
                    log.warning("Skipping malformed PAF line %d", line_num)
                    continue
                raise


def _parse_paf_line(line: str, line_num: int) -> PafRecord:
    """Parse one PAF line into a :class:`PafRecord`."""
    fields = line.split("\t")
    if len(fields) < PAF_MANDATORY_COLUMNS:
        raise ValueError(
            f"PAF line {line_num} has {len(fields)} columns; "
            f"at least {PAF_MANDATORY_COLUMNS} are required."
        )
    try:
        record = PafRecord(
            query_name=fields[0],
            query_length=int(fields[1]),
            query_start=int(fields[2]),
            query_end=int(fields[3]),
            strand=fields[4],
            target_name=fields[5],
            target_length=int(fields[6]),
            target_start=int(fields[7]),
            target_end=int(fields[8]),
            residue_matches=int(fields[9]),
            alignment_block_length=int(fields[10]),
            mapping_quality=int(fields[11]),
            tags=_parse_tags(fields[PAF_MANDATORY_COLUMNS:], line_num),
        )
    except ValueError as exc:
        raise ValueError(f"PAF line {line_num}: {exc}") from exc
    if record.strand not in ("+", "-"):
        raise ValueError(
            f"PAF line {line_num}: strand must be '+' or '-', got {record.strand!r}."
        )
    return record


def format_paf_record(record: PafRecord) -> str:
    """Render a :class:`PafRecord` to a single PAF line (no trailing newline)."""
    columns = [
        record.query_name,
        str(record.query_length),
        str(record.query_start),
        str(record.query_end),
        record.strand,
        record.target_name,
        str(record.target_length),
        str(record.target_start),
        str(record.target_end),
        str(record.residue_matches),
        str(record.alignment_block_length),
        str(record.mapping_quality),
    ]
    columns.extend(str(tag) for tag in record.tags.values())
    return "\t".join(columns)


def write_paf(records: Iterable[PafRecord], paf_path: Path) -> int:
    """Write *records* to a (optionally gzipped) PAF file.

    A ``.gz`` suffix triggers gzip compression.

    Returns:
        The number of records written.
    """
    paf_path = Path(paf_path)
    opener = (
        gzip.open(paf_path, "wt", encoding="utf-8")
        if paf_path.suffix == ".gz"
        else open(paf_path, "w", encoding="utf-8")
    )
    count = 0
    with opener as handle:
        for record in records:
            handle.write(format_paf_record(record))
            handle.write("\n")
            count += 1
    return count


# ---------------------------------------------------------------------------
# BEDPE parsing
# ---------------------------------------------------------------------------

BEDPE_MANDATORY_COLUMNS = 6


def parse_bedpe(bedpe_path: Path, *, skip_malformed: bool = False) -> Iterator[BedpeRecord]:
    """Stream :class:`BedpeRecord` objects from a (optionally gzipped) BEDPE file.

    Reads the standard layout ``chrom1 start1 end1 chrom2 start2 end2`` with
    optional ``name score strand1 strand2`` and any further columns captured in
    :attr:`BedpeRecord.extra`.  Coordinates are 0-based half-open.
    """
    with _open_text(bedpe_path) as handle:
        for line_num, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n\r")
            if not line or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            if len(fields) < BEDPE_MANDATORY_COLUMNS:
                if skip_malformed:
                    log.warning("Skipping malformed BEDPE line %d", line_num)
                    continue
                raise ValueError(
                    f"BEDPE line {line_num} has {len(fields)} columns; "
                    f"at least {BEDPE_MANDATORY_COLUMNS} are required."
                )
            try:
                yield BedpeRecord(
                    chrom1=fields[0],
                    start1=int(fields[1]),
                    end1=int(fields[2]),
                    chrom2=fields[3],
                    start2=int(fields[4]),
                    end2=int(fields[5]),
                    name=fields[6] if len(fields) > 6 else None,
                    score=fields[7] if len(fields) > 7 else None,
                    strand1=fields[8] if len(fields) > 8 else None,
                    strand2=fields[9] if len(fields) > 9 else None,
                    extra=fields[10:],
                )
            except ValueError:
                if skip_malformed:
                    log.warning("Skipping malformed BEDPE line %d", line_num)
                    continue
                raise
