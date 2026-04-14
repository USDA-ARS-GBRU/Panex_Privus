"""Locus — represents a genomic site or interval under evaluation.

A :class:`Locus` is the fundamental unit of analysis in Panex Privus.  It
can represent a single VCF record, a merged region of nearby variants, or an
interval derived from BAM/GFA/XMFA evidence.

Coordinate convention:
    0-based, half-open: [start, end).  This matches pysam and BED format.
    VCF 1-based POS values must be converted to 0-based when creating a Locus.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Optional


class LocusType(str, Enum):
    """Classification of the variant or region type represented by a Locus."""

    SNP = "snp"
    INDEL = "indel"
    SV = "sv"
    REGION = "region"
    GRAPH_REGION = "graph_region"
    ALIGNMENT_REGION = "alignment_region"


class PrimarySource(str, Enum):
    """Data source that originated a Locus."""

    VCF = "vcf"
    BAM = "bam"
    GFA = "gfa"
    XMFA = "xmfa"


@dataclass
class Locus:
    """A discrete genomic site or interval under evaluation.

    Attributes:
        locus_id: Unique identifier for this locus within a run (e.g., ``PPX000001``).
        contig: Chromosome or contig name (must match VCF/BAM header).
        start: 0-based start position, inclusive.
        end: 0-based end position, exclusive.
        locus_type: Classification of this locus.
        primary_source: Evidence source that first produced this locus.
        source_ids: Raw identifiers from the source (VCF ID field, read names, …).
        metadata: Optional key-value annotations attached to this locus.
    """

    locus_id: str
    contig: str
    start: int
    end: int
    locus_type: LocusType
    primary_source: PrimarySource
    source_ids: list[str] = field(default_factory=list)
    metadata: dict[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.start < 0:
            raise ValueError(f"Locus.start must be >= 0 (got {self.start}).")
        if self.end < self.start:
            raise ValueError(
                f"Locus.end must be >= start (start={self.start}, end={self.end})."
            )

    # ---------------------------------------------------------------- geometry

    @property
    def length(self) -> int:
        """Length of the locus in base pairs."""
        return self.end - self.start

    def overlaps(self, other: "Locus") -> bool:
        """Return True if this locus overlaps *other* on the same contig."""
        if self.contig != other.contig:
            return False
        return self.start < other.end and other.start < self.end

    def contains(self, other: "Locus") -> bool:
        """Return True if this locus fully contains *other*."""
        if self.contig != other.contig:
            return False
        return self.start <= other.start and self.end >= other.end

    def distance_to(self, other: "Locus") -> Optional[int]:
        """Return bp gap between non-overlapping same-contig loci.

        Returns:
            0 if overlapping, positive integer for gap, None for different contigs.
        """
        if self.contig != other.contig:
            return None
        if self.overlaps(other):
            return 0
        left, right = sorted([self, other], key=lambda loc: loc.start)
        return right.start - left.end

    def merge_with(self, other: "Locus", merged_id: Optional[str] = None) -> "Locus":
        """Return a new :class:`Locus` spanning both loci.

        The merged locus inherits this locus's ``primary_source`` and gets
        ``locus_type = LocusType.REGION``.

        Raises:
            ValueError: If loci are on different contigs.
        """
        if self.contig != other.contig:
            raise ValueError(
                f"Cannot merge loci on different contigs: "
                f"{self.contig!r} vs {other.contig!r}."
            )
        return Locus(
            locus_id=merged_id or f"{self.locus_id}__{other.locus_id}",
            contig=self.contig,
            start=min(self.start, other.start),
            end=max(self.end, other.end),
            locus_type=LocusType.REGION,
            primary_source=self.primary_source,
            source_ids=self.source_ids + other.source_ids,
        )

    # ------------------------------------------------------------------ coord

    @property
    def region_string(self) -> str:
        """samtools-style region string: ``contig:start+1-end`` (1-based, inclusive)."""
        return f"{self.contig}:{self.start + 1}-{self.end}"

    def __repr__(self) -> str:
        return (
            f"Locus({self.locus_id!r}, "
            f"{self.contig}:{self.start}-{self.end}, "
            f"{self.locus_type.value})"
        )
