"""Core data model for the comparative-pangenome (synteny) layer.

These typed, frozen dataclasses are the shared vocabulary used across synteny
construction, coordinate projection, microhaplotypes, and visualization:

* :class:`GenomeInterval` — a half-open interval on one genome's contig.
* :class:`Anchor` — a query interval mapped to a target interval (the unifying
  unit; mirrors a PAF row, also produced by graph co-traversal or gene pairs).
* :class:`SyntenyBlock` — a collinear / inverted / translocated / duplicated run
  of anchors between two genomes (SyRI-style typing).
* :class:`SyntenyRegion` — overlapping blocks merged on a reference.
* :class:`ReferenceRange` — a PHG-style reference-anchored tiling interval whose
  members are local haplotypes (content-hashed IDs).
* :class:`ProjectionMap` — a region projected onto many references at once.

Coordinate convention: 0-based, half-open ``[start, end)`` — matching
:class:`~privy.core.locus.Locus`, the GFA parser, and PAF.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from enum import Enum

from privy.io.paf import PafRecord

# ---------------------------------------------------------------------------
# Enums
# ---------------------------------------------------------------------------


class AnchorSource(str, Enum):
    """Where an anchor came from."""

    GRAPH = "graph"          # shared segment co-traversal
    PAF = "paf"              # ingested alignment (odgi untangle / minimap2 / wfmash)
    GENE = "gene"            # orthogroup-constrained gene pair


class BlockType(str, Enum):
    """SyRI-style classification of a synteny block by alignment conformation."""

    COLLINEAR = "collinear"
    INVERSION = "inversion"
    TRANSLOCATION = "translocation"
    DUPLICATION = "duplication"
    UNALIGNED = "unaligned"


# ---------------------------------------------------------------------------
# Intervals and anchors
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class GenomeInterval:
    """A 0-based half-open interval on one genome's contig."""

    genome: str
    contig: str
    start: int
    end: int

    def __post_init__(self) -> None:
        if self.start < 0:
            raise ValueError(f"start must be >= 0, got {self.start}.")
        if self.end < self.start:
            raise ValueError(f"end ({self.end}) must be >= start ({self.start}).")

    @property
    def length(self) -> int:
        """Interval length in bp."""
        return self.end - self.start

    def overlaps(self, other: GenomeInterval) -> bool:
        """True when *other* is on the same contig and the intervals overlap."""
        return (
            self.contig == other.contig
            and self.start < other.end
            and other.start < self.end
        )


@dataclass(frozen=True)
class Anchor:
    """A query interval mapped to a target interval — the unifying synteny unit."""

    query: GenomeInterval
    target: GenomeInterval
    strand: str = "+"               # relative orientation, "+" or "-"
    score: float | None = None
    identity: float | None = None   # fraction in [0, 1]
    source: AnchorSource = AnchorSource.GRAPH
    name: str | None = None

    def __post_init__(self) -> None:
        if self.strand not in ("+", "-"):
            raise ValueError(f"strand must be '+' or '-', got {self.strand!r}.")

    @property
    def is_reverse(self) -> bool:
        """True for a reverse-strand (inverted) anchor."""
        return self.strand == "-"

    @classmethod
    def from_paf(cls, record: PafRecord, *, pansn_delimiter: str = "#") -> Anchor:
        """Build an :class:`Anchor` from a parsed :class:`~privy.io.paf.PafRecord`.

        PanSN-style names (``sample#hap#contig``) populate ``genome`` (sample) and
        ``contig`` (final component); plain names use the name for both.
        """
        q_genome, _q_hap, q_contig = split_pansn(record.query_name, pansn_delimiter)
        t_genome, _t_hap, t_contig = split_pansn(record.target_name, pansn_delimiter)
        return cls(
            query=GenomeInterval(q_genome, q_contig, record.query_start, record.query_end),
            target=GenomeInterval(t_genome, t_contig, record.target_start, record.target_end),
            strand=record.strand,
            score=float(record.mapping_quality),
            identity=record.blast_identity,
            source=AnchorSource.PAF,
        )


# ---------------------------------------------------------------------------
# Blocks and regions
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class SyntenyBlock:
    """A collinear (or rearranged) run of anchors between two genomes."""

    block_id: str
    query: GenomeInterval
    target: GenomeInterval
    strand: str
    block_type: BlockType = BlockType.COLLINEAR
    anchors: tuple[Anchor, ...] = ()
    score: float | None = None
    e_value: float | None = None

    def __post_init__(self) -> None:
        if self.strand not in ("+", "-"):
            raise ValueError(f"strand must be '+' or '-', got {self.strand!r}.")

    @property
    def n_anchors(self) -> int:
        """Number of anchors supporting the block."""
        return len(self.anchors)


@dataclass(frozen=True)
class SyntenyRegion:
    """Overlapping synteny blocks merged on a reference interval."""

    region_id: str
    reference: GenomeInterval
    blocks: tuple[SyntenyBlock, ...] = ()

    @property
    def n_blocks(self) -> int:
        """Number of member blocks."""
        return len(self.blocks)

    @property
    def genomes(self) -> tuple[str, ...]:
        """Distinct genomes participating in this region, in first-seen order."""
        seen: dict[str, None] = {}
        for block in self.blocks:
            seen.setdefault(block.query.genome, None)
            seen.setdefault(block.target.genome, None)
        return tuple(seen)


# ---------------------------------------------------------------------------
# Reference ranges (PHG concept) and projection results
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ReferenceRange:
    """A reference-anchored tiling interval whose members are local haplotypes.

    ``haplotypes`` maps a genome/path id to a content-hash (MD5-of-sequence)
    haplotype ID — the PHG/hVCF convention adopted for free dedup and stable
    cross-run identity.
    """

    range_id: str
    contig: str
    start: int
    end: int
    haplotypes: Mapping[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.end < self.start:
            raise ValueError(f"end ({self.end}) must be >= start ({self.start}).")

    @property
    def length(self) -> int:
        """Range length in bp."""
        return self.end - self.start

    @property
    def n_distinct_haplotypes(self) -> int:
        """Count of distinct haplotype IDs across members (allelic richness)."""
        return len(set(self.haplotypes.values()))


@dataclass(frozen=True)
class ProjectionMap:
    """A source region projected onto many reference coordinate systems at once.

    ``projections`` maps a target genome/path id to its projected interval, or
    ``None`` when the source region is absent from that genome.
    """

    source: str
    projections: Mapping[str, GenomeInterval | None] = field(default_factory=dict)

    def present_in(self) -> tuple[str, ...]:
        """Genome/path ids where the region projects (non-None), in insertion order."""
        return tuple(k for k, v in self.projections.items() if v is not None)

    def absent_in(self) -> tuple[str, ...]:
        """Genome/path ids where the region is absent (None), in insertion order."""
        return tuple(k for k, v in self.projections.items() if v is None)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def split_pansn(name: str, delimiter: str = "#") -> tuple[str, int | None, str]:
    """Split a PanSN-style ``sample<delim>hap<delim>contig`` name.

    Returns ``(sample, haplotype, contig)``.  For non-PanSN names the whole name
    is returned as both sample and contig with ``haplotype=None``.  A two-field
    ``sample<delim>contig`` name is also accepted (haplotype ``None``).
    """
    parts = name.split(delimiter)
    if len(parts) >= 3:
        hap: int | None
        try:
            hap = int(parts[1])
        except ValueError:
            hap = None
        return parts[0], hap, parts[-1]
    if len(parts) == 2:
        return parts[0], None, parts[1]
    return name, None, name
