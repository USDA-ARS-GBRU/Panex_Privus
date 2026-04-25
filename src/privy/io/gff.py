"""GFF3 parser and interval index for privy annotate.

Coordinates are stored internally as 0-based half-open [start, end) intervals.
GFF3 uses 1-based closed [start, end] — conversion applied at parse time.
"""

from __future__ import annotations

import gzip
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path
from typing import IO

# Feature types used for classification; mRNA and intron are structural only
CLASSIFY_FEATURES: frozenset[str] = frozenset(
    {"gene", "CDS", "exon", "five_prime_UTR", "three_prime_UTR"}
)

# Feature types that indicate a sub-exonic (transcript-body) position
UTR_FEATURES: frozenset[str] = frozenset({"five_prime_UTR", "three_prime_UTR"})


@dataclass(slots=True)
class GffRecord:
    """A single parsed GFF3 record with 0-based half-open coordinates."""

    seqid: str
    feature_type: str
    start: int  # 0-based inclusive
    end: int    # 0-based exclusive
    strand: str
    gene_id: str  # Name= attribute value (empty string if absent)
    attrs: dict[str, str]


@dataclass
class AnnotationIndex:
    """In-memory lookup structure for GFF3 annotation.

    All interval lists are sorted by ``start`` to allow bisect-based queries.
    """

    # contig → sorted [(start, end, gene_name, strand)]
    genes: dict[str, list[tuple[int, int, str, str]]] = field(default_factory=dict)

    # contig → feature_type → sorted [(start, end)]
    sub_features: dict[str, dict[str, list[tuple[int, int]]]] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def _parse_attrs(raw: str) -> dict[str, str]:
    """Parse GFF3 attribute string into a dict."""
    attrs: dict[str, str] = {}
    for item in raw.strip().split(";"):
        item = item.strip()
        if not item or "=" not in item:
            continue
        key, _, val = item.partition("=")
        attrs[key.strip()] = val.strip()
    return attrs


def _open_gff(path: Path) -> IO[str]:
    """Return a file handle, handling .gz transparently."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, encoding="utf-8")


def parse_gff3(
    path: Path,
    feature_types: frozenset[str] | None = None,
) -> Iterator[GffRecord]:
    """Yield :class:`GffRecord` objects from a GFF3 file.

    Args:
        path: Path to GFF3 file (plain or .gz).
        feature_types: If given, only yield records whose feature type is in
            this set.  Pass ``None`` to yield all features.

    Yields:
        :class:`GffRecord` with 0-based half-open coordinates.
    """
    wanted = feature_types  # None → accept all
    with _open_gff(path) as fh:
        for raw_line in fh:
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            feat = parts[2]
            if wanted is not None and feat not in wanted:
                continue
            try:
                gff_start = int(parts[3])
                gff_end = int(parts[4])
            except ValueError:
                continue
            attrs = _parse_attrs(parts[8])
            gene_id = attrs.get("Name", attrs.get("ID", ""))
            yield GffRecord(
                seqid=parts[0],
                feature_type=feat,
                start=gff_start - 1,  # convert to 0-based
                end=gff_end,          # GFF3 end is inclusive → exclusive in 0-based
                strand=parts[6],
                gene_id=gene_id,
                attrs=attrs,
            )


# ---------------------------------------------------------------------------
# Index building
# ---------------------------------------------------------------------------

def build_annotation_index(path: Path) -> AnnotationIndex:
    """Parse a GFF3 file and return a queryable :class:`AnnotationIndex`.

    Only features in :data:`CLASSIFY_FEATURES` are indexed.

    Args:
        path: Path to GFF3 file (plain or .gz).

    Returns:
        Populated :class:`AnnotationIndex` with all interval lists sorted by
        start position.
    """
    idx = AnnotationIndex()

    for rec in parse_gff3(path, feature_types=CLASSIFY_FEATURES):
        contig = rec.seqid

        if rec.feature_type == "gene":
            bucket = idx.genes.setdefault(contig, [])
            bucket.append((rec.start, rec.end, rec.gene_id, rec.strand))
        else:
            contig_dict = idx.sub_features.setdefault(contig, {})
            sub_bucket = contig_dict.setdefault(rec.feature_type, [])
            sub_bucket.append((rec.start, rec.end))

    # Sort all buckets by start position for bisect queries
    for contig in idx.genes:
        idx.genes[contig].sort(key=lambda t: t[0])
    for contig in idx.sub_features:
        for feat in idx.sub_features[contig]:
            idx.sub_features[contig][feat].sort(key=lambda t: t[0])

    return idx


# ---------------------------------------------------------------------------
# Querying
# ---------------------------------------------------------------------------

def _overlapping_gene_intervals(
    intervals: list[tuple[int, int, str, str]],
    q_start: int,
    q_end: int,
) -> list[tuple[int, int, str, str]]:
    """Return gene intervals overlapping [q_start, q_end)."""
    return [iv for iv in intervals if iv[0] < q_end and iv[1] > q_start]


def _overlaps_sub_feature(
    intervals: list[tuple[int, int]],
    q_start: int,
    q_end: int,
) -> bool:
    """Return True if any sub-feature interval overlaps [q_start, q_end)."""
    for iv in intervals:
        if iv[0] < q_end and iv[1] > q_start:
            return True
    return False


def query_genes(
    idx: AnnotationIndex,
    contig: str,
    q_start: int,
    q_end: int,
) -> list[tuple[int, int, str, str]]:
    """Return all gene records overlapping [q_start, q_end) on *contig*.

    Returns:
        List of ``(start, end, gene_name, strand)`` tuples, possibly empty.
    """
    bucket = idx.genes.get(contig)
    if not bucket:
        return []
    return _overlapping_gene_intervals(bucket, q_start, q_end)


def query_sub_feature(
    idx: AnnotationIndex,
    contig: str,
    feature_type: str,
    q_start: int,
    q_end: int,
) -> bool:
    """Return True if any interval of *feature_type* overlaps the query locus."""
    contig_dict = idx.sub_features.get(contig)
    if not contig_dict:
        return False
    bucket = contig_dict.get(feature_type)
    if not bucket:
        return False
    return _overlaps_sub_feature(bucket, q_start, q_end)


# ---------------------------------------------------------------------------
# Contig alias normalisation
# ---------------------------------------------------------------------------

def load_contig_alias(path: Path) -> dict[str, str]:
    """Load a two-column TSV mapping source contig names to canonical names.

    Each line: ``<source_name>\\t<canonical_name>``
    Lines beginning with ``#`` are ignored.

    Args:
        path: Path to alias file.

    Returns:
        Dict mapping source names → canonical names.
    """
    alias: dict[str, str] = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                alias[parts[0]] = parts[1]
    return alias
