"""GFA (Graphical Fragment Assembly) graph reader for Panex Privus.

Parses GFA1 and GFA1.1 files into an in-memory graph with prebuilt inverted
indices for efficient cohort-support queries.

Supported record types
----------------------
H   Header — version and metadata tags
S   Segment — sequence node with optional coordinate tags
L   Link — directed edge between two segments
P   Path — ordered traversal (GFA1 classic, one per sample/haplotype)
W   Walk — ordered traversal with explicit reference coordinates (GFA1.1)

Coordinate conventions
----------------------
Segment positions come from optional tags on S-lines (minigraph/PGGB output):
    SN:Z:<contig>   reference sequence name
    SO:i:<offset>   0-based start offset on the reference
    LN:i:<length>   segment length in bp

W-line coordinates (``seq_start``, ``seq_end``) are 0-based half-open, matching
the Locus domain object.

Path-name conventions
---------------------
Pangenome tools (minigraph-cactus, PGGB) name P-lines as
``SAMPLE#HAP_INDEX#CONTIG``.  Plain names (no ``#``) are also supported.
W-lines carry the sample name as a separate field, so no parsing is needed.
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

log = logging.getLogger("privy.io.gfa")


# ---------------------------------------------------------------------------
# Domain objects
# ---------------------------------------------------------------------------


@dataclass
class GfaSegment:
    """One S-line from a GFA file."""

    name: str
    sequence: str          # raw bases, or "*" if not stored
    length: int            # from LN:i: tag, or len(sequence) if bases present
    tags: dict[str, Any] = field(default_factory=dict)
    # Linear reference coordinates (from SN / SO / LN tags if present)
    ref_contig: str | None = None   # SN:Z: value
    ref_start: int | None = None    # SO:i: value, 0-based
    ref_end: int | None = None      # SO + LN, 0-based exclusive


@dataclass
class GfaLink:
    """One L-line from a GFA file."""

    from_seg: str
    from_orient: str    # "+" or "-"
    to_seg: str
    to_orient: str
    overlap: str        # CIGAR string or "*"
    tags: dict[str, Any] = field(default_factory=dict)


@dataclass
class GfaPath:
    """One P-line from a GFA file (GFA1 classic paths)."""

    name: str                     # full path name, e.g. "HG002#1#chr1"
    sample: str                   # extracted sample name
    haplotype: int | None         # extracted haplotype index (None if not encoded)
    segment_names: list[str]
    orientations: list[str]       # "+" or "-" per segment, parallel to segment_names
    tags: dict[str, Any] = field(default_factory=dict)


@dataclass
class GfaWalkStep:
    """One step in a W-line walk path."""

    segment: str
    orient: str   # "+" or "-"


@dataclass
class GfaWalk:
    """One W-line from a GFA file (GFA1.1 pangenome walks)."""

    sample: str
    hap_index: int
    seq_id: str      # reference contig name
    seq_start: int   # 0-based inclusive
    seq_end: int     # 0-based exclusive
    steps: list[GfaWalkStep]
    tags: dict[str, Any] = field(default_factory=dict)


@dataclass
class GfaGraph:
    """Parsed GFA graph with prebuilt inverted indices for cohort queries.

    Build via :func:`parse_gfa`; do not construct directly.

    Attributes:
        segments: Mapping of segment name → :class:`GfaSegment`.
        links: All L-lines in file order.
        paths: Mapping of path name → :class:`GfaPath` (from P-lines).
        walks: All W-lines in file order.
        header_tags: Key→value dict from the H-line.
        segment_to_paths: Inverted index — segment name → path names.
        segment_to_walks: Inverted index — segment name → walk indices.
        sample_to_paths: sample name → path names.
        sample_to_walks: sample name → walk indices.
    """

    segments: dict[str, GfaSegment]
    links: list[GfaLink]
    paths: dict[str, GfaPath]
    walks: list[GfaWalk]
    header_tags: dict[str, str]

    # Inverted indices — populated by _build_indices()
    segment_to_paths: dict[str, list[str]] = field(default_factory=dict)
    segment_to_walks: dict[str, list[int]] = field(default_factory=dict)
    sample_to_paths: dict[str, list[str]] = field(default_factory=dict)
    sample_to_walks: dict[str, list[int]] = field(default_factory=dict)

    # Position index — contig → sorted list of (start, end, seg_name)
    # Populated only for segments that carry SN/SO/LN tags.
    _contig_segments: dict[str, list[tuple[int, int, str]]] = field(
        default_factory=dict
    )


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def parse_gfa(gfa_path: Path) -> GfaGraph:
    """Parse a GFA1 or GFA1.1 file and return a fully indexed :class:`GfaGraph`.

    Reads S, L, P, W, and H lines.  Unknown line types (J, C, E, F) are
    silently skipped.  Compressed GFA (.gfa.gz) is not supported — decompress
    first.

    Args:
        gfa_path: Path to a GFA file.

    Returns:
        A :class:`GfaGraph` with all inverted indices built.

    Raises:
        FileNotFoundError: If *gfa_path* does not exist.
        ValueError: If the file declares GFA version 2 (not supported).
    """
    if not gfa_path.exists():
        raise FileNotFoundError(f"GFA file not found: {gfa_path}")

    segments: dict[str, GfaSegment] = {}
    links: list[GfaLink] = []
    paths: dict[str, GfaPath] = {}
    walks: list[GfaWalk] = []
    header_tags: dict[str, str] = {}

    with open(gfa_path) as fh:
        for line_num, raw_line in enumerate(fh, 1):
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            rec_type = fields[0]

            if rec_type == "H":
                _parse_h_line(fields, header_tags)
            elif rec_type == "S":
                seg = _parse_s_line(fields, line_num)
                segments[seg.name] = seg
            elif rec_type == "L":
                links.append(_parse_l_line(fields, line_num))
            elif rec_type == "P":
                path = _parse_p_line(fields, line_num)
                paths[path.name] = path
            elif rec_type == "W":
                walks.append(_parse_w_line(fields, line_num))
            # GFA2 / extension lines silently skipped

    if header_tags.get("VN", "").startswith("2"):
        raise ValueError(
            f"GFA version 2 is not supported (found VN:{header_tags['VN']}). "
            "Only GFA1 and GFA1.1 are supported."
        )

    graph = GfaGraph(
        segments=segments,
        links=links,
        paths=paths,
        walks=walks,
        header_tags=header_tags,
    )
    _build_indices(graph)

    log.info(
        "Parsed GFA: %d segments, %d links, %d paths, %d walks, %d samples",
        len(segments),
        len(links),
        len(paths),
        len(walks),
        len(set(graph.sample_to_paths) | set(graph.sample_to_walks)),
    )
    return graph


def get_gfa_samples(gfa_path: Path) -> list[str]:
    """Return a sorted list of sample names present in a GFA file.

    Sample names are extracted from:
    - P-line path names (``SAMPLE#hap#contig`` or plain ``SAMPLE``).
    - W-line ``sample_id`` fields.

    Args:
        gfa_path: Path to a GFA1/1.1 file.

    Returns:
        Sorted list of unique sample name strings.
    """
    graph = parse_gfa(gfa_path)
    return sorted(set(graph.sample_to_paths) | set(graph.sample_to_walks))


def query_segments_at_locus(
    graph: GfaGraph,
    contig: str,
    start: int,
    end: int,
) -> list[GfaSegment]:
    """Return segments whose reference coordinates overlap ``[start, end)``.

    Only segments that carry SN/SO/LN tags are indexed.  Returns an empty list
    for contigs with no coordinate-tagged segments.

    Args:
        graph: Parsed :class:`GfaGraph`.
        contig: Contig/chromosome name.
        start: Query start (0-based inclusive).
        end: Query end (0-based exclusive).

    Returns:
        List of overlapping :class:`GfaSegment` objects (unsorted).
    """
    if contig not in graph._contig_segments:
        return []
    result: list[GfaSegment] = []
    for seg_start, seg_end, seg_name in graph._contig_segments[contig]:
        if seg_start >= end:
            break
        if seg_end > start:
            seg = graph.segments.get(seg_name)
            if seg is not None:
                result.append(seg)
    return result


def get_samples_traversing_segment(
    graph: GfaGraph,
    seg_name: str,
) -> set[str]:
    """Return the set of sample names whose paths or walks traverse *seg_name*.

    Checks both P-lines (via ``segment_to_paths``) and W-lines (via
    ``segment_to_walks``).

    Args:
        graph: Parsed :class:`GfaGraph`.
        seg_name: Segment name to query.

    Returns:
        Set of sample name strings.
    """
    samples: set[str] = set()
    for path_name in graph.segment_to_paths.get(seg_name, []):
        path = graph.paths.get(path_name)
        if path is not None:
            samples.add(path.sample)
    for walk_idx in graph.segment_to_walks.get(seg_name, []):
        samples.add(graph.walks[walk_idx].sample)
    return samples


def get_samples_present_at_locus(
    graph: GfaGraph,
    contig: str,
    start: int,
    end: int,
) -> set[str]:
    """Return samples that have ANY graph coverage overlapping ``[start, end)``.

    "Present" means the sample has at least one path or walk that touches the
    locus — regardless of which segment they traverse there.  Samples that are
    present but do not traverse a specific segment are counted as *absent* for
    that segment; samples with no coverage at all are counted as *missing*.

    Coverage is detected via:
    1. W-lines whose ``seq_start``/``seq_end`` interval overlaps the locus.
    2. P-lines whose traversed segments have SN/SO tags overlapping the locus.

    If neither W-lines nor coordinate-tagged segments are available for the
    contig, returns an empty set (all samples treated as missing).

    Args:
        graph: Parsed :class:`GfaGraph`.
        contig: Contig/chromosome name.
        start: Locus start (0-based inclusive).
        end: Locus end (0-based exclusive).

    Returns:
        Set of sample name strings with coverage at the locus.
    """
    samples: set[str] = set()

    # W-lines: check interval overlap directly
    for walk in graph.walks:
        if walk.seq_id == contig and walk.seq_start < end and walk.seq_end > start:
            samples.add(walk.sample)

    # P-lines: check if any traversed segment overlaps the locus
    if contig in graph._contig_segments:
        overlapping_segs: set[str] = set()
        for seg_start, seg_end, seg_name in graph._contig_segments[contig]:
            if seg_start >= end:
                break
            if seg_end > start:
                overlapping_segs.add(seg_name)
        for seg_name in overlapping_segs:
            for path_name in graph.segment_to_paths.get(seg_name, []):
                path = graph.paths.get(path_name)
                if path is not None:
                    samples.add(path.sample)

    return samples


def extract_cohort_segment_counts(
    graph: GfaGraph,
    seg_name: str,
    seg_contig: str | None,
    seg_start: int | None,
    seg_end: int | None,
    target_samples: list[str],
    offtarget_samples: list[str],
) -> tuple[int, int, int, int, int, int]:
    """Count cohort support and missingness for one graph segment.

    Returns the same six-tuple shape as
    :func:`~privy.io.vcf.extract_cohort_counts`:
    ``(target_support_n, target_total_n, offtarget_support_n,
    offtarget_total_n, target_missing_n, offtarget_missing_n)``

    - **support**: sample traverses this segment.
    - **missing**: sample has no path/walk coverage at the locus at all.
      Only detectable when coordinate information is available.
    - **absent** (not returned): sample is present at the locus but traverses
      a different segment.  Counted as "called but not supporting."

    Args:
        graph: Parsed :class:`GfaGraph`.
        seg_name: Segment to evaluate.
        seg_contig: Reference contig for missingness detection
            (``None`` → can't detect missing).
        seg_start: Segment start for missingness detection
            (``None`` → can't detect missing).
        seg_end: Segment end for missingness detection
            (``None`` → can't detect missing).
        target_samples: Full list of target sample names from the cohort
            definition (including samples absent from the GFA).
        offtarget_samples: Full list of off-target sample names.

    Returns:
        ``(ts_n, tt_n, os_n, ot_n, tm_n, om_n)``
    """
    traversing = get_samples_traversing_segment(graph, seg_name)

    present: set[str] | None = None
    if seg_contig is not None and seg_start is not None and seg_end is not None:
        present = get_samples_present_at_locus(graph, seg_contig, seg_start, seg_end)

    ts_n = 0
    tm_n = 0
    os_n = 0
    om_n = 0

    for sample in target_samples:
        if sample in traversing:
            ts_n += 1
        elif present is not None and sample not in present:
            tm_n += 1
        # else: present at locus, doesn't traverse this segment → absent, not missing

    for sample in offtarget_samples:
        if sample in traversing:
            os_n += 1
        elif present is not None and sample not in present:
            om_n += 1

    return ts_n, len(target_samples), os_n, len(offtarget_samples), tm_n, om_n


# ---------------------------------------------------------------------------
# GFA line parsers
# ---------------------------------------------------------------------------


def _parse_optional_tags(fields: list[str], start_idx: int) -> dict[str, Any]:
    """Parse ``TAG:TYPE:VALUE`` optional fields from *fields[start_idx:]*.

    Numeric types (``i``, ``f``) are converted; all others are kept as strings.
    Malformed tags are silently skipped.
    """
    tags: dict[str, Any] = {}
    for f in fields[start_idx:]:
        m = re.match(r"^([A-Za-z][A-Za-z0-9]):([AiZfJHB]):(.+)$", f)
        if not m:
            continue
        tag, type_code, value = m.group(1), m.group(2), m.group(3)
        if type_code == "i":
            try:
                tags[tag] = int(value)
            except ValueError:
                tags[tag] = value
        elif type_code == "f":
            try:
                tags[tag] = float(value)
            except ValueError:
                tags[tag] = value
        else:
            tags[tag] = value
    return tags


def _parse_h_line(fields: list[str], header_tags: dict[str, str]) -> None:
    tags = _parse_optional_tags(fields, 1)
    for k, v in tags.items():
        header_tags[k] = str(v)


def _parse_s_line(fields: list[str], line_num: int) -> GfaSegment:
    if len(fields) < 3:
        raise ValueError(
            f"Line {line_num}: S-line requires at least 3 tab-separated fields, "
            f"got {len(fields)}: {fields!r}"
        )
    name = fields[1]
    sequence = fields[2]
    tags = _parse_optional_tags(fields, 3)

    if sequence != "*":
        length = len(sequence)
    elif "LN" in tags:
        length = int(tags["LN"])
    else:
        length = 0

    ref_contig = str(tags["SN"]) if "SN" in tags else None
    ref_start: int | None = int(tags["SO"]) if "SO" in tags else None
    ref_end: int | None = (ref_start + length) if (ref_start is not None and length > 0) else None

    return GfaSegment(
        name=name,
        sequence=sequence,
        length=length,
        tags=tags,
        ref_contig=ref_contig,
        ref_start=ref_start,
        ref_end=ref_end,
    )


def _parse_l_line(fields: list[str], line_num: int) -> GfaLink:
    if len(fields) < 6:
        raise ValueError(
            f"Line {line_num}: L-line requires at least 6 fields, "
            f"got {len(fields)}: {fields!r}"
        )
    return GfaLink(
        from_seg=fields[1],
        from_orient=fields[2],
        to_seg=fields[3],
        to_orient=fields[4],
        overlap=fields[5],
        tags=_parse_optional_tags(fields, 6),
    )


def _parse_p_line(fields: list[str], line_num: int) -> GfaPath:
    if len(fields) < 4:
        raise ValueError(
            f"Line {line_num}: P-line requires at least 4 fields, "
            f"got {len(fields)}: {fields!r}"
        )
    path_name = fields[1]
    sample, haplotype = _parse_sample_from_path_name(path_name)

    segment_names: list[str] = []
    orientations: list[str] = []
    for so in fields[2].split(","):
        so = so.strip()
        if not so:
            continue
        if so.endswith("+") or so.endswith("-"):
            segment_names.append(so[:-1])
            orientations.append(so[-1])
        else:
            segment_names.append(so)
            orientations.append("+")

    return GfaPath(
        name=path_name,
        sample=sample,
        haplotype=haplotype,
        segment_names=segment_names,
        orientations=orientations,
        tags=_parse_optional_tags(fields, 4),
    )


def _parse_w_line(fields: list[str], line_num: int) -> GfaWalk:
    """Parse one W-line.

    Format: ``W <sample_id> <hap_index> <seq_id> <seq_start> <seq_end> <walk>``

    Walk field uses ``>`` (forward) and ``<`` (reverse) before each segment name.
    """
    if len(fields) < 7:
        raise ValueError(
            f"Line {line_num}: W-line requires at least 7 fields, "
            f"got {len(fields)}: {fields!r}"
        )
    try:
        hap_index = int(fields[2])
        seq_start = int(fields[4])
        seq_end = int(fields[5])
    except ValueError as exc:
        raise ValueError(
            f"Line {line_num}: W-line numeric field parse error: {exc}"
        ) from exc

    steps: list[GfaWalkStep] = []
    for m in re.finditer(r"([><])([^><]+)", fields[6]):
        orient = "+" if m.group(1) == ">" else "-"
        steps.append(GfaWalkStep(segment=m.group(2), orient=orient))

    return GfaWalk(
        sample=fields[1],
        hap_index=hap_index,
        seq_id=fields[3],
        seq_start=seq_start,
        seq_end=seq_end,
        steps=steps,
        tags=_parse_optional_tags(fields, 7),
    )


def _parse_sample_from_path_name(path_name: str) -> tuple[str, int | None]:
    """Extract (sample_name, haplotype_index) from a P-line path name.

    Supports the pangenome convention ``SAMPLE#HAP_INDEX#CONTIG``
    (e.g. ``HG002#1#chr1``) and plain names (e.g. ``HG002``).
    """
    parts = path_name.split("#")
    if len(parts) >= 2:
        sample = parts[0]
        try:
            haplotype: int | None = int(parts[1])
        except ValueError:
            haplotype = None
        return sample, haplotype
    return path_name, None


# ---------------------------------------------------------------------------
# Index building
# ---------------------------------------------------------------------------


def _build_indices(graph: GfaGraph) -> None:
    """Build all inverted indices and the position index on *graph* in-place."""
    # P-line indices
    for path_name, path in graph.paths.items():
        graph.sample_to_paths.setdefault(path.sample, []).append(path_name)
        for seg_name in path.segment_names:
            graph.segment_to_paths.setdefault(seg_name, []).append(path_name)

    # W-line indices
    for walk_idx, walk in enumerate(graph.walks):
        graph.sample_to_walks.setdefault(walk.sample, []).append(walk_idx)
        for step in walk.steps:
            graph.segment_to_walks.setdefault(step.segment, []).append(walk_idx)

    # Position index from SN/SO/LN tags
    raw: dict[str, list[tuple[int, int, str]]] = {}
    for seg_name, seg in graph.segments.items():
        if (
            seg.ref_contig is not None
            and seg.ref_start is not None
            and seg.ref_end is not None
        ):
            raw.setdefault(seg.ref_contig, []).append(
                (seg.ref_start, seg.ref_end, seg_name)
            )
    for contig, entries in raw.items():
        graph._contig_segments[contig] = sorted(entries)
