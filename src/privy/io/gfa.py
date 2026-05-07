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

import gzip
import logging
import pickle
import re
import sqlite3
import time
from bisect import bisect_right
from collections import defaultdict
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, TextIO

log = logging.getLogger("privy.io.gfa")
_WALK_SEGMENT_RE = re.compile(r"[><]([^><\t\r\n]+)")
_SQLITE_HEADER = b"SQLite format 3\0"
_LEGACY_GFA_SCAN_INDEX_MAGIC = b"PRIVY_GFA_SCAN_INDEX\0"
_GFA_SCAN_INDEX_SCHEMA_VERSION = 3
_SQLITE_INSERT_BATCH_SIZE = 100_000


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


@dataclass(slots=True)
class GfaScanSegment:
    """Lightweight coordinate record for scan-time GFA segment evaluation."""

    contig: str
    start: int
    end: int
    length: int
    sample_mask: int = 0


@dataclass(frozen=True, slots=True)
class GfaPresenceIntervals:
    """Merged half-open intervals for one sample on one contig."""

    starts: tuple[int, ...]
    ends: tuple[int, ...]

    def overlaps(self, start: int, end: int) -> bool:
        """Return True when any stored interval overlaps ``[start, end)``."""
        idx = bisect_right(self.starts, start) - 1
        if idx >= 0 and self.ends[idx] > start:
            return True

        next_idx = idx + 1
        return next_idx < len(self.starts) and self.starts[next_idx] < end


@dataclass(slots=True)
class GfaScanIndex:
    """Memory-conscious GFA index for private-segment scans.

    Unlike :class:`GfaGraph`, this object does not retain sequences, links, full
    walks, or full paths.  It keeps only the information needed by the GFA scan:
    coordinate-tagged segments, cohort-sample traversal bitmasks, and compact
    per-sample coverage intervals for missingness checks.
    """

    segments: dict[str, GfaScanSegment]
    contig_segments: dict[str, list[tuple[int, int, str]]]
    segment_sample_mask: dict[str, int]
    sample_intervals: list[dict[str, GfaPresenceIntervals]]
    sample_order: tuple[str, ...]
    sample_to_index: dict[str, int]
    samples_seen: set[str]
    sqlite_index_path: Path | None = None
    sqlite_contigs: tuple[str, ...] = ()
    sqlite_contig_counts: dict[str, int] = field(default_factory=dict)
    sqlite_contig_ranges: dict[str, tuple[int, int]] = field(default_factory=dict)
    n_segments: int = 0
    n_links: int = 0
    n_paths: int = 0
    n_walks: int = 0
    metadata: dict[str, Any] = field(default_factory=dict)

    def sample_mask(self, samples: list[str] | tuple[str, ...]) -> int:
        """Return a bitmask for samples present in this index's sample order."""
        mask = 0
        for sample in samples:
            idx = self.sample_to_index.get(sample)
            if idx is not None:
                mask |= 1 << idx
        return mask

    def present_mask(self, contig: str, start: int, end: int) -> int:
        """Return a cohort-sample bitmask for samples with coverage at a locus."""
        mask = 0
        for idx, intervals_by_contig in enumerate(self.sample_intervals):
            intervals = intervals_by_contig.get(contig)
            if intervals is not None and intervals.overlaps(start, end):
                mask |= 1 << idx
        return mask

    def mask_to_statuses(
        self,
        support_mask: int,
        present_mask: int,
        samples: list[str] | tuple[str, ...],
    ) -> dict[str, str]:
        """Return sample → ``traverses``/``absent``/``missing`` statuses."""
        statuses: dict[str, str] = {}
        for sample in samples:
            idx = self.sample_to_index.get(sample)
            if idx is None:
                statuses[sample] = "missing"
                continue

            bit = 1 << idx
            if support_mask & bit:
                statuses[sample] = "traverses"
            elif present_mask & bit:
                statuses[sample] = "absent"
            else:
                statuses[sample] = "missing"
        return statuses

    def coordinate_segment_count(self) -> int:
        """Return the number of coordinate-tagged segments available to scan."""
        if self.sqlite_contig_counts:
            return sum(self.sqlite_contig_counts.values())
        return sum(len(entries) for entries in self.contig_segments.values())

    def contig_names(self) -> tuple[str, ...]:
        """Return contigs with coordinate-tagged segments."""
        if self.sqlite_contigs:
            return self.sqlite_contigs
        return tuple(self.contig_segments.keys())

    def has_contig(self, contig: str) -> bool:
        """Return True when *contig* has coordinate-tagged segments."""
        if self.sqlite_contig_counts:
            return contig in self.sqlite_contig_counts
        return contig in self.contig_segments

    def contig_segment_count(self, contig: str) -> int:
        """Return the number of coordinate-tagged segments for one contig."""
        if self.sqlite_contig_counts:
            return self.sqlite_contig_counts.get(contig, 0)
        return len(self.contig_segments.get(contig, []))

    def iter_contig_segments(self, contig: str) -> Iterator[tuple[int, int, str, int]]:
        """Yield ``(start, end, segment_name, sample_mask)`` rows for *contig*."""
        if self.sqlite_index_path is not None:
            row_range = self.sqlite_contig_ranges.get(contig)
            if row_range is None:
                return
            first_segment_id, last_segment_id = row_range
            with _connect_gfa_scan_index_db(self.sqlite_index_path, readonly=True) as conn:
                rows = conn.execute(
                    "SELECT start, end, name, sample_mask FROM segments "
                    "WHERE segment_id BETWEEN ? AND ? ORDER BY segment_id",
                    (first_segment_id, last_segment_id),
                )
                for start, end, name, mask_blob in rows:
                    yield int(start), int(end), str(name), _decode_sample_mask(mask_blob)
            return

        for entry in self.contig_segments.get(contig, []):
            start, end, seg_name = entry
            yield start, end, seg_name, self.segment_sample_mask.get(seg_name, 0)

    def support_mask_for_segment(self, segment_name: str) -> int:
        """Return the traversal mask for one segment name."""
        mask = self.segment_sample_mask.get(segment_name)
        if mask is not None:
            return mask
        if self.sqlite_index_path is None:
            return 0
        with _connect_gfa_scan_index_db(self.sqlite_index_path, readonly=True) as conn:
            row = conn.execute(
                "SELECT sample_mask FROM segments WHERE name = ? LIMIT 1",
                (segment_name,),
            ).fetchone()
        if row is None:
            return 0
        return _decode_sample_mask(row[0])


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def parse_gfa(gfa_path: Path) -> GfaGraph:
    """Parse a GFA1 or GFA1.1 file and return a fully indexed :class:`GfaGraph`.

    Reads S, L, P, W, and H lines.  Unknown line types (J, C, E, F) are
    silently skipped.  Plain-text ``.gfa`` and gzip-compressed ``.gfa.gz``
    inputs are supported.

    Args:
        gfa_path: Path to a GFA or GFA.GZ file.

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

    with _open_gfa_text(gfa_path) as fh:
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


def _open_gfa_text(gfa_path: Path) -> TextIO:
    """Open a plain or gzip-compressed GFA path as text."""
    if gfa_path.suffix.lower() == ".gz":
        return gzip.open(gfa_path, "rt", encoding="utf-8")
    return open(gfa_path, encoding="utf-8")


def get_gfa_samples(gfa_path: Path) -> list[str]:
    """Return a sorted list of sample names present in a GFA file.

    Sample names are extracted from:
    - P-line path names (``SAMPLE#hap#contig`` or plain ``SAMPLE``).
    - W-line ``sample_id`` fields.

    Args:
        gfa_path: Path to a GFA1/1.1 file, plain-text or gzip-compressed.

    Returns:
        Sorted list of unique sample name strings.
    """
    graph = parse_gfa(gfa_path)
    return sorted(set(graph.sample_to_paths) | set(graph.sample_to_walks))


def default_gfa_index_path(gfa_path: Path) -> Path:
    """Return the sidecar Privy GFA index path for *gfa_path*."""
    return Path(f"{gfa_path}.privy.gfaidx")


def write_gfa_scan_index(
    scan_index: GfaScanIndex,
    index_path: Path,
    gfa_path: Path,
) -> None:
    """Write a reusable Privy GFA scan index.

    The index is a SQLite sidecar storing exactly the scan structures Privy
    needs. Segment rows are streamed by contig during ``privy scan`` so large
    indexes do not need to be deserialized into memory all at once.
    """
    index_path.parent.mkdir(parents=True, exist_ok=True)
    metadata = {
        "schema_version": _GFA_SCAN_INDEX_SCHEMA_VERSION,
        "source": _gfa_file_fingerprint(gfa_path),
        "sample_count": len(scan_index.sample_order),
        "coordinate_segments": scan_index.coordinate_segment_count(),
        "contigs": len(scan_index.contig_names()),
        "n_segments": scan_index.n_segments,
        "n_links": scan_index.n_links,
        "n_paths": scan_index.n_paths,
        "n_walks": scan_index.n_walks,
        "created_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
    }
    scan_index.metadata = metadata
    if index_path.exists():
        index_path.unlink()

    log.info("Writing SQLite GFA scan index: %s", index_path)
    with _connect_gfa_scan_index_db(index_path, readonly=False) as conn:
        _init_gfa_scan_index_db(conn)
        conn.execute(
            "INSERT INTO metadata(key, value) VALUES (?, ?)",
            ("metadata", sqlite3.Binary(_pickle_payload(metadata))),
        )
        conn.executemany(
            "INSERT INTO samples(sample_idx, sample) VALUES (?, ?)",
            enumerate(scan_index.sample_order),
        )
        conn.executemany(
            "INSERT INTO samples_seen(sample) VALUES (?)",
            ((sample,) for sample in sorted(scan_index.samples_seen)),
        )

        interval_rows: list[tuple[int, str, sqlite3.Binary]] = []
        for sample_idx, intervals_by_contig in enumerate(scan_index.sample_intervals):
            for contig, intervals in intervals_by_contig.items():
                interval_rows.append((
                    sample_idx,
                    contig,
                    sqlite3.Binary(_pickle_payload(intervals)),
                ))
        conn.executemany(
            "INSERT INTO sample_intervals(sample_idx, contig, payload) VALUES (?, ?, ?)",
            interval_rows,
        )

        next_segment_id = 1
        written_segments = 0
        for contig, entries in scan_index.contig_segments.items():
            first_segment_id = next_segment_id
            batch: list[tuple[int, str, int, int, str, sqlite3.Binary]] = []
            contig_written = 0
            for start, end, seg_name in entries:
                batch.append((
                    next_segment_id,
                    contig,
                    start,
                    end,
                    seg_name,
                    sqlite3.Binary(_encode_sample_mask(
                        scan_index.segment_sample_mask.get(seg_name, 0)
                    )),
                ))
                next_segment_id += 1
                if len(batch) >= _SQLITE_INSERT_BATCH_SIZE:
                    batch_size = len(batch)
                    _insert_segment_batch(conn, batch)
                    contig_written += batch_size
                    written_segments += batch_size
                    batch.clear()
            if batch:
                batch_size = len(batch)
                _insert_segment_batch(conn, batch)
                contig_written += batch_size
                written_segments += batch_size

            last_segment_id = next_segment_id - 1
            conn.execute(
                "INSERT INTO contigs(contig, n_segments, first_segment_id, "
                "last_segment_id) VALUES (?, ?, ?, ?)",
                (contig, len(entries), first_segment_id, last_segment_id),
            )
            log.info(
                "  indexed %s: %d segments (%d total written)",
                contig,
                contig_written,
                written_segments,
            )

        conn.commit()


def load_gfa_scan_index(
    index_path: Path,
    gfa_path: Path | None = None,
) -> GfaScanIndex:
    """Load a reusable Privy GFA scan index and optionally validate its source."""
    if not index_path.exists():
        raise FileNotFoundError(f"GFA index file not found: {index_path}")

    with open(index_path, "rb") as fh:
        header = fh.read(max(len(_SQLITE_HEADER), len(_LEGACY_GFA_SCAN_INDEX_MAGIC)))
        if header.startswith(_LEGACY_GFA_SCAN_INDEX_MAGIC):
            raise ValueError(
                f"{index_path} is a legacy pickle-based Privy GFA index. "
                "Rebuild it with 'privy index gfa --gfa <GFA> --force' to "
                "create the SQLite-backed index used by this Privy version."
            )
        if not header.startswith(_SQLITE_HEADER):
            raise ValueError(
                f"{index_path} is not a Privy GFA scan index. "
                "Rebuild it with 'privy index gfa'."
            )

    with _connect_gfa_scan_index_db(index_path, readonly=True) as conn:
        metadata_row = conn.execute(
            "SELECT value FROM metadata WHERE key = 'metadata'"
        ).fetchone()
        if metadata_row is None:
            raise ValueError(f"{index_path} has an invalid Privy GFA index payload.")
        metadata = _unpickle_payload(metadata_row[0])

        if not isinstance(metadata, dict):
            raise ValueError(f"{index_path} has an invalid Privy GFA index payload.")

        schema_version = metadata.get("schema_version")
        if schema_version != _GFA_SCAN_INDEX_SCHEMA_VERSION:
            raise ValueError(
                f"{index_path} uses GFA index schema {schema_version!r}; "
                f"this Privy version expects {_GFA_SCAN_INDEX_SCHEMA_VERSION}. "
                "Rebuild it with 'privy index gfa --gfa <GFA> --force'."
            )

        sample_rows = conn.execute(
            "SELECT sample_idx, sample FROM samples ORDER BY sample_idx"
        ).fetchall()
        sample_order = tuple(str(row[1]) for row in sample_rows)
        sample_to_index = {sample: idx for idx, sample in enumerate(sample_order)}
        samples_seen = {
            str(row[0]) for row in conn.execute("SELECT sample FROM samples_seen")
        }
        sample_intervals: list[dict[str, GfaPresenceIntervals]] = [
            {} for _ in sample_order
        ]
        for sample_idx, contig, payload in conn.execute(
            "SELECT sample_idx, contig, payload FROM sample_intervals"
        ):
            intervals = _unpickle_payload(payload)
            if not isinstance(intervals, GfaPresenceIntervals):
                raise ValueError(
                    f"{index_path} has invalid sample interval data for {contig!r}."
                )
            sample_intervals[int(sample_idx)][str(contig)] = intervals

        contig_counts: dict[str, int] = {}
        contig_ranges: dict[str, tuple[int, int]] = {}
        for contig, n_segments, first_segment_id, last_segment_id in conn.execute(
            "SELECT contig, n_segments, first_segment_id, last_segment_id FROM contigs"
        ):
            contig_name = str(contig)
            contig_counts[contig_name] = int(n_segments)
            contig_ranges[contig_name] = (int(first_segment_id), int(last_segment_id))

    if gfa_path is not None:
        _validate_gfa_index_source(metadata, gfa_path, index_path)

    return GfaScanIndex(
        segments={},
        contig_segments={},
        segment_sample_mask={},
        sample_intervals=sample_intervals,
        sample_order=sample_order,
        sample_to_index=sample_to_index,
        samples_seen=samples_seen,
        sqlite_index_path=index_path,
        sqlite_contigs=tuple(sorted(contig_counts)),
        sqlite_contig_counts=contig_counts,
        sqlite_contig_ranges=contig_ranges,
        n_segments=int(metadata.get("n_segments", 0) or 0),
        n_links=int(metadata.get("n_links", 0) or 0),
        n_paths=int(metadata.get("n_paths", 0) or 0),
        n_walks=int(metadata.get("n_walks", 0) or 0),
        metadata=metadata,
    )


def _connect_gfa_scan_index_db(index_path: Path, *, readonly: bool) -> sqlite3.Connection:
    if readonly:
        conn = sqlite3.connect(f"file:{index_path}?mode=ro", uri=True)
        conn.execute("PRAGMA query_only=ON")
        conn.execute("PRAGMA temp_store=MEMORY")
        conn.execute("PRAGMA cache_size=-200000")
        conn.execute("PRAGMA mmap_size=268435456")
        return conn

    conn = sqlite3.connect(index_path)
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("PRAGMA synchronous=OFF")
    conn.execute("PRAGMA temp_store=MEMORY")
    conn.execute("PRAGMA cache_size=-200000")
    return conn


def _init_gfa_scan_index_db(conn: sqlite3.Connection) -> None:
    conn.executescript(
        """
        CREATE TABLE metadata (
            key TEXT PRIMARY KEY,
            value BLOB NOT NULL
        );
        CREATE TABLE samples (
            sample_idx INTEGER PRIMARY KEY,
            sample TEXT NOT NULL
        );
        CREATE TABLE samples_seen (
            sample TEXT PRIMARY KEY
        );
        CREATE TABLE sample_intervals (
            sample_idx INTEGER NOT NULL,
            contig TEXT NOT NULL,
            payload BLOB NOT NULL,
            PRIMARY KEY (sample_idx, contig)
        );
        CREATE TABLE contigs (
            contig TEXT PRIMARY KEY,
            n_segments INTEGER NOT NULL,
            first_segment_id INTEGER NOT NULL,
            last_segment_id INTEGER NOT NULL
        );
        CREATE TABLE segments (
            segment_id INTEGER PRIMARY KEY,
            contig TEXT NOT NULL,
            start INTEGER NOT NULL,
            end INTEGER NOT NULL,
            name TEXT NOT NULL,
            sample_mask BLOB NOT NULL
        );
        """
    )


def _insert_segment_batch(
    conn: sqlite3.Connection,
    batch: list[tuple[int, str, int, int, str, sqlite3.Binary]],
) -> None:
    conn.executemany(
        "INSERT INTO segments(segment_id, contig, start, end, name, sample_mask) "
        "VALUES (?, ?, ?, ?, ?, ?)",
        batch,
    )


def _encode_sample_mask(mask: int) -> bytes:
    n_bytes = max(1, (mask.bit_length() + 7) // 8)
    return mask.to_bytes(n_bytes, "little")


def _decode_sample_mask(mask_blob: bytes) -> int:
    return int.from_bytes(mask_blob, "little")


def _pickle_payload(payload: Any) -> bytes:
    return pickle.dumps(payload, protocol=pickle.HIGHEST_PROTOCOL)


def _unpickle_payload(payload: bytes) -> Any:
    return pickle.loads(payload)  # noqa: S301 - trusted local cache file


def _gfa_file_fingerprint(gfa_path: Path) -> dict[str, Any]:
    """Return stable-enough source metadata for stale-index detection."""
    stat = gfa_path.stat()
    return {
        "path": str(gfa_path),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }


def _validate_gfa_index_source(
    metadata: dict[str, Any],
    gfa_path: Path,
    index_path: Path,
) -> None:
    """Raise if *index_path* does not appear to match *gfa_path*."""
    expected = metadata.get("source")
    if not isinstance(expected, dict):
        raise ValueError(
            f"{index_path} does not contain source GFA metadata. "
            "Rebuild it with 'privy index gfa'."
        )

    observed = _gfa_file_fingerprint(gfa_path)
    if (
        expected.get("size") != observed["size"]
        or expected.get("mtime_ns") != observed["mtime_ns"]
    ):
        raise ValueError(
            f"{index_path} does not match {gfa_path}. "
            "The GFA file size or modification time changed; rebuild the index "
            "with 'privy index gfa'."
        )


def build_gfa_scan_index(
    gfa_path: Path,
    sample_names: list[str] | tuple[str, ...] | None,
    filter_contig: str | None = None,
    filter_start: int | None = None,
    filter_end: int | None = None,
    progress_interval_seconds: float = 30.0,
) -> GfaScanIndex:
    """Build a memory-conscious index for GFA scanning.

    This performs a single streaming pass over the GFA. It records
    coordinate-tagged segments, cohort sample traversal bitmasks, and compact
    coverage intervals needed for missingness checks.
    """
    if not gfa_path.exists():
        raise FileNotFoundError(f"GFA file not found: {gfa_path}")

    sample_order_list = list(dict.fromkeys(sample_names or ()))
    sample_to_idx = {sample: idx for idx, sample in enumerate(sample_order_list)}
    raw_intervals: list[defaultdict[str, list[tuple[int, int]]]] = [
        defaultdict(list) for _ in sample_order_list
    ]
    index_all_samples = sample_names is None
    segments: dict[str, GfaScanSegment] = {}
    contig_segments: dict[str, list[tuple[int, int, str]]] = {}
    contig_pool: dict[str, str] = {}
    samples_seen: set[str] = set()
    n_segments = 0
    n_links = 0
    n_paths = 0
    n_walks = 0
    n_cohort_paths = 0
    n_cohort_walks = 0
    n_cohort_segment_refs = 0
    n_p_segment_refs_without_coords = 0
    header_tags: dict[str, str] = {}
    record_count = 0
    uncompressed_bytes = 0
    started_at = time.monotonic()
    last_progress_at = started_at

    def maybe_log_progress(context: str = "indexing") -> None:
        nonlocal last_progress_at
        if progress_interval_seconds <= 0:
            return
        now = time.monotonic()
        if now - last_progress_at < progress_interval_seconds:
            return
        last_progress_at = now
        log.info(
            "GFA %s progress | records=%d | read=%.2f GB | "
            "coordinate_segments=%d | paths=%d | walks=%d | "
            "cohort_paths=%d | cohort_walks=%d | cohort_segment_refs=%d | "
            "samples_seen=%d | elapsed=%.1fs",
            context,
            record_count,
            uncompressed_bytes / 1_000_000_000,
            len(segments),
            n_paths,
            n_walks,
            n_cohort_paths,
            n_cohort_walks,
            n_cohort_segment_refs,
            len(samples_seen),
            now - started_at,
        )

    def sample_index_for(sample: str) -> int | None:
        idx = sample_to_idx.get(sample)
        if idx is not None or not index_all_samples:
            return idx
        idx = len(sample_order_list)
        sample_order_list.append(sample)
        sample_to_idx[sample] = idx
        raw_intervals.append(defaultdict(list))
        return idx

    log.info("Building streaming GFA scan index: %s", gfa_path)

    pending_segment_sample_mask: dict[str, int] = {}

    with _open_gfa_text(gfa_path) as fh:
        for line_num, raw_line in enumerate(fh, 1):
            if not raw_line or raw_line.startswith("#"):
                continue
            record_count += 1
            uncompressed_bytes += len(raw_line)
            rec_type = raw_line[0]
            if rec_type == "H":
                _parse_h_line(raw_line.rstrip("\n").split("\t"), header_tags)
                if header_tags.get("VN", "").startswith("2"):
                    raise ValueError(
                        f"GFA version 2 is not supported (found VN:{header_tags['VN']}). "
                        "Only GFA1 and GFA1.1 are supported."
                    )
            elif rec_type == "S":
                n_segments += 1
                parsed = _parse_scan_s_line(raw_line, line_num, contig_pool)
                if parsed is None:
                    continue
                seg_name, segment_record = parsed
                if filter_contig is not None and segment_record.contig != filter_contig:
                    pending_segment_sample_mask.pop(seg_name, None)
                    continue
                if filter_start is not None and segment_record.end <= filter_start:
                    pending_segment_sample_mask.pop(seg_name, None)
                    continue
                if filter_end is not None and segment_record.start >= filter_end:
                    pending_segment_sample_mask.pop(seg_name, None)
                    continue
                segment_record.sample_mask = pending_segment_sample_mask.pop(
                    seg_name, 0
                )
                segments[seg_name] = segment_record
                contig_segments.setdefault(segment_record.contig, []).append(
                    (segment_record.start, segment_record.end, seg_name)
                )
            elif rec_type == "L":
                n_links += 1
            elif rec_type == "W":
                n_walks += 1
                tabs = _required_tab_positions(raw_line, 6, line_num, "W")
                sample = raw_line[tabs[0] + 1:tabs[1]]
                samples_seen.add(sample)
                sample_idx = sample_index_for(sample)
                if sample_idx is None:
                    maybe_log_progress()
                    continue

                n_cohort_walks += 1
                try:
                    walk_contig = _intern_text(
                        raw_line[tabs[2] + 1:tabs[3]], contig_pool
                    )
                    walk_start = int(raw_line[tabs[3] + 1:tabs[4]])
                    walk_end = int(raw_line[tabs[4] + 1:tabs[5]])
                except ValueError as exc:
                    raise ValueError(
                        f"Line {line_num}: W-line numeric field parse error: {exc}"
                    ) from exc

                if filter_contig is not None and walk_contig != filter_contig:
                    maybe_log_progress()
                    continue
                if filter_start is not None and walk_end <= filter_start:
                    maybe_log_progress()
                    continue
                if filter_end is not None and walk_start >= filter_end:
                    maybe_log_progress()
                    continue

                raw_intervals[sample_idx][walk_contig].append((walk_start, walk_end))

                bit = 1 << sample_idx
                walk_field_start = tabs[5] + 1
                walk_field_end = _next_field_end(raw_line, walk_field_start)
                segment_refs_in_record = 0
                for seg_name in _iter_walk_segment_names(
                    raw_line, walk_field_start, walk_field_end
                ):
                    segment_refs_in_record += 1
                    n_cohort_segment_refs += 1
                    segment = segments.get(seg_name)
                    if segment is not None:
                        segment.sample_mask |= bit
                    else:
                        pending_segment_sample_mask[seg_name] = (
                            pending_segment_sample_mask.get(seg_name, 0) | bit
                        )
                    if segment_refs_in_record % 500_000 == 0:
                        maybe_log_progress("cohort walk")
                maybe_log_progress()

            elif rec_type == "P":
                n_paths += 1
                tabs = _required_tab_positions(raw_line, 3, line_num, "P")
                path_name = raw_line[tabs[0] + 1:tabs[1]]
                sample, _haplotype = _parse_sample_from_path_name(path_name)
                samples_seen.add(sample)
                sample_idx = sample_index_for(sample)
                if sample_idx is None:
                    maybe_log_progress()
                    continue

                n_cohort_paths += 1
                bit = 1 << sample_idx
                segment_refs_in_record = 0
                for seg_name in _iter_path_segment_names(
                    raw_line, tabs[1] + 1, tabs[2]
                ):
                    segment_refs_in_record += 1
                    n_cohort_segment_refs += 1
                    segment = segments.get(seg_name)
                    if segment is None:
                        pending_segment_sample_mask[seg_name] = (
                            pending_segment_sample_mask.get(seg_name, 0) | bit
                        )
                        n_p_segment_refs_without_coords += 1
                        if segment_refs_in_record % 500_000 == 0:
                            maybe_log_progress("cohort path")
                        continue
                    segment.sample_mask |= bit
                    raw_intervals[sample_idx][segment.contig].append(
                        (segment.start, segment.end)
                    )
                    if segment_refs_in_record % 500_000 == 0:
                        maybe_log_progress("cohort path")
                maybe_log_progress()

            if record_count % 100_000 == 0:
                maybe_log_progress()

    if header_tags.get("VN", "").startswith("2"):
        raise ValueError(
            f"GFA version 2 is not supported (found VN:{header_tags['VN']}). "
            "Only GFA1 and GFA1.1 are supported."
        )

    log.info(
        "Sorting GFA scan index | coordinate_segments=%d | contigs=%d",
        len(segments),
        len(contig_segments),
    )
    for entries in contig_segments.values():
        entries.sort()

    if n_p_segment_refs_without_coords:
        log.warning(
            "Skipped %d cohort P-line segment references while building presence "
            "intervals because their coordinates were unavailable. This usually "
            "means P-lines appeared before matching S-lines or those segments lack "
            "SN/SO/LN tags.",
            n_p_segment_refs_without_coords,
        )

    log.info("Merging GFA sample presence intervals")
    sample_intervals = [
        _merge_sample_intervals(intervals_by_contig)
        for intervals_by_contig in raw_intervals
    ]
    segment_sample_mask = {
        name: segment.sample_mask
        for name, segment in segments.items()
        if segment.sample_mask
    }

    log.info(
        "Built GFA scan index: %d coordinate segments, %d paths, %d walks, "
        "%d cohort segment references, %d samples, %.1fs elapsed",
        len(segments),
        n_paths,
        n_walks,
        n_cohort_segment_refs,
        len(samples_seen),
        time.monotonic() - started_at,
    )

    return GfaScanIndex(
        segments=segments,
        contig_segments=contig_segments,
        segment_sample_mask=segment_sample_mask,
        sample_intervals=sample_intervals,
        sample_order=tuple(sample_order_list),
        sample_to_index=sample_to_idx,
        samples_seen=samples_seen,
        n_segments=n_segments,
        n_links=n_links,
        n_paths=n_paths,
        n_walks=n_walks,
    )


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


def _parse_scan_s_line(
    raw_line: str,
    line_num: int,
    contig_pool: dict[str, str],
) -> tuple[str, GfaScanSegment] | None:
    """Parse the coordinate fields needed for scan-time segment indexing."""
    first_tab = raw_line.find("\t")
    second_tab = raw_line.find("\t", first_tab + 1)
    if first_tab < 0 or second_tab < 0:
        raise ValueError(
            f"Line {line_num}: S-line requires at least 3 tab-separated fields, "
            "but fewer were found"
        )
    third_tab = raw_line.find("\t", second_tab + 1)

    line_end = _line_text_end(raw_line)
    name = raw_line[first_tab + 1:second_tab]
    sequence_start = second_tab + 1
    sequence_end = third_tab if third_tab >= 0 else line_end
    tag_text = raw_line[third_tab + 1:line_end] if third_tab >= 0 else ""
    sn, so, ln = _parse_scan_coordinate_tags(tag_text)

    if sn is None or so is None:
        return None

    if ln is not None:
        length = ln
    elif raw_line[sequence_start:sequence_end] != "*":
        length = sequence_end - sequence_start
    else:
        return None

    return name, GfaScanSegment(
        contig=_intern_text(sn, contig_pool),
        start=so,
        end=so + length,
        length=length,
    )


def _parse_scan_coordinate_tags(tag_text: str) -> tuple[str | None, int | None, int | None]:
    """Parse only the S-line coordinate tags needed by the scan index."""
    sn: str | None = None
    so: int | None = None
    ln: int | None = None
    field_start = 0
    text_len = len(tag_text)

    while field_start < text_len:
        field_end = tag_text.find("\t", field_start)
        if field_end < 0:
            field_end = text_len

        if tag_text.startswith("SN:Z:", field_start, field_end):
            sn = tag_text[field_start + 5:field_end]
        elif tag_text.startswith("SO:i:", field_start, field_end):
            try:
                so = int(tag_text[field_start + 5:field_end])
            except ValueError:
                so = None
        elif tag_text.startswith("LN:i:", field_start, field_end):
            try:
                ln = int(tag_text[field_start + 5:field_end])
            except ValueError:
                ln = None

        field_start = field_end + 1

    return sn, so, ln


def _iter_walk_segment_names(
    walk: str,
    start_idx: int = 0,
    end_idx: int | None = None,
) -> Iterator[str]:
    """Yield segment names from a W-line walk string without building a list."""
    if end_idx is None:
        end_idx = len(walk)
    for match in _WALK_SEGMENT_RE.finditer(walk, start_idx, end_idx):
        yield match.group(1)


def _iter_path_segment_names(
    path_segments: str,
    start_idx: int = 0,
    end_idx: int | None = None,
) -> Iterator[str]:
    """Yield segment names from a P-line segment list without building a list."""
    if end_idx is None:
        end_idx = len(path_segments)
    start = start_idx
    for idx in range(start_idx, end_idx):
        if path_segments[idx] != ",":
            continue
        if start < idx:
            yield _strip_path_orientation(path_segments[start:idx])
        start = idx + 1
    if start < end_idx:
        yield _strip_path_orientation(path_segments[start:end_idx])


def _strip_path_orientation(segment: str) -> str:
    segment = segment.strip()
    if segment.endswith("+") or segment.endswith("-"):
        return segment[:-1]
    return segment


def _merge_sample_intervals(
    raw: dict[str, list[tuple[int, int]]],
) -> dict[str, GfaPresenceIntervals]:
    """Sort and merge raw per-contig intervals for one sample."""
    merged_by_contig: dict[str, GfaPresenceIntervals] = {}
    for contig, intervals in raw.items():
        if not intervals:
            continue

        merged: list[tuple[int, int]] = []
        for start, end in sorted(intervals):
            if end <= start:
                continue
            if not merged or start > merged[-1][1]:
                merged.append((start, end))
            else:
                prev_start, prev_end = merged[-1]
                merged[-1] = (prev_start, max(prev_end, end))

        if merged:
            starts, ends = zip(*merged, strict=False)
            merged_by_contig[contig] = GfaPresenceIntervals(
                starts=tuple(starts),
                ends=tuple(ends),
            )
    return merged_by_contig


def _required_tab_positions(
    raw_line: str,
    n_tabs: int,
    line_num: int,
    record_type: str,
) -> list[int]:
    """Return positions for the first ``n_tabs`` tabs or raise a parse error."""
    tabs: list[int] = []
    pos = -1
    for _ in range(n_tabs):
        pos = raw_line.find("\t", pos + 1)
        if pos < 0:
            raise ValueError(
                f"Line {line_num}: {record_type}-line requires at least "
                f"{n_tabs + 1} tab-separated fields"
            )
        tabs.append(pos)
    return tabs


def _next_field_end(raw_line: str, field_start: int) -> int:
    """Return the exclusive end index of a tab-delimited field."""
    tab = raw_line.find("\t", field_start)
    if tab >= 0:
        return tab
    return _line_text_end(raw_line)


def _line_text_end(raw_line: str) -> int:
    """Return line length excluding trailing newline characters."""
    end = len(raw_line)
    while end > 0 and raw_line[end - 1] in {"\n", "\r"}:
        end -= 1
    return end


def _intern_text(value: str, pool: dict[str, str]) -> str:
    """Reuse repeated small strings, mainly contig names, in large GFA indices."""
    existing = pool.get(value)
    if existing is not None:
        return existing
    pool[value] = value
    return value


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
