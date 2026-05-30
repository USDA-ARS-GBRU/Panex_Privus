"""Detect multi-allelic microhaplotypes from a pangenome graph.

Graph-native, snarl-free MVP: the invariant backbone is the set of *core* segments
present in (almost) every genome.  Between two consecutive core segments, each
genome's traversed interior is its local allele; hashing the orientation-aware
interior sequence (MD5) gives an allele id, and loci with more than one distinct
allele are emitted as microhaplotypes.

This deliberately mirrors the PHG reference-range / haplotype-id model and the
crop-microhaplotype framing (Zhao et al. 2026): linked variants within a short
locus, inherited as a unit, encoded as a multi-allelic marker.  Full snarl/bubble
decomposition (``vg deconstruct``) is a future ingestion path; this needs no
external tools.

All coordinates are 0-based half-open.
"""

from __future__ import annotations

import hashlib

from privy.io.gfa import GfaGraph
from privy.microhap.model import Microhaplotype
from privy.synteny.coordinates import PathCoordinateModel

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def _segment_token(graph: GfaGraph, segment: str, orientation: str) -> str:
    """Orientation-aware sequence (or a stable ``<name>`` token if bases absent)."""
    seg = graph.segments.get(segment)
    if seg is None or seg.sequence in ("", "*"):
        return f"<{segment}{orientation}>"
    return seg.sequence if orientation == "+" else _revcomp(seg.sequence)


def detect_microhaplotypes(
    graph: GfaGraph,
    model: PathCoordinateModel,
    ref_path: str,
    *,
    min_core_fraction: float = 1.0,
    multiallelic_only: bool = True,
) -> list[Microhaplotype]:
    """Detect microhaplotype loci between consecutive core segments.

    Args:
        graph: Parsed GFA (for segment sequences).
        model: Coordinate model over the same graph.
        ref_path: Path id whose coordinates anchor the loci (and defines the
            reference allele).
        min_core_fraction: A segment is "core" when present on at least this
            fraction of paths (1.0 = present in every genome).
        multiallelic_only: Emit only loci with >1 distinct allele (default).

    Returns:
        Microhaplotypes sorted by reference start.
    """
    if ref_path not in model:
        raise KeyError(f"unknown reference path {ref_path!r}")

    paths = model.path_ids()
    n_paths = len(paths)
    threshold = max(1, int(round(min_core_fraction * n_paths)))

    # Core segments: present on >= threshold paths AND on the reference.
    core: list[str] = []
    seen: set[str] = set()
    for step in model.iter_steps(ref_path):
        seg = step.segment
        if seg in seen:
            continue
        seen.add(seg)
        present = sum(1 for p in paths if model.occurrences(p, seg))
        if present >= threshold:
            core.append(seg)   # already in reference order

    loci: list[Microhaplotype] = []
    for left, right in zip(core, core[1:], strict=False):
        locus = _build_locus(graph, model, ref_path, paths, left, right)
        if locus is None:
            continue
        if multiallelic_only and not locus.is_multiallelic:
            continue
        loci.append(locus)

    loci.sort(key=lambda m: (m.contig, m.start, m.end))
    return loci


def _build_locus(
    graph: GfaGraph,
    model: PathCoordinateModel,
    ref_path: str,
    paths: list[str],
    left: str,
    right: str,
) -> Microhaplotype | None:
    """Build one microhaplotype between core segments *left* and *right*."""
    ref_left = model.occurrences(ref_path, left)
    ref_right = model.occurrences(ref_path, right)
    if not ref_left or not ref_right:
        return None
    contig, start = model.to_stable(ref_path, ref_left[0].start)
    _, end = model.to_stable(ref_path, ref_right[0].end - 1)
    end += 1

    alleles: dict[str, str] = {}
    tokens: dict[str, str] = {}
    for path in paths:
        occ_l = model.occurrences(path, left)
        occ_r = model.occurrences(path, right)
        if not occ_l or not occ_r:
            continue  # genome missing a flank → no call (missing data)
        lo = occ_l[0].end
        hi = occ_r[0].start
        if hi < lo:
            continue  # flanks out of order on this genome (rearranged) → skip
        interior = [
            s for s in model.segments_in_range(path, lo, hi)
            if lo <= s.start and s.end <= hi
        ]
        token = "".join(_segment_token(graph, s.segment, s.orientation) for s in interior)
        allele_id = hashlib.md5(token.encode("utf-8")).hexdigest()  # noqa: S324
        alleles[path] = allele_id
        tokens.setdefault(allele_id, "".join(
            f"{s.segment}{s.orientation}" for s in interior
        ) or "(empty)")

    if not alleles:
        return None

    return Microhaplotype(
        locus_id=f"{contig}_{start:010d}",
        contig=contig,
        start=start,
        end=end,
        alleles=alleles,
        ref_allele=alleles.get(ref_path),
        allele_tokens=tokens,
    )
