"""Multi-genome synteny orchestration, region grouping, and private-region tagging.

Builds typed synteny blocks for every query genome against a chosen reference
(:func:`build_synteny`), groups overlapping blocks into reference-anchored regions
(:func:`group_regions`), and — tying the comparative layer back to Privy's core —
flags regions whose structure is **target-private**: present in the target cohort
and absent from the off-target cohort (:func:`tag_region_privacy`).

Region presence is measured at the segment level (does a genome traverse the
region's reference segments?), so a deletion is detected even when a genome's
spanning block overlaps the region.

All coordinates are 0-based half-open.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

from privy.synteny.coordinates import PathCoordinateModel
from privy.synteny.graph_blocks import build_pairwise_blocks
from privy.synteny.model import GenomeInterval, SyntenyBlock, SyntenyRegion, split_pansn


@dataclass(frozen=True)
class SyntenyResult:
    """All blocks (every query vs the reference) plus grouped regions."""

    reference: str
    blocks: tuple[SyntenyBlock, ...]
    regions: tuple[SyntenyRegion, ...]


@dataclass(frozen=True)
class RegionPrivacy:
    """Cohort presence + private-region verdict for one synteny region."""

    region_id: str
    present: tuple[str, ...]          # genomes traversing the region's ref segments
    target_present: tuple[str, ...]
    offtarget_present: tuple[str, ...]
    target_private: bool              # present in ≥1 target and no off-target


def build_synteny(
    model: PathCoordinateModel,
    ref_path: str,
    query_paths: Sequence[str] | None = None,
) -> SyntenyResult:
    """Build typed synteny blocks for each query vs *ref_path*, grouped into regions.

    Args:
        query_paths: Genomes to compare against the reference (default: all paths
            except the reference).
    """
    queries = list(query_paths) if query_paths is not None else [
        p for p in model.path_ids() if p != ref_path
    ]
    blocks: list[SyntenyBlock] = []
    for query in queries:
        blocks.extend(build_pairwise_blocks(model, query, ref_path))
    regions = group_regions(blocks, ref_genome=split_pansn(ref_path)[0])
    return SyntenyResult(reference=ref_path, blocks=tuple(blocks), regions=tuple(regions))


def group_regions(
    blocks: Sequence[SyntenyBlock],
    *,
    ref_genome: str | None = None,
) -> list[SyntenyRegion]:
    """Merge blocks whose reference (target) intervals overlap into regions.

    Blocks are grouped per reference contig; overlapping reference spans become one
    region carrying all contributing blocks.  The region's reference genome label is
    *ref_genome* when given, else taken from the blocks on that contig (so PAF-mode
    blocks, which have no single reference path, still group correctly).
    """
    by_contig: dict[str, list[SyntenyBlock]] = {}
    for block in blocks:
        by_contig.setdefault(block.target.contig, []).append(block)

    regions: list[SyntenyRegion] = []
    region_idx = 0
    for contig in sorted(by_contig):
        ordered = sorted(by_contig[contig], key=lambda b: (b.target.start, b.target.end))
        contig_genome: str = (
            ref_genome if ref_genome is not None else ordered[0].target.genome
        )
        cur: list[SyntenyBlock] = []
        cur_start = cur_end = 0
        for block in ordered:
            if not cur:
                cur = [block]
                cur_start, cur_end = block.target.start, block.target.end
                continue
            if block.target.start < cur_end:   # overlaps current region
                cur.append(block)
                cur_end = max(cur_end, block.target.end)
            else:
                regions.append(
                    _make_region(contig_genome, contig, cur_start, cur_end, cur, region_idx)
                )
                region_idx += 1
                cur = [block]
                cur_start, cur_end = block.target.start, block.target.end
        if cur:
            regions.append(
                _make_region(contig_genome, contig, cur_start, cur_end, cur, region_idx)
            )
            region_idx += 1
    return regions


def _make_region(
    ref_genome: str,
    contig: str,
    start: int,
    end: int,
    blocks: list[SyntenyBlock],
    idx: int,
) -> SyntenyRegion:
    return SyntenyRegion(
        region_id=f"SR{idx:06d}",
        reference=GenomeInterval(ref_genome, contig, start, end),
        blocks=tuple(blocks),
    )


def tag_region_privacy(
    model: PathCoordinateModel,
    regions: Sequence[SyntenyRegion],
    targets: Sequence[str],
    off_targets: Sequence[str],
    *,
    min_present_fraction: float = 1.0,
) -> dict[str, RegionPrivacy]:
    """Flag regions whose structure is private to the target cohort.

    A genome is "present" in a region when it traverses at least
    *min_present_fraction* of the region's reference segments.  A region is
    ``target_private`` when ≥1 target genome is present and no off-target is.

    *targets* / *off_targets* are path ids (e.g. from PanSN cohort resolution).
    """
    target_set = list(dict.fromkeys(targets))
    offtarget_set = list(dict.fromkeys(off_targets))
    result: dict[str, RegionPrivacy] = {}

    for region in regions:
        segments = {a.name for b in region.blocks for a in b.anchors if a.name}
        n_segs = len(segments)

        def _present(path: str, segs: set[str] = segments, total: int = n_segs) -> bool:
            if total == 0 or path not in model:
                return False
            hit = sum(1 for s in segs if model.occurrences(path, s))
            return hit / total >= min_present_fraction

        target_present = tuple(p for p in target_set if _present(p))
        offtarget_present = tuple(p for p in offtarget_set if _present(p))
        present = target_present + offtarget_present
        result[region.region_id] = RegionPrivacy(
            region_id=region.region_id,
            present=present,
            target_present=target_present,
            offtarget_present=offtarget_present,
            target_private=bool(target_present) and not offtarget_present,
        )
    return result
