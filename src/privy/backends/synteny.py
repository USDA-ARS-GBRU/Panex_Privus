"""Backend for ``privy synteny``: graph-native synteny blocks + private regions.

Loads a GFA pangenome, builds typed synteny blocks for each query genome against a
chosen reference, groups them into regions, optionally flags target-private
structural regions, and writes ``synteny_blocks.tsv``, ``synteny_regions.tsv``, and
``synteny.json``.
"""

from __future__ import annotations

import csv
import json
import logging
from collections.abc import Sequence
from pathlib import Path

from privy.io.gfa import parse_gfa
from privy.synteny.build import (
    RegionPrivacy,
    SyntenyResult,
    build_synteny,
    tag_region_privacy,
)
from privy.synteny.coordinates import PathCoordinateModel
from privy.synteny.model import SyntenyBlock, SyntenyRegion, split_pansn

log = logging.getLogger("privy.backends.synteny")

BLOCK_COLUMNS = [
    "block_id", "query_genome", "query_contig", "query_start", "query_end",
    "ref_genome", "ref_contig", "ref_start", "ref_end",
    "strand", "block_type", "n_anchors", "score",
]
REGION_COLUMNS = [
    "region_id", "ref_contig", "ref_start", "ref_end", "n_blocks", "genomes",
    "target_private", "target_present", "offtarget_present",
]


def run_synteny(
    gfa_path: Path,
    *,
    reference: str,
    query_paths: Sequence[str] | None = None,
    targets: Sequence[str] | None = None,
    off_targets: Sequence[str] | None = None,
    min_block_anchors: int = 1,
    outdir: Path,
) -> list[Path]:
    """Build synteny for a graph and write block/region tables.

    Args:
        reference: Path id used as the coordinate reference.
        query_paths: Paths to compare (default: all except the reference).
        targets / off_targets: Cohort sample names or path ids; when given,
            regions are tagged target-private.
        min_block_anchors: Minimum anchors for a reported collinear/inversion block.

    Raises:
        ValueError: If the reference path is unknown.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    model = PathCoordinateModel.from_graph(parse_gfa(Path(gfa_path)))
    if reference not in model:
        raise ValueError(f"reference path {reference!r} not found in graph")

    result = build_synteny(model, reference, query_paths)

    privacy: dict[str, RegionPrivacy] = {}
    if targets or off_targets:
        target_paths = _resolve_paths(model, targets)
        offtarget_paths = _resolve_paths(model, off_targets)
        privacy = tag_region_privacy(model, result.regions, target_paths, offtarget_paths)

    if min_block_anchors > 1:
        kept = {b.block_id for b in result.blocks if b.n_anchors >= min_block_anchors}
    else:
        kept = {b.block_id for b in result.blocks}

    blocks_path = _write_blocks(result.blocks, kept, outdir)
    regions_path = _write_regions(result.regions, privacy, outdir)
    meta_path = _write_metadata(result, privacy, Path(gfa_path), outdir)

    n_private = sum(1 for v in privacy.values() if v.target_private)
    log.info(
        "Synteny complete | blocks=%d regions=%d private=%d",
        len(kept), len(result.regions), n_private,
    )
    return [blocks_path, regions_path, meta_path]


def _resolve_paths(model: PathCoordinateModel, names: Sequence[str] | None) -> list[str]:
    """Resolve cohort tokens (path ids or PanSN sample names) to path ids."""
    if not names:
        return []
    path_ids = model.path_ids()
    path_set = set(path_ids)
    out: list[str] = []
    for name in names:
        if name in path_set:
            out.append(name)
        else:
            out.extend(p for p in path_ids if split_pansn(p)[0] == name)
    return list(dict.fromkeys(out))


def _write_blocks(blocks: Sequence[SyntenyBlock], kept: set[str], outdir: Path) -> Path:
    path = outdir / "synteny_blocks.tsv"
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(BLOCK_COLUMNS)
        for b in blocks:
            if b.block_id not in kept:
                continue
            writer.writerow([
                b.block_id,
                b.query.genome, b.query.contig, b.query.start, b.query.end,
                b.target.genome, b.target.contig, b.target.start, b.target.end,
                b.strand, b.block_type.value, b.n_anchors,
                "" if b.score is None else b.score,
            ])
    return path


def _write_regions(
    regions: Sequence[SyntenyRegion],
    privacy: dict[str, RegionPrivacy],
    outdir: Path,
) -> Path:
    path = outdir / "synteny_regions.tsv"
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(REGION_COLUMNS)
        for region in regions:
            verdict = privacy.get(region.region_id)
            if verdict is None:
                tp, t_present, o_present = "NA", "", ""
            else:
                tp = str(verdict.target_private)
                t_present = ",".join(verdict.target_present)
                o_present = ",".join(verdict.offtarget_present)
            writer.writerow([
                region.region_id,
                region.reference.contig, region.reference.start, region.reference.end,
                region.n_blocks, ",".join(region.genomes),
                tp, t_present, o_present,
            ])
    return path


def _write_metadata(
    result: SyntenyResult,
    privacy: dict[str, RegionPrivacy],
    gfa_path: Path,
    outdir: Path,
) -> Path:
    path = outdir / "synteny.json"
    metadata = {
        "tool": "privy synteny",
        "gfa": str(gfa_path),
        "reference": result.reference,
        "n_blocks": len(result.blocks),
        "n_regions": len(result.regions),
        "n_target_private_regions": sum(1 for v in privacy.values() if v.target_private),
        "block_type_counts": _block_type_counts(result.blocks),
    }
    path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    return path


def _block_type_counts(blocks: Sequence[SyntenyBlock]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for b in blocks:
        counts[b.block_type.value] = counts.get(b.block_type.value, 0) + 1
    return counts
