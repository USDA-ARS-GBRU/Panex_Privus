"""Backend for ``privy microhap``: detect microhaplotypes and write allele tables.

Loads a GFA pangenome, detects multi-allelic microhaplotype loci between core
segments, optionally flags target-private alleles, and writes
``microhaplotypes.tsv`` (one row per locus), ``allele_matrix.tsv`` (loci × genomes,
integer allele indices), and ``microhap.json``.
"""

from __future__ import annotations

import csv
import json
import logging
from collections.abc import Sequence
from pathlib import Path

from privy.io.gfa import parse_gfa
from privy.microhap.detect import detect_microhaplotypes
from privy.microhap.model import Microhaplotype
from privy.synteny.coordinates import PathCoordinateModel
from privy.synteny.model import split_pansn

log = logging.getLogger("privy.backends.microhap")

MICROHAP_COLUMNS = [
    "locus_id", "contig", "start", "end", "n_genomes", "n_alleles",
    "aaf", "ref_allele", "target_private", "n_private_alleles",
]


def run_microhap(
    gfa_path: Path,
    *,
    reference: str,
    targets: Sequence[str] | None = None,
    off_targets: Sequence[str] | None = None,
    min_core_fraction: float = 1.0,
    multiallelic_only: bool = True,
    outdir: Path,
) -> list[Path]:
    """Detect microhaplotypes for a graph and write loci + allele-matrix tables.

    Raises:
        KeyError / ValueError: If the reference path is unknown.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    graph = parse_gfa(Path(gfa_path))
    model = PathCoordinateModel.from_graph(graph)
    if reference not in model:
        raise ValueError(f"reference path {reference!r} not found in graph")

    loci = detect_microhaplotypes(
        graph, model, reference,
        min_core_fraction=min_core_fraction,
        multiallelic_only=multiallelic_only,
    )

    target_paths = _resolve_paths(model, targets)
    offtarget_paths = _resolve_paths(model, off_targets)
    have_cohorts = bool(target_paths or off_targets)

    loci_path = _write_loci(loci, target_paths, offtarget_paths, have_cohorts, outdir)
    matrix_path = _write_allele_matrix(loci, model.path_ids(), outdir)
    meta_path = _write_metadata(
        loci, target_paths, offtarget_paths, have_cohorts, Path(gfa_path), reference, outdir
    )

    n_private = (
        sum(1 for m in loci if m.is_target_private(target_paths, offtarget_paths))
        if have_cohorts else 0
    )
    log.info(
        "Microhap complete | loci=%d target_private=%d", len(loci), n_private
    )
    return [loci_path, matrix_path, meta_path]


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


def _allele_index_map(mh: Microhaplotype, path_order: Sequence[str]) -> dict[str, int]:
    """Assign integer indices to alleles: reference allele 0, then first-seen order."""
    order: list[str] = []
    if mh.ref_allele is not None:
        order.append(mh.ref_allele)
    for path in path_order:
        allele = mh.alleles.get(path)
        if allele is not None and allele not in order:
            order.append(allele)
    return {allele: idx for idx, allele in enumerate(order)}


def _write_loci(
    loci: Sequence[Microhaplotype],
    targets: Sequence[str],
    off_targets: Sequence[str],
    have_cohorts: bool,
    outdir: Path,
) -> Path:
    path = outdir / "microhaplotypes.tsv"
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(MICROHAP_COLUMNS)
        for mh in loci:
            if have_cohorts:
                private = mh.private_alleles(targets, off_targets)
                tp, n_priv = str(bool(private)), len(private)
            else:
                tp, n_priv = "NA", 0
            writer.writerow([
                mh.locus_id, mh.contig, mh.start, mh.end,
                mh.n_genomes, mh.n_alleles, f"{mh.aaf():.4f}",
                "" if mh.ref_allele is None else mh.ref_allele[:12],
                tp, n_priv,
            ])
    return path


def _write_allele_matrix(
    loci: Sequence[Microhaplotype],
    path_order: Sequence[str],
    outdir: Path,
) -> Path:
    path = outdir / "allele_matrix.tsv"
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["locus_id", *path_order])
        for mh in loci:
            idx_map = _allele_index_map(mh, path_order)
            row = [mh.locus_id]
            for genome in path_order:
                allele = mh.alleles.get(genome)
                row.append("." if allele is None else str(idx_map[allele]))
            writer.writerow(row)
    return path


def _write_metadata(
    loci: Sequence[Microhaplotype],
    targets: Sequence[str],
    off_targets: Sequence[str],
    have_cohorts: bool,
    gfa_path: Path,
    reference: str,
    outdir: Path,
) -> Path:
    path = outdir / "microhap.json"
    n_private = (
        sum(1 for m in loci if m.is_target_private(targets, off_targets))
        if have_cohorts else 0
    )
    metadata = {
        "tool": "privy microhap",
        "gfa": str(gfa_path),
        "reference": reference,
        "n_loci": len(loci),
        "n_target_private_loci": n_private,
        "max_alleles": max((m.n_alleles for m in loci), default=0),
    }
    path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    return path
