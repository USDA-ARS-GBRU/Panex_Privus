"""Backend for ``privy popgen``: breeder population-genetics summaries.

Detects microhaplotype loci from a pangenome graph, then writes per-locus allelic
diversity and target-vs-off-target differentiation, surfacing fully diagnostic
(target-private) markers — the breeder-actionable output.  Pure-Python.
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
from privy.polyploid.dosage import group_paths_by_sample
from privy.popgen.differentiation import genome_wide_fst, locus_differentiation
from privy.popgen.diversity import (
    allele_frequencies,
    locus_diversity,
    nei_gene_diversity,
)
from privy.popgen.relationship import DosageMatrix, build_dosage_matrix, vanraden_grm
from privy.synteny.coordinates import PathCoordinateModel
from privy.synteny.model import split_pansn

log = logging.getLogger("privy.backends.popgen")

POPGEN_COLUMNS = [
    "locus_id", "contig", "start", "end", "n_genomes", "n_alleles",
    "gene_diversity", "effective_alleles", "aaf",
    "target_he", "offtarget_he", "gst", "jost_d", "is_diagnostic",
]


def run_popgen(
    gfa_path: Path,
    *,
    reference: str,
    targets: Sequence[str],
    off_targets: Sequence[str],
    min_core_fraction: float = 1.0,
    outdir: Path,
) -> list[Path]:
    """Compute population-genetics summaries for a graph and write tables.

    Raises:
        ValueError: If the reference path is unknown or no cohorts are given.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    graph = parse_gfa(Path(gfa_path))
    model = PathCoordinateModel.from_graph(graph)
    if reference not in model:
        raise ValueError(f"reference path {reference!r} not found in graph")

    target_paths = _resolve_paths(model, targets)
    offtarget_paths = _resolve_paths(model, off_targets)
    if not target_paths or not offtarget_paths:
        raise ValueError("both --targets and --off-targets must resolve to genomes in the graph")

    loci = detect_microhaplotypes(graph, model, reference, min_core_fraction=min_core_fraction)

    loci_path = _write_loci(loci, target_paths, offtarget_paths, outdir)
    fst = genome_wide_fst(loci, target_paths, offtarget_paths)
    n_diagnostic = sum(
        1 for mh in loci
        if (d := locus_differentiation(mh, target_paths, offtarget_paths)) and d.is_diagnostic
    )

    # GP-ready exports: per-sample dosage matrix + VanRaden genomic relationship matrix.
    samples_to_paths = group_paths_by_sample(model.path_ids())
    written = [loci_path]
    if loci and samples_to_paths:
        dm = build_dosage_matrix(loci, samples_to_paths)
        written.append(_write_dosage_matrix(dm, outdir))
        written.append(_write_grm(dm, outdir))

    meta_path = _write_metadata(
        loci, fst, n_diagnostic, target_paths, offtarget_paths, Path(gfa_path), reference, outdir
    )
    written.append(meta_path)

    log.info(
        "Popgen complete | loci=%d genome_wide_fst=%.3f diagnostic=%d",
        len(loci), fst, n_diagnostic,
    )
    return written


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


def _write_loci(
    loci: Sequence[Microhaplotype],
    targets: Sequence[str],
    off_targets: Sequence[str],
    outdir: Path,
) -> Path:
    path = outdir / "popgen_loci.tsv"
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(POPGEN_COLUMNS)
        for mh in loci:
            div = locus_diversity(mh)
            diff = locus_differentiation(mh, targets, off_targets)
            if diff is None:
                t_he = o_he = gst = jd = ""
                diag = "NA"
            else:
                t_he = f"{nei_gene_diversity(allele_frequencies(mh.alleles, targets)):.4f}"
                o_he = f"{nei_gene_diversity(allele_frequencies(mh.alleles, off_targets)):.4f}"
                gst = f"{diff.gst:.4f}"
                jd = f"{diff.jost_d:.4f}"
                diag = str(diff.is_diagnostic)
            writer.writerow([
                div.locus_id, mh.contig, mh.start, mh.end, div.n_genomes, div.n_alleles,
                f"{div.gene_diversity:.4f}", f"{div.effective_alleles:.4f}", f"{div.aaf:.4f}",
                t_he, o_he, gst, jd, diag,
            ])
    return path


def _write_dosage_matrix(dm: DosageMatrix, outdir: Path) -> Path:
    """Write the samples × loci alt-allele dosage matrix (GP input)."""
    path = outdir / "dosage_matrix.tsv"
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", *dm.locus_ids])
        for sample, row in zip(dm.samples, dm.matrix, strict=True):
            writer.writerow([sample, *["." if v is None else v for v in row]])
    return path


def _write_grm(dm: DosageMatrix, outdir: Path) -> Path:
    """Write the VanRaden genomic relationship matrix (labelled square matrix)."""
    path = outdir / "grm.tsv"
    samples, grm = vanraden_grm(dm)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", *samples])
        for sample, row in zip(samples, grm, strict=True):
            writer.writerow([sample, *row])
    return path


def _write_metadata(
    loci: Sequence[Microhaplotype],
    fst: float,
    n_diagnostic: int,
    targets: Sequence[str],
    off_targets: Sequence[str],
    gfa_path: Path,
    reference: str,
    outdir: Path,
) -> Path:
    path = outdir / "popgen.json"
    metadata = {
        "tool": "privy popgen",
        "gfa": str(gfa_path),
        "reference": reference,
        "n_loci": len(loci),
        "genome_wide_fst": round(fst, 4),
        "n_diagnostic_loci": n_diagnostic,
        "n_targets": len(targets),
        "n_off_targets": len(off_targets),
    }
    path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    return path
