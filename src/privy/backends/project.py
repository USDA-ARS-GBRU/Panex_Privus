"""Backend for ``privy project``: project a region or node-set to references.

Loads a GFA pangenome, builds a :class:`~privy.synteny.coordinates.PathCoordinateModel`,
projects a source region (given as ``path:start-end`` in the source path's stable
coordinates) or a raw set of graph segments onto every requested target genome,
and writes ``projection.tsv`` + ``project.json``.
"""

from __future__ import annotations

import csv
import json
import logging
from collections.abc import Sequence
from pathlib import Path

from privy.io.gfa import parse_gfa
from privy.synteny.coordinates import PathCoordinateModel
from privy.synteny.model import ProjectionMap
from privy.synteny.projection import project_node_set, project_region

log = logging.getLogger("privy.backends.project")

PROJECTION_COLUMNS = ["target_path", "genome", "contig", "start", "end", "present"]


def run_project(
    gfa_path: Path,
    *,
    source_path: str | None = None,
    region: tuple[int, int] | None = None,
    node_set: Sequence[str] | None = None,
    to_genomes: Sequence[str] | None = None,
    outdir: Path,
) -> list[Path]:
    """Project a region/node-set across a pangenome graph and write outputs.

    Provide either *node_set* (segments defining the region in graph node space)
    or both *source_path* and *region* (a ``(start, end)`` interval in that path's
    stable coordinates).

    Returns:
        Paths written: ``projection.tsv`` and ``project.json``.

    Raises:
        ValueError: On missing/contradictory inputs or an unknown source path.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    graph = parse_gfa(Path(gfa_path))
    model = PathCoordinateModel.from_graph(graph)
    targets = list(to_genomes) if to_genomes else None

    if node_set:
        pm = project_node_set(
            model, list(node_set), targets=targets,
            source_label="node-set:" + ",".join(node_set),
        )
    elif region is not None and source_path is not None:
        if source_path not in model:
            raise ValueError(f"source path {source_path!r} not found in graph")
        start, end = region
        if end <= start:
            raise ValueError(f"region end ({end}) must be > start ({start})")
        local_start = model.to_path_local(source_path, start)
        local_end = model.to_path_local(source_path, end - 1) + 1
        pm = project_region(model, source_path, local_start, local_end, targets=targets)
    else:
        raise ValueError("provide either --node-set or both --source-path and --region")

    proj_path = _write_projection_tsv(pm, outdir)
    meta_path = _write_metadata(pm, Path(gfa_path), outdir)
    log.info(
        "Projection complete | source=%s present=%d absent=%d",
        pm.source, len(pm.present_in()), len(pm.absent_in()),
    )
    return [proj_path, meta_path]


def _write_projection_tsv(pm: ProjectionMap, outdir: Path) -> Path:
    path = outdir / "projection.tsv"
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(PROJECTION_COLUMNS)
        for target_path, interval in pm.projections.items():
            if interval is None:
                writer.writerow([target_path, "", "", "", "", "False"])
            else:
                writer.writerow([
                    target_path,
                    interval.genome,
                    interval.contig,
                    interval.start,
                    interval.end,
                    "True",
                ])
    return path


def _write_metadata(pm: ProjectionMap, gfa_path: Path, outdir: Path) -> Path:
    path = outdir / "project.json"
    metadata = {
        "tool": "privy project",
        "gfa": str(gfa_path),
        "source": pm.source,
        "n_targets": len(pm.projections),
        "n_present": len(pm.present_in()),
        "n_absent": len(pm.absent_in()),
        "present_in": list(pm.present_in()),
        "absent_in": list(pm.absent_in()),
    }
    path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    return path
