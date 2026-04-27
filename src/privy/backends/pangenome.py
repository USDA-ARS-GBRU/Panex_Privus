"""Pangenome analysis backend.

This backend turns source-specific inputs into the shared pangenome feature
matrix, then writes source-independent analysis tables.  GFA is implemented
first; VCF can join by producing the same matrix shape.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any

from privy.io.gfa import parse_gfa
from privy.io.tsv import TsvWriter
from privy.pangenome import (
    build_composition_rows,
    build_coverage_histogram_rows,
    build_feature_summary_rows,
    build_gfa_feature_matrix,
    build_growth_curve_rows,
    resolve_pangenome_groups,
)
from privy.pangenome.model import FeatureMatrix, PangenomeGroups
from privy.utils.misc import now_iso

log = logging.getLogger("privy.backends.pangenome")

FEATURE_SUMMARY_COLUMNS = [
    "feature_id",
    "source_type",
    "feature_type",
    "contig",
    "start",
    "end",
    "length",
    "total_present_n",
    "target_present_n",
    "target_total_n",
    "offtarget_present_n",
    "offtarget_total_n",
    "full_category",
    "target_category",
    "offtarget_category",
    "target_private",
    "offtarget_private",
]

COVERAGE_HISTOGRAM_COLUMNS = ["group", "coverage", "n_features", "n_bp"]
COMPOSITION_COLUMNS = ["group", "category", "n_features", "n_bp"]
GROWTH_CURVE_COLUMNS = [
    "group",
    "trial",
    "n",
    "sample_added",
    "features",
    "bp",
    "new_features",
    "new_bp",
    "singleton_features",
    "singleton_bp",
]


def run_pangenome_gfa(
    gfa: Path,
    targets: list[str],
    outdir: Path,
    off_targets: list[str] | None = None,
    ignored_samples: list[str] | None = None,
    permutations: int = 100,
    seed: int = 42,
    write_plots: bool = True,
    plot_format: str = "png",
) -> None:
    """Run the first-pass GFA pangenome analysis."""
    if not gfa.exists():
        raise FileNotFoundError(f"GFA file not found: {gfa}")

    outdir.mkdir(parents=True, exist_ok=True)
    log.info("Parsing GFA for pangenome analysis: %s", gfa)
    graph = parse_gfa(gfa)
    matrix = build_gfa_feature_matrix(graph)
    groups = resolve_pangenome_groups(
        all_samples=matrix.samples,
        targets=targets,
        off_targets=off_targets,
        ignored_samples=ignored_samples,
    )
    write_pangenome_outputs(
        matrix=matrix,
        groups=groups,
        outdir=outdir,
        input_path=gfa,
        permutations=permutations,
        seed=seed,
        write_plots=write_plots,
        plot_format=plot_format,
    )


def write_pangenome_outputs(
    matrix: FeatureMatrix,
    groups: PangenomeGroups,
    outdir: Path,
    input_path: Path,
    permutations: int,
    seed: int,
    write_plots: bool = True,
    plot_format: str = "png",
) -> None:
    """Write shared pangenome analysis outputs."""
    feature_rows = build_feature_summary_rows(matrix, groups)
    coverage_rows = build_coverage_histogram_rows(matrix, groups)
    composition_rows = build_composition_rows(matrix, groups)
    growth_rows = build_growth_curve_rows(matrix, groups, permutations=permutations, seed=seed)

    with TsvWriter(outdir / "feature_summary.tsv", FEATURE_SUMMARY_COLUMNS) as writer:
        writer.write_rows(feature_rows)
    with TsvWriter(outdir / "coverage_histogram.tsv", COVERAGE_HISTOGRAM_COLUMNS) as writer:
        writer.write_rows(coverage_rows)
    with TsvWriter(outdir / "composition.tsv", COMPOSITION_COLUMNS) as writer:
        writer.write_rows(composition_rows)
    with TsvWriter(outdir / "growth_curves.tsv", GROWTH_CURVE_COLUMNS) as writer:
        writer.write_rows(growth_rows)

    plot_paths: list[str] = []
    if write_plots:
        from privy.plot.pangenome import plot_all_pangenome  # noqa: PLC0415

        plot_paths = [
            str(path.name)
            for path in plot_all_pangenome(
                coverage_rows=coverage_rows,
                composition_rows=composition_rows,
                growth_rows=growth_rows,
                outdir=outdir,
                output_format=plot_format,
            )
        ]

    metadata: dict[str, Any] = {
        "created_at": now_iso(),
        "analysis": "pangenome",
        "source_type": matrix.source_type,
        "inputs": {
            matrix.source_type: str(input_path),
        },
        "parameters": {
            "permutations": permutations,
            "plot_format": plot_format,
            "seed": seed,
            "write_plots": write_plots,
        },
        "samples": {
            "full": list(groups.full),
            "target": list(groups.target),
            "off_target": list(groups.off_target),
            "ignored": list(groups.ignored),
        },
        "summary": {
            "n_features": len(matrix.features),
            "n_samples": len(matrix.samples),
            "n_target_samples": len(groups.target),
            "n_offtarget_samples": len(groups.off_target),
        },
        "outputs": [
            "feature_summary.tsv",
            "coverage_histogram.tsv",
            "composition.tsv",
            "growth_curves.tsv",
            *plot_paths,
            "pangenome.json",
        ],
    }
    (outdir / "pangenome.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
