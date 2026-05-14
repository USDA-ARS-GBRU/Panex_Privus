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
from privy.io.tsv import TsvWriter, read_tsv
from privy.pangenome import (
    build_composition_rows,
    build_coverage_histogram_rows,
    build_feature_summary_rows,
    build_gfa_feature_matrix,
    build_growth_curve_rows,
    build_vcf_feature_matrix,
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
    write_plots: bool = False,
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


def run_pangenome_vcf(
    vcf: Path,
    targets: list[str],
    outdir: Path,
    off_targets: list[str] | None = None,
    ignored_samples: list[str] | None = None,
    permutations: int = 100,
    seed: int = 42,
    write_plots: bool = False,
    plot_format: str = "png",
) -> None:
    """Run VCF allele-level pangenome analysis."""
    if not vcf.exists():
        raise FileNotFoundError(f"VCF file not found: {vcf}")

    outdir.mkdir(parents=True, exist_ok=True)
    log.info("Parsing VCF for pangenome analysis: %s", vcf)
    matrix = build_vcf_feature_matrix(vcf)
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
        input_path=vcf,
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
    write_plots: bool = False,
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

        log.info("Rendering pangenome plots | format=%s", plot_format)
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
        log.info("Rendered pangenome plots | count=%d", len(plot_paths))

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


def run_pangenome_plots(
    input_dir: Path,
    plot_format: str = "png",
    outdir: Path | None = None,
) -> list[Path]:
    """Render pangenome plots from existing pangenome TSV outputs."""
    plot_outdir = outdir or input_dir
    coverage_path = input_dir / "coverage_histogram.tsv"
    composition_path = input_dir / "composition.tsv"
    growth_path = input_dir / "growth_curves.tsv"
    for path in (coverage_path, composition_path, growth_path):
        if not path.exists():
            raise FileNotFoundError(f"Pangenome table not found: {path}")

    log.info(
        "Loading existing pangenome tables for plotting | input_dir=%s | outdir=%s",
        input_dir,
        plot_outdir,
    )
    coverage_rows = read_tsv(coverage_path)
    composition_rows = read_tsv(composition_path)
    growth_rows = read_tsv(growth_path)

    from privy.plot.pangenome import plot_all_pangenome  # noqa: PLC0415

    log.info(
        "Rendering pangenome plots from existing tables | format=%s | "
        "coverage_rows=%d | composition_rows=%d | growth_rows=%d",
        plot_format,
        len(coverage_rows),
        len(composition_rows),
        len(growth_rows),
    )
    paths = plot_all_pangenome(
        coverage_rows=coverage_rows,
        composition_rows=composition_rows,
        growth_rows=growth_rows,
        outdir=plot_outdir,
        output_format=plot_format,
    )
    if plot_outdir == input_dir:
        _update_plot_metadata(input_dir, plot_format, paths)
    log.info("Rendered pangenome plots from existing tables | count=%d", len(paths))
    return paths


def _update_plot_metadata(
    outdir: Path,
    plot_format: str,
    plot_paths: list[Path],
) -> None:
    metadata_path = outdir / "pangenome.json"
    if not metadata_path.exists():
        return
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    parameters = metadata.setdefault("parameters", {})
    if isinstance(parameters, dict):
        parameters["write_plots"] = True
        parameters["plot_format"] = plot_format
    outputs = metadata.setdefault("outputs", [])
    if isinstance(outputs, list):
        existing = {str(item) for item in outputs}
        for path in plot_paths:
            if path.name not in existing:
                outputs.append(path.name)
                existing.add(path.name)
    metadata_path.write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
