"""Windowed VCF landscape backend."""

from __future__ import annotations

import json
import logging
from contextlib import ExitStack
from pathlib import Path
from typing import Any

from privy.io.tsv import TsvWriter, read_tsv
from privy.landscape import (
    BACKGROUND_BLOCK_COLUMNS,
    CANDIDATE_INTROGRESSION_BLOCK_COLUMNS,
    LANDSCAPE_LOCAL_PCA_COLUMNS,
    LANDSCAPE_SAMPLE_WINDOW_COLUMNS,
    LANDSCAPE_SIMILARITY_COLUMNS,
    LANDSCAPE_WINDOW_COLUMNS,
    run_vcf_landscape,
)
from privy.utils.misc import now_iso

log = logging.getLogger("privy.backends.landscape")

SIMILARITY_OUTPUT_MODES = {"full", "summary", "none"}
VCF_ENGINES = {"auto", "pysam", "cyvcf2"}
LANDSCAPE_PLOT_SCOPES = {"chromosome", "genome", "both"}
LANDSCAPE_PLOT_INDEX_COLUMNS = [
    "plot_type",
    "plot_scope",
    "contig",
    "path",
    "n_windows",
    "start",
    "end",
    "output_format",
]


def run_landscape_vcf(
    vcf: Path,
    targets: list[str],
    outdir: Path,
    off_targets: list[str] | None = None,
    ignored_samples: list[str] | None = None,
    window_records: int = 200,
    step_records: int = 50,
    window_bp: int | None = None,
    step_bp: int | None = None,
    pass_only: bool = True,
    min_qual: float | None = None,
    rare_max_count: int = 1,
    rare_max_freq: float = 0.05,
    min_called_for_freq: int = 10,
    min_freq_values: int = 10,
    min_background_similarity: float = 0.65,
    min_introgression_similarity: float | None = None,
    min_introgression_delta: float = 0.05,
    max_introgression_missing_rate: float = 0.5,
    min_introgression_windows: int = 10,
    similarity_output: str = "full",
    vcf_engine: str = "auto",
    local_pca: bool = False,
    write_plots: bool = False,
    plot_format: str = "png",
) -> None:
    """Run a VCF landscape analysis and write tables, plots, and metadata."""
    if similarity_output not in SIMILARITY_OUTPUT_MODES:
        raise ValueError(
            "--similarity-output must be one of: full, summary, none."
        )
    if vcf_engine not in VCF_ENGINES:
        raise ValueError("--vcf-engine must be one of: auto, pysam, cyvcf2.")

    outdir.mkdir(parents=True, exist_ok=True)
    log.info("Running VCF landscape analysis: %s", vcf)
    log.info(
        "Streaming large landscape tables while analyzing | outdir=%s | "
        "similarity_output=%s | vcf_engine=%s | local_pca=%s",
        outdir,
        similarity_output,
        vcf_engine,
        local_pca,
    )
    with ExitStack() as stack:
        sample_writer = stack.enter_context(
            TsvWriter(outdir / "sample_windows.tsv", LANDSCAPE_SAMPLE_WINDOW_COLUMNS)
        )
        window_writer = stack.enter_context(
            TsvWriter(outdir / "windows.tsv", LANDSCAPE_WINDOW_COLUMNS)
        )
        similarity_writer = (
            stack.enter_context(
                TsvWriter(outdir / "similarity.tsv", LANDSCAPE_SIMILARITY_COLUMNS)
            )
            if similarity_output == "full"
            else None
        )
        local_pca_writer = (
            stack.enter_context(
                TsvWriter(outdir / "local_pca.tsv", LANDSCAPE_LOCAL_PCA_COLUMNS)
            )
            if local_pca
            else None
        )
        result = run_vcf_landscape(
            vcf_path=vcf,
            targets=targets,
            off_targets=off_targets,
            ignored_samples=ignored_samples,
            window_records=window_records,
            step_records=step_records,
            window_bp=window_bp,
            step_bp=step_bp,
            pass_only=pass_only,
            min_qual=min_qual,
            rare_max_count=rare_max_count,
            rare_max_freq=rare_max_freq,
            min_called_for_freq=min_called_for_freq,
            min_freq_values=min_freq_values,
            min_background_similarity=min_background_similarity,
            min_introgression_similarity=min_introgression_similarity,
            min_introgression_delta=min_introgression_delta,
            max_introgression_missing_rate=max_introgression_missing_rate,
            min_introgression_windows=min_introgression_windows,
            sample_row_writer=sample_writer,
            window_row_writer=window_writer,
            similarity_row_writer=similarity_writer,
            local_pca_row_writer=local_pca_writer,
            emit_similarity_rows=similarity_output == "full",
            retain_output_rows=write_plots,
            retain_similarity_rows=False,
            retain_local_pca_rows=False,
            local_pca=local_pca,
            vcf_engine=vcf_engine,
        )

    log.info(
        "Landscape analysis complete | windows=%d | sample_window_rows=%d | "
        "background_blocks=%d | candidate_introgression_blocks=%d | "
        "similarity_rows=%d | local_pca_rows=%d",
        result.n_window_rows,
        result.n_sample_rows,
        len(result.background_block_rows),
        len(result.candidate_introgression_rows),
        result.n_similarity_rows,
        result.n_local_pca_rows,
    )

    log.info("Writing landscape block tables | outdir=%s", outdir)
    with TsvWriter(outdir / "background_blocks.tsv", BACKGROUND_BLOCK_COLUMNS) as writer:
        writer.write_rows(result.background_block_rows)
    with TsvWriter(
        outdir / "candidate_introgression_blocks.tsv",
        CANDIDATE_INTROGRESSION_BLOCK_COLUMNS,
    ) as writer:
        writer.write_rows(result.candidate_introgression_rows)
    if similarity_output == "summary":
        with TsvWriter(outdir / "similarity.tsv", LANDSCAPE_SIMILARITY_COLUMNS) as writer:
            writer.write_rows(result.similarity_summary_rows)

    plot_paths: list[str] = []
    if write_plots:
        from privy.plot.landscape import plot_all_landscape  # noqa: PLC0415

        log.info("Rendering landscape plots | format=%s", plot_format)
        plot_paths = [
            str(path.name)
            for path in plot_all_landscape(
                sample_rows=result.sample_rows,
                window_rows=result.window_rows,
                similarity_rows=(
                    result.similarity_rows or result.similarity_summary_rows
                ),
                outdir=outdir,
                output_format=plot_format,
            )
        ]
        log.info("Rendered landscape plots | count=%d", len(plot_paths))
    else:
        log.info("Skipping landscape plots (--no-plots)")

    window_parameters: dict[str, Any] = {
        "window_mode": result.window_mode,
        "window_records": window_records,
        "step_records": step_records,
        "window_bp": window_bp,
        "step_bp": step_bp,
    }
    metadata: dict[str, Any] = {
        "created_at": now_iso(),
        "analysis": "landscape",
        "source_type": "vcf",
        "inputs": {"vcf": str(vcf)},
        "parameters": {
            **window_parameters,
            "pass_only": pass_only,
            "min_qual": min_qual,
            "rare_max_count": rare_max_count,
            "rare_max_freq": rare_max_freq,
            "min_called_for_freq": min_called_for_freq,
            "min_freq_values": min_freq_values,
            "min_background_similarity": min_background_similarity,
            "min_introgression_similarity": (
                min_background_similarity
                if min_introgression_similarity is None
                else min_introgression_similarity
            ),
            "min_introgression_delta": min_introgression_delta,
            "max_introgression_missing_rate": max_introgression_missing_rate,
            "min_introgression_windows": min_introgression_windows,
            "similarity_output": similarity_output,
            "vcf_engine": result.vcf_engine,
            "local_pca": local_pca,
            "write_plots": write_plots,
            "plot_format": plot_format,
        },
        "samples": {
            "full": list(result.groups.full),
            "target": list(result.groups.target),
            "off_target": list(result.groups.off_target),
            "ignored": list(result.groups.ignored),
        },
        "summary": {
            "n_samples": len(result.groups.full),
            "n_windows": result.n_window_rows,
            "n_sample_window_rows": result.n_sample_rows,
            "n_background_blocks": len(result.background_block_rows),
            "n_candidate_introgression_blocks": len(
                result.candidate_introgression_rows
            ),
            "n_pairwise_window_similarity_rows": result.n_similarity_rows,
            "n_local_pca_rows": result.n_local_pca_rows,
            "n_similarity_rows": (
                result.n_similarity_rows
                if similarity_output == "full"
                else len(result.similarity_summary_rows)
                if similarity_output == "summary"
                else 0
            ),
        },
        "outputs": [
            "sample_windows.tsv",
            "windows.tsv",
            "background_blocks.tsv",
            "candidate_introgression_blocks.tsv",
            *(["similarity.tsv"] if similarity_output != "none" else []),
            *(["local_pca.tsv"] if local_pca else []),
            *plot_paths,
            "landscape.json",
        ],
    }
    (outdir / "landscape.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    log.info("Wrote landscape metadata: %s", outdir / "landscape.json")


def run_landscape_plots(
    input_dir: Path,
    plot_format: str = "png",
    outdir: Path | None = None,
    plot_scope: str = "chromosome",
    contigs: list[str] | None = None,
) -> list[Path]:
    """Render landscape plots from existing TSV outputs."""
    if plot_scope not in LANDSCAPE_PLOT_SCOPES:
        raise ValueError("--plot-scope must be one of: chromosome, genome, both.")
    plot_outdir = outdir or input_dir / "plots"
    sample_windows = input_dir / "sample_windows.tsv"
    windows = input_dir / "windows.tsv"
    similarity = input_dir / "similarity.tsv"
    if not sample_windows.exists():
        raise FileNotFoundError(f"Landscape table not found: {sample_windows}")
    if not windows.exists():
        raise FileNotFoundError(f"Landscape table not found: {windows}")

    log.info(
        "Loading existing landscape tables for plotting | input_dir=%s | outdir=%s",
        input_dir,
        plot_outdir,
    )
    sample_rows = read_tsv(sample_windows)
    window_rows = read_tsv(windows)
    similarity_rows = read_tsv(similarity) if similarity.exists() else []

    from privy.plot.landscape import plot_landscape_set  # noqa: PLC0415

    log.info(
        "Rendering landscape plots from existing tables | format=%s | "
        "scope=%s | contigs=%s | sample_rows=%d | windows=%d | similarity_rows=%d",
        plot_format,
        plot_scope,
        ",".join(contigs) if contigs else "all",
        len(sample_rows),
        len(window_rows),
        len(similarity_rows),
    )
    paths, index_rows = plot_landscape_set(
        sample_rows=sample_rows,
        window_rows=window_rows,
        similarity_rows=similarity_rows,
        outdir=plot_outdir,
        output_format=plot_format,
        plot_scope=plot_scope,
        contigs=contigs,
    )
    if index_rows:
        index_path = plot_outdir / "landscape_plot_index.tsv"
        with TsvWriter(index_path, LANDSCAPE_PLOT_INDEX_COLUMNS) as writer:
            writer.write_rows(index_rows)
        paths.append(index_path)

    if _is_relative_to(plot_outdir, input_dir):
        _update_plot_metadata(
            input_dir,
            plot_format,
            paths,
            plot_scope=plot_scope,
            contigs=contigs,
        )
    log.info("Rendered landscape plots from existing tables | count=%d", len(paths))
    return paths


def _update_plot_metadata(
    outdir: Path,
    plot_format: str,
    plot_paths: list[Path],
    plot_scope: str | None = None,
    contigs: list[str] | None = None,
) -> None:
    metadata_path = outdir / "landscape.json"
    if not metadata_path.exists():
        return
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    parameters = metadata.setdefault("parameters", {})
    if isinstance(parameters, dict):
        parameters["write_plots"] = True
        parameters["plot_format"] = plot_format
        if plot_scope is not None:
            parameters["plot_scope"] = plot_scope
        if contigs is not None:
            parameters["plot_contigs"] = list(contigs)
    outputs = metadata.setdefault("outputs", [])
    if isinstance(outputs, list):
        existing = {str(item) for item in outputs}
        for path in plot_paths:
            output_name = _output_name(outdir, path)
            if output_name not in existing:
                outputs.append(output_name)
                existing.add(output_name)
    metadata_path.write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _is_relative_to(path: Path, parent: Path) -> bool:
    try:
        path.resolve().relative_to(parent.resolve())
    except ValueError:
        return False
    return True


def _output_name(outdir: Path, path: Path) -> str:
    try:
        return str(path.resolve().relative_to(outdir.resolve()))
    except ValueError:
        return path.name
