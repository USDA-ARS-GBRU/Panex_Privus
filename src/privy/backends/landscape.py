"""Windowed VCF landscape backend."""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any

from privy.io.tsv import TsvWriter
from privy.landscape import (
    BACKGROUND_BLOCK_COLUMNS,
    CANDIDATE_INTROGRESSION_BLOCK_COLUMNS,
    LANDSCAPE_SAMPLE_WINDOW_COLUMNS,
    LANDSCAPE_SIMILARITY_COLUMNS,
    LANDSCAPE_WINDOW_COLUMNS,
    run_vcf_landscape,
)
from privy.utils.misc import now_iso

log = logging.getLogger("privy.backends.landscape")


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
    min_introgression_delta: float = 0.0,
    max_introgression_missing_rate: float = 0.5,
    min_introgression_windows: int = 1,
    write_plots: bool = True,
    plot_format: str = "png",
) -> None:
    """Run a VCF landscape analysis and write tables, plots, and metadata."""
    outdir.mkdir(parents=True, exist_ok=True)
    log.info("Running VCF landscape analysis: %s", vcf)
    log.info("Streaming large landscape tables while analyzing | outdir=%s", outdir)
    with (
        TsvWriter(outdir / "sample_windows.tsv", LANDSCAPE_SAMPLE_WINDOW_COLUMNS) as sample_writer,
        TsvWriter(outdir / "windows.tsv", LANDSCAPE_WINDOW_COLUMNS) as window_writer,
        TsvWriter(outdir / "similarity.tsv", LANDSCAPE_SIMILARITY_COLUMNS) as similarity_writer,
    ):
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
            retain_output_rows=write_plots,
            retain_similarity_rows=False,
        )
    log.info(
        "Landscape analysis complete | windows=%d | sample_window_rows=%d | "
        "background_blocks=%d | candidate_introgression_blocks=%d | "
        "similarity_rows=%d",
        result.n_window_rows,
        result.n_sample_rows,
        len(result.background_block_rows),
        len(result.candidate_introgression_rows),
        result.n_similarity_rows,
    )

    log.info("Writing landscape block tables | outdir=%s", outdir)
    with TsvWriter(outdir / "background_blocks.tsv", BACKGROUND_BLOCK_COLUMNS) as writer:
        writer.write_rows(result.background_block_rows)
    with TsvWriter(
        outdir / "candidate_introgression_blocks.tsv",
        CANDIDATE_INTROGRESSION_BLOCK_COLUMNS,
    ) as writer:
        writer.write_rows(result.candidate_introgression_rows)

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
            "n_similarity_rows": result.n_similarity_rows,
        },
        "outputs": [
            "sample_windows.tsv",
            "windows.tsv",
            "background_blocks.tsv",
            "candidate_introgression_blocks.tsv",
            "similarity.tsv",
            *plot_paths,
            "landscape.json",
        ],
    }
    (outdir / "landscape.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    log.info("Wrote landscape metadata: %s", outdir / "landscape.json")
