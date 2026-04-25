"""``privy report`` — interpretation and ranking engine.

Converts raw scan and compare outputs into ranked summaries, QC tables,
and human-readable Markdown or HTML reports suitable for collaborators.
"""

from __future__ import annotations

import logging
from pathlib import Path

import typer

from privy.cli.context import get_state
from privy.core.config import PrivyConfig, load_config

log = logging.getLogger("privy.cli.report")

app = typer.Typer(
    name="report",
    help=(
        "Generate ranked summaries and human-readable reports.\n\n"
        "[bold]Outputs:[/bold] summary.tsv, ranked_hits.tsv, "
        "strictness_summary.tsv, support_summary.tsv, "
        "contradiction_summary.tsv, report.md, report.html (optional)"
    ),
    rich_markup_mode="rich",
    no_args_is_help=True,
)


@app.callback(invoke_without_command=True)
def report(
    # --------------------------------------------------------------- inputs
    hits: Path | None = typer.Option(
        None, "--hits", metavar="PATH",
        help="hits.tsv from privy scan.",
    ),
    regions: Path | None = typer.Option(
        None, "--regions", metavar="PATH",
        help="regions.tsv from privy scan.",
    ),
    evidence: Path | None = typer.Option(
        None, "--evidence", metavar="PATH",
        help="evidence.tsv from privy scan.",
    ),
    compare: Path | None = typer.Option(
        None, "--compare", metavar="PATH",
        help="compare.tsv from privy compare.",
    ),
    qc: Path | None = typer.Option(
        None, "--qc", metavar="PATH",
        help="qc.tsv from privy scan.",
    ),
    run_json: Path | None = typer.Option(
        None, "--run-json", metavar="PATH",
        help="run.json from privy scan.",
    ),
    # ---------------------------------------------------------- report options
    fmt: str = typer.Option(
        "markdown", "--format", metavar="TEXT",
        help="Report format: markdown | html | both.",
    ),
    top_n: int = typer.Option(
        20, "--top-n", metavar="INTEGER", min=1,
        help="Number of top loci/regions to summarize.",
    ),
    include_qc: bool = typer.Option(
        True, "--include-qc/--no-include-qc",
        help="Include QC section.",
    ),
    include_strictness: bool = typer.Option(
        True, "--include-strictness/--no-include-strictness",
        help="Include strictness class summary.",
    ),
    include_compare: bool = typer.Option(
        True, "--include-compare/--no-include-compare",
        help="Include compare summary.",
    ),
    include_regions: bool = typer.Option(
        True, "--include-regions/--no-include-regions",
        help="Include candidate region summary.",
    ),
    title: str | None = typer.Option(
        None, "--title", metavar="TEXT",
        help="Optional report title.",
    ),
    outdir: Path | None = typer.Option(
        None, "--outdir", metavar="PATH",
        help="Output directory. Overrides global --outdir.",
    ),
) -> None:
    """Generate ranked summaries and human-readable reports.

    Designed to produce collaborator-ready summaries from privy scan and
    privy compare outputs, including strictness class distributions and
    cross-evidence support summaries.
    """
    from privy.report.summary import run_report  # noqa: PLC0415

    state = get_state()
    effective_outdir = outdir or state.outdir

    cfg: PrivyConfig
    if state.config_path is not None:
        cfg = load_config(state.config_path)
    else:
        from privy.core.config import default_config
        cfg = default_config()

    if hits is None:
        typer.echo("[error] --hits is required.", err=True)
        raise typer.Exit(code=1)

    if not hits.exists():
        typer.echo(f"[error] hits.tsv not found: {hits}", err=True)
        raise typer.Exit(code=1)

    effective_outdir.mkdir(parents=True, exist_ok=True)
    log.info("Starting report | format=%s top_n=%d", fmt, top_n)

    run_report(
        hits=hits,
        regions=regions,
        evidence=evidence,
        compare=compare,
        qc=qc,
        run_json=run_json,
        cfg=cfg,
        fmt=fmt,
        top_n=top_n,
        include_qc=include_qc,
        include_strictness=include_strictness,
        include_compare=include_compare,
        include_regions=include_regions,
        title=title or cfg.project_name,
        outdir=effective_outdir,
    )
