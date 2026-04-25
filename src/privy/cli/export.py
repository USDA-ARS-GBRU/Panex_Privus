"""``privy export`` — write scan results to downstream genome-tool formats."""

from __future__ import annotations

import logging
from pathlib import Path

import typer

from privy.cli.context import get_state

log = logging.getLogger("privy.cli.export")

app = typer.Typer(
    name="export",
    help=(
        "Export [bold]privy scan[/bold] outputs to downstream genome-tool formats.\n\n"
        "Supported formats: BED and GFF3.\n\n"
        "[bold]Outputs:[/bold] hits.bed/gff3, regions.bed/gff3, export.json"
    ),
    rich_markup_mode="rich",
    no_args_is_help=True,
)


@app.callback(invoke_without_command=True)
def export(
    # ---------------------------------------------------------------- inputs
    hits: Path | None = typer.Option(
        None, "--hits", metavar="PATH",
        help="hits.tsv produced by privy scan.",
    ),
    regions: Path | None = typer.Option(
        None, "--regions", metavar="PATH",
        help="regions.tsv produced by privy scan.",
    ),
    # --------------------------------------------------------------- options
    fmt: str = typer.Option(
        "bed", "--format", metavar="TEXT",
        help="Export format: bed | gff3.",
    ),
    kind: str = typer.Option(
        "both", "--kind", metavar="TEXT",
        help="What to export: hits | regions | both.",
    ),
    track_name: str = typer.Option(
        "Panex Privus", "--track-name", metavar="TEXT",
        help="BED track name when writing track headers.",
    ),
    include_header: bool = typer.Option(
        True, "--include-header/--no-include-header",
        help="Include a UCSC-style BED track header.",
    ),
    # --------------------------------------------------------------- outputs
    outdir: Path | None = typer.Option(
        None, "--outdir", metavar="PATH",
        help="Output directory (overrides global --outdir).",
    ),
) -> None:
    """Export scan result tables to BED or GFF3."""
    from privy.backends.export import ExportFormat, ExportKind, run_export  # noqa: PLC0415

    state = get_state()
    effective_outdir = outdir or state.outdir

    if hits is not None and not hits.exists():
        typer.echo(f"[error] --hits not found: {hits}", err=True)
        raise typer.Exit(code=1)
    if regions is not None and not regions.exists():
        typer.echo(f"[error] --regions not found: {regions}", err=True)
        raise typer.Exit(code=1)

    if fmt == "bed":
        export_format: ExportFormat = "bed"
    elif fmt == "gff3":
        export_format = "gff3"
    else:
        typer.echo(
            f"[error] Unsupported export format: {fmt!r}. Use 'bed' or 'gff3'.",
            err=True,
        )
        raise typer.Exit(code=1)
    if kind == "hits":
        export_kind: ExportKind = "hits"
    elif kind == "regions":
        export_kind = "regions"
    elif kind == "both":
        export_kind = "both"
    else:
        typer.echo(
            f"[error] Unsupported export kind: {kind!r}. Use 'hits', 'regions', or 'both'.",
            err=True,
        )
        raise typer.Exit(code=1)

    try:
        written = run_export(
            hits_path=hits,
            regions_path=regions,
            outdir=effective_outdir,
            export_format=export_format,
            export_kind=export_kind,
            track_name=track_name,
            include_header=include_header,
        )
    except ValueError as exc:
        typer.echo(f"[error] {exc}", err=True)
        raise typer.Exit(code=1) from exc

    log.info("Export complete | outputs=%d outdir=%s", len(written), effective_outdir)
    if not state.quiet:
        for path in written:
            typer.echo(f"  {path}")
