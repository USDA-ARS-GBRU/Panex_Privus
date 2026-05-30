"""``privy project`` — project a region or graph node-set onto any reference(s)."""

from __future__ import annotations

import logging
from pathlib import Path

import typer

from privy.cli.context import get_state

log = logging.getLogger("privy.cli.project")

app = typer.Typer(
    name="project",
    help=(
        "Project a region — or a raw set of graph segments — onto [bold]any[/bold] "
        "reference genome embedded in a pangenome graph.\n\n"
        "Give a source region as [italic]PATH:START-END[/italic] (stable coordinates "
        "on that path), or define it in graph node space with --node-set.\n\n"
        "[bold]Outputs:[/bold] projection.tsv, project.json"
    ),
    rich_markup_mode="rich",
    no_args_is_help=True,
)


def _parse_region(region: str) -> tuple[str, int, int]:
    """Parse ``PATH:START-END`` into ``(path, start, end)`` (0-based half-open)."""
    if ":" not in region:
        raise ValueError("region must be PATH:START-END, e.g. 'sample0#0#chr1:100-500'")
    path, span = region.rsplit(":", 1)
    if "-" not in span:
        raise ValueError("region span must be START-END, e.g. '100-500'")
    start_s, end_s = span.split("-", 1)
    try:
        start, end = int(start_s), int(end_s)
    except ValueError as exc:
        raise ValueError(f"region coordinates must be integers: {span!r}") from exc
    if not path:
        raise ValueError("region must include a source path before ':'")
    return path, start, end


@app.callback(invoke_without_command=True)
def project(
    gfa: Path = typer.Option(
        ..., "--gfa", metavar="PATH",
        help="Pangenome graph (GFA / GFA.gz) to project within.",
    ),
    region: str | None = typer.Option(
        None, "--region", metavar="PATH:START-END",
        help="Source region in a path's stable coordinates, e.g. sample0#0#chr1:100-500.",
    ),
    node_set: str | None = typer.Option(
        None, "--node-set", metavar="SEG,SEG,...",
        help="Define the region in graph node space as a comma-separated segment list.",
    ),
    to_genomes: str | None = typer.Option(
        None, "--to-genomes", metavar="PATH,PATH,...",
        help="Comma-separated target path ids (default: all paths in the graph).",
    ),
    outdir: Path | None = typer.Option(
        None, "--outdir", metavar="PATH",
        help="Output directory (overrides global --outdir).",
    ),
) -> None:
    """Project a region or node-set across a pangenome graph to chosen references."""
    from privy.backends.project import run_project  # noqa: PLC0415

    state = get_state()
    effective_outdir = outdir or state.outdir

    if not gfa.exists():
        typer.echo(f"[error] --gfa not found: {gfa}", err=True)
        raise typer.Exit(code=1)
    if (region is None) == (node_set is None):
        typer.echo(
            "[error] provide exactly one of --region or --node-set.", err=True
        )
        raise typer.Exit(code=1)

    source_path: str | None = None
    span: tuple[int, int] | None = None
    segments: list[str] | None = None
    if region is not None:
        try:
            source_path, start, end = _parse_region(region)
        except ValueError as exc:
            typer.echo(f"[error] {exc}", err=True)
            raise typer.Exit(code=1) from exc
        span = (start, end)
    else:
        assert node_set is not None
        segments = [s for s in node_set.split(",") if s]
        if not segments:
            typer.echo("[error] --node-set is empty.", err=True)
            raise typer.Exit(code=1)

    targets = [t for t in to_genomes.split(",") if t] if to_genomes else None

    try:
        written = run_project(
            gfa,
            source_path=source_path,
            region=span,
            node_set=segments,
            to_genomes=targets,
            outdir=effective_outdir,
        )
    except (ValueError, KeyError, IndexError) as exc:
        typer.echo(f"[error] {exc}", err=True)
        raise typer.Exit(code=1) from exc

    log.info("Project complete | outputs=%d outdir=%s", len(written), effective_outdir)
    if not state.quiet:
        for path in written:
            typer.echo(f"  {path}")
