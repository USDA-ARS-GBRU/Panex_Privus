"""``privy dashboard`` — build self-contained interactive comparative dashboards."""

from __future__ import annotations

import logging
from pathlib import Path

import typer

from privy.cli.context import get_state

log = logging.getLogger("privy.cli.dashboard")

app = typer.Typer(
    name="dashboard",
    help=(
        "Build a [bold]self-contained interactive HTML dashboard[/bold] (no server, "
        "offline) from comparative-synteny outputs.\n\n"
        "Point [bold]--synteny[/bold] at a [italic]privy synteny[/italic] output "
        "directory to render linked riparian + dotplot views with target-private "
        "highlighting.\n\n"
        "[bold]Output:[/bold] synteny_dashboard.html"
    ),
    rich_markup_mode="rich",
    no_args_is_help=True,
)


@app.callback(invoke_without_command=True)
def dashboard(
    synteny: Path = typer.Option(
        ..., "--synteny", metavar="PATH",
        help="A `privy synteny` output directory (contains synteny_blocks.tsv).",
    ),
    outdir: Path | None = typer.Option(
        None, "--outdir", metavar="PATH",
        help="Output directory (default: the --synteny directory).",
    ),
) -> None:
    """Build an interactive synteny dashboard from a synteny output directory."""
    from privy.interactive.synteny_dashboard import build_synteny_dashboard  # noqa: PLC0415

    state = get_state()
    if not synteny.exists():
        typer.echo(f"[error] --synteny directory not found: {synteny}", err=True)
        raise typer.Exit(code=1)

    try:
        out_path = build_synteny_dashboard(synteny, outdir=outdir)
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        typer.echo(f"[error] {exc}", err=True)
        raise typer.Exit(code=1) from exc

    log.info("Dashboard complete | %s", out_path)
    if not state.quiet:
        typer.echo(f"  {out_path}")
