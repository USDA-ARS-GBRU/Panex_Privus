"""``privy index`` — build reusable indexes for large primary inputs."""

from __future__ import annotations

import logging
from pathlib import Path

import typer

from privy.io.gfa import (
    build_gfa_scan_index,
    default_gfa_index_path,
    write_gfa_scan_index,
)

log = logging.getLogger("privy.cli.index")

app = typer.Typer(
    name="index",
    help=(
        "Build reusable input indexes. GFA indexes are optional but recommended "
        "for large pangenome graphs because they let later scans skip the slow "
        "GFA walk-parsing step."
    ),
    rich_markup_mode="rich",
    no_args_is_help=True,
)


@app.command("gfa")
def index_gfa(
    gfa: Path = typer.Option(
        ...,
        "--gfa",
        metavar="PATH",
        help="Input GFA graph file (.gfa or .gfa.gz).",
    ),
    out: Path | None = typer.Option(
        None,
        "--out",
        metavar="PATH",
        help="Output index path. Defaults to <GFA>.privy.gfaidx.",
    ),
    force: bool = typer.Option(
        False,
        "--force",
        help="Overwrite an existing GFA index.",
    ),
) -> None:
    """Build a reusable Privy GFA scan index."""
    if not gfa.exists():
        typer.echo(f"[error] GFA file not found: {gfa}", err=True)
        raise typer.Exit(code=1)

    index_path = out or default_gfa_index_path(gfa)
    if index_path.exists() and not force:
        typer.echo(
            f"[error] GFA index already exists: {index_path}\n"
            "Use --force to rebuild it.",
            err=True,
        )
        raise typer.Exit(code=1)

    log.info("Building GFA index | gfa=%s | out=%s", gfa, index_path)
    scan_index = build_gfa_scan_index(gfa_path=gfa, sample_names=None)
    write_gfa_scan_index(scan_index=scan_index, index_path=index_path, gfa_path=gfa)
    typer.echo(
        f"Wrote {index_path} "
        f"({len(scan_index.segments)} coordinate segments, "
        f"{len(scan_index.sample_order)} samples)"
    )
