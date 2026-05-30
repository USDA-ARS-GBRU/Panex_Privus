"""``privy synteny`` — derive graph-native synteny blocks and private regions."""

from __future__ import annotations

import logging
from pathlib import Path

import typer

from privy.cli.context import get_state

log = logging.getLogger("privy.cli.synteny")

app = typer.Typer(
    name="synteny",
    help=(
        "Derive [bold]graph-native synteny[/bold] from a pangenome graph: collinear "
        "blocks plus typed rearrangements (inversion / translocation / duplication), "
        "and flag [bold]target-private[/bold] structural regions.\n\n"
        "Provide --targets / --off-targets (comma-separated PanSN sample names or path "
        "ids) to tag regions present in the target cohort and absent from off-targets.\n\n"
        "[bold]Outputs:[/bold] synteny_blocks.tsv, synteny_regions.tsv, synteny.json"
    ),
    rich_markup_mode="rich",
    no_args_is_help=True,
)


def _split(value: str | None) -> list[str] | None:
    if value is None:
        return None
    return [v for v in value.split(",") if v]


@app.callback(invoke_without_command=True)
def synteny(
    gfa: Path = typer.Option(
        ..., "--gfa", metavar="PATH",
        help="Pangenome graph (GFA / GFA.gz).",
    ),
    reference: str = typer.Option(
        ..., "--reference", metavar="PATH_ID",
        help="Path id used as the coordinate reference, e.g. sample0#0#chr1.",
    ),
    queries: str | None = typer.Option(
        None, "--queries", metavar="PATH,PATH,...",
        help="Query path ids to compare (default: all paths except the reference).",
    ),
    targets: str | None = typer.Option(
        None, "--targets", metavar="NAME,NAME,...",
        help="Target cohort: PanSN sample names or path ids (for private-region tagging).",
    ),
    off_targets: str | None = typer.Option(
        None, "--off-targets", metavar="NAME,NAME,...",
        help="Off-target cohort: PanSN sample names or path ids.",
    ),
    min_block_anchors: int = typer.Option(
        1, "--min-block-anchors", metavar="INTEGER", min=1,
        help="Minimum anchors for a reported collinear/inversion block.",
    ),
    outdir: Path | None = typer.Option(
        None, "--outdir", metavar="PATH",
        help="Output directory (overrides global --outdir).",
    ),
) -> None:
    """Build synteny blocks/regions from a pangenome graph."""
    from privy.backends.synteny import run_synteny  # noqa: PLC0415

    state = get_state()
    effective_outdir = outdir or state.outdir

    if not gfa.exists():
        typer.echo(f"[error] --gfa not found: {gfa}", err=True)
        raise typer.Exit(code=1)

    try:
        written = run_synteny(
            gfa,
            reference=reference,
            query_paths=_split(queries),
            targets=_split(targets),
            off_targets=_split(off_targets),
            min_block_anchors=min_block_anchors,
            outdir=effective_outdir,
        )
    except (ValueError, KeyError) as exc:
        typer.echo(f"[error] {exc}", err=True)
        raise typer.Exit(code=1) from exc

    log.info("Synteny complete | outputs=%d outdir=%s", len(written), effective_outdir)
    if not state.quiet:
        for path in written:
            typer.echo(f"  {path}")
