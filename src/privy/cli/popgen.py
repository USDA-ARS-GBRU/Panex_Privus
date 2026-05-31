"""``privy popgen`` — breeder population-genetics summaries from a pangenome graph."""

from __future__ import annotations

import logging
from pathlib import Path

import typer

from privy.cli.context import get_state

log = logging.getLogger("privy.cli.popgen")

app = typer.Typer(
    name="popgen",
    help=(
        "[bold]Population-genetics summaries[/bold] for breeders, computed from "
        "multi-allelic microhaplotypes in a pangenome graph.\n\n"
        "Reports per-locus allelic diversity and [bold]target-vs-off-target "
        "differentiation[/bold] (Nei G_ST / Jost's D), surfacing fully diagnostic "
        "(target-private) markers. Cohorts are comma-separated PanSN sample names "
        "or path ids.\n\n"
        "[bold]Outputs:[/bold] popgen_loci.tsv, popgen.json"
    ),
    rich_markup_mode="rich",
    no_args_is_help=True,
)


def _split(value: str | None) -> list[str]:
    if value is None:
        return []
    return [v for v in value.split(",") if v]


@app.callback(invoke_without_command=True)
def popgen(
    gfa: Path = typer.Option(
        ..., "--gfa", metavar="PATH",
        help="Pangenome graph (GFA / GFA.gz).",
    ),
    reference: str = typer.Option(
        ..., "--reference", metavar="PATH_ID",
        help="Path id whose coordinates anchor loci, e.g. sample0#0#chr1.",
    ),
    targets: str = typer.Option(
        ..., "--targets", metavar="NAME,NAME,...",
        help="Target cohort: PanSN sample names or path ids.",
    ),
    off_targets: str = typer.Option(
        ..., "--off-targets", metavar="NAME,NAME,...",
        help="Off-target cohort: PanSN sample names or path ids.",
    ),
    min_core_fraction: float = typer.Option(
        1.0, "--min-core-fraction", metavar="FLOAT", min=0.0, max=1.0,
        help="A segment is 'core' (backbone) when present on >= this fraction of genomes.",
    ),
    outdir: Path | None = typer.Option(
        None, "--outdir", metavar="PATH",
        help="Output directory (overrides global --outdir).",
    ),
) -> None:
    """Compute diversity + cohort differentiation and write breeder tables."""
    from privy.backends.popgen import run_popgen  # noqa: PLC0415

    state = get_state()
    effective_outdir = outdir or state.outdir

    if not gfa.exists():
        typer.echo(f"[error] --gfa not found: {gfa}", err=True)
        raise typer.Exit(code=1)

    try:
        written = run_popgen(
            gfa,
            reference=reference,
            targets=_split(targets),
            off_targets=_split(off_targets),
            min_core_fraction=min_core_fraction,
            outdir=effective_outdir,
        )
    except (ValueError, KeyError) as exc:
        typer.echo(f"[error] {exc}", err=True)
        raise typer.Exit(code=1) from exc

    log.info("Popgen complete | outputs=%d outdir=%s", len(written), effective_outdir)
    if not state.quiet:
        for path in written:
            typer.echo(f"  {path}")
