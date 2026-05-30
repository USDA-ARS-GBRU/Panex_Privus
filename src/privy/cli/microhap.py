"""``privy microhap`` — detect multi-allelic microhaplotypes from a pangenome graph."""

from __future__ import annotations

import logging
from pathlib import Path

import typer

from privy.cli.context import get_state

log = logging.getLogger("privy.cli.microhap")

app = typer.Typer(
    name="microhap",
    help=(
        "Detect [bold]microhaplotypes[/bold] — local multi-allelic loci — from a "
        "pangenome graph, and flag [bold]target-private[/bold] alleles.\n\n"
        "Each genome's local sequence between shared flanks is its allele "
        "(content-hashed, PHG-style). Provide --targets / --off-targets (comma-"
        "separated PanSN sample names or path ids) to flag alleles present in the "
        "target cohort and absent from off-targets.\n\n"
        "[bold]Outputs:[/bold] microhaplotypes.tsv, allele_matrix.tsv, microhap.json"
    ),
    rich_markup_mode="rich",
    no_args_is_help=True,
)


def _split(value: str | None) -> list[str] | None:
    if value is None:
        return None
    return [v for v in value.split(",") if v]


@app.callback(invoke_without_command=True)
def microhap(
    gfa: Path = typer.Option(
        ..., "--gfa", metavar="PATH",
        help="Pangenome graph (GFA / GFA.gz).",
    ),
    reference: str = typer.Option(
        ..., "--reference", metavar="PATH_ID",
        help="Path id whose coordinates anchor loci, e.g. sample0#0#chr1.",
    ),
    targets: str | None = typer.Option(
        None, "--targets", metavar="NAME,NAME,...",
        help="Target cohort: PanSN sample names or path ids (for private-allele flagging).",
    ),
    off_targets: str | None = typer.Option(
        None, "--off-targets", metavar="NAME,NAME,...",
        help="Off-target cohort: PanSN sample names or path ids.",
    ),
    min_core_fraction: float = typer.Option(
        1.0, "--min-core-fraction", metavar="FLOAT", min=0.0, max=1.0,
        help="A segment is 'core' (backbone) when present on >= this fraction of genomes.",
    ),
    all_loci: bool = typer.Option(
        False, "--all-loci/--multiallelic-only",
        help="Emit all loci, including monomorphic ones (default: multi-allelic only).",
    ),
    outdir: Path | None = typer.Option(
        None, "--outdir", metavar="PATH",
        help="Output directory (overrides global --outdir).",
    ),
) -> None:
    """Detect microhaplotypes and write allele tables."""
    from privy.backends.microhap import run_microhap  # noqa: PLC0415

    state = get_state()
    effective_outdir = outdir or state.outdir

    if not gfa.exists():
        typer.echo(f"[error] --gfa not found: {gfa}", err=True)
        raise typer.Exit(code=1)

    try:
        written = run_microhap(
            gfa,
            reference=reference,
            targets=_split(targets),
            off_targets=_split(off_targets),
            min_core_fraction=min_core_fraction,
            multiallelic_only=not all_loci,
            outdir=effective_outdir,
        )
    except (ValueError, KeyError) as exc:
        typer.echo(f"[error] {exc}", err=True)
        raise typer.Exit(code=1) from exc

    log.info("Microhap complete | outputs=%d outdir=%s", len(written), effective_outdir)
    if not state.quiet:
        for path in written:
            typer.echo(f"  {path}")
