"""``privy annotate`` — intersect private loci with GFF3 gene annotations.

Reads a hits.tsv from ``privy scan`` and a GFF3 annotation file, classifies
each hit as CDS / UTR / exonic / intronic / intergenic, and writes:

  - ``annotated_hits.tsv`` — all hits columns plus annotation columns
  - ``annotation_summary.tsv`` — class counts and percentages
  - ``annotate.json`` — run metadata
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import typer

from privy.cli.context import get_state

log = logging.getLogger("privy.cli.annotate")

app = typer.Typer(
    name="annotate",
    help=(
        "Intersect private loci from [bold]privy scan[/bold] with a GFF3 gene "
        "annotation.\n\n"
        "Each hit is classified using a feature hierarchy: "
        "[italic]CDS → UTR → exonic → intronic → intergenic[/italic].\n\n"
        "[bold]Outputs:[/bold] annotated_hits.tsv, annotation_summary.tsv, annotate.json"
    ),
    rich_markup_mode="rich",
    no_args_is_help=True,
)


@app.callback(invoke_without_command=True)
def annotate(
    # ---------------------------------------------------------------- inputs
    hits: Path = typer.Option(
        ..., "--hits", metavar="PATH",
        help="hits.tsv produced by privy scan.",
    ),
    gff: Path = typer.Option(
        ..., "--gff", metavar="PATH",
        help="GFF3 annotation file (plain or .gz).",
    ),
    # -------------------------------------------------- contig name handling
    contig_alias: Optional[Path] = typer.Option(
        None, "--contig-alias", metavar="PATH",
        help=(
            "Two-column TSV mapping contig names.  Default direction: "
            "GFF3 name → hits name.  Flip with --hits-to-gff."
        ),
    ),
    hits_to_gff: bool = typer.Option(
        False, "--hits-to-gff/--gff-to-hits",
        help=(
            "Direction of --contig-alias.  "
            "--hits-to-gff: hits name → GFF3 name (default: GFF3 → hits)."
        ),
    ),
    # --------------------------------------------------------------- outputs
    outdir: Optional[Path] = typer.Option(
        None, "--outdir", metavar="PATH",
        help="Output directory (overrides global --outdir).",
    ),
) -> None:
    """Annotate private loci with GFF3 gene features."""
    from privy.backends.annotate import run_annotate  # noqa: PLC0415

    state = get_state()
    effective_outdir = outdir or state.outdir

    if not hits.exists():
        typer.echo(f"[error] --hits not found: {hits}", err=True)
        raise typer.Exit(code=1)
    if not gff.exists():
        typer.echo(f"[error] --gff not found: {gff}", err=True)
        raise typer.Exit(code=1)
    if contig_alias is not None and not contig_alias.exists():
        typer.echo(f"[error] --contig-alias not found: {contig_alias}", err=True)
        raise typer.Exit(code=1)

    log.info("Starting annotate | hits=%s gff=%s", hits, gff)

    run_annotate(
        hits_path=hits,
        gff_path=gff,
        outdir=effective_outdir,
        contig_alias_path=contig_alias,
        hits_contig_to_gff=hits_to_gff,
    )

    typer.echo(
        f"privy annotate complete. Outputs in: {effective_outdir}"
    )
