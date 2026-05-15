"""``privy plot`` — focused evidence visualization.

Generates publication-quality diagnostic figures from existing Privy output
tables.  Designed to explain findings, not to replace a genome browser.
"""

from __future__ import annotations

import logging
from pathlib import Path

import typer

from privy.cli.context import get_state
from privy.core.config import PrivyConfig, load_config

log = logging.getLogger("privy.cli.plot")

PLOT_SETS = {"scan", "landscape", "pangenome"}

app = typer.Typer(
    name="plot",
    help=(
        "Generate focused diagnostic plots from existing Privy outputs.\n\n"
        "[bold]Plot sets:[/bold] scan | landscape | pangenome\n\n"
        "[bold]Scan plot types:[/bold] "
        "locus_panel | strictness_bar | score_distribution | "
        "support_bar | compare_summary | all\n\n"
        "Use [bold]--plot-set scan[/bold] for the current scan/compare plots, "
        "or point [bold]--input-dir[/bold] at a landscape or pangenome result "
        "directory."
    ),
    rich_markup_mode="rich",
    no_args_is_help=True,
)


@app.callback(invoke_without_command=True)
def plot(
    # ---------------------------------------------------------------- inputs
    plot_set: str = typer.Option(
        "scan", "--plot-set", metavar="TEXT",
        help="Plot set to render: scan, landscape, or pangenome.",
    ),
    input_dir: Path | None = typer.Option(
        None, "--input-dir", metavar="PATH",
        help="Existing landscape or pangenome result directory to plot.",
    ),
    hits: Path | None = typer.Option(
        None, "--hits", metavar="PATH",
        help="hits.tsv from privy scan (required for --plot-set scan).",
    ),
    regions: Path | None = typer.Option(
        None, "--regions", metavar="PATH",
        help="regions.tsv from privy scan (optional).",
    ),
    evidence: Path | None = typer.Option(
        None, "--evidence", metavar="PATH",
        help="evidence.tsv from privy scan (enables support_bar plot).",
    ),
    compare: Path | None = typer.Option(
        None, "--compare", metavar="PATH",
        help="compare.tsv from privy compare (enables compare_summary plot).",
    ),
    # ----------------------------------------------------------- plot type
    plot_type: str = typer.Option(
        "all", "--plot-type", metavar="TEXT",
        help=(
            "Which plot to generate: "
            "locus_panel | strictness_bar | score_distribution | "
            "support_bar | compare_summary | all."
        ),
    ),
    # ------------------------------------------------------------- selection
    top_n: int | None = typer.Option(
        None, "--top-n", metavar="INTEGER", min=1,
        help="Number of top loci to show in locus_panel [default: 30].",
    ),
    show_labels: bool = typer.Option(
        True, "--show-labels/--no-show-labels",
        help="Annotate x-axis with locus IDs in locus_panel (only when top_n ≤ 20).",
    ),
    # --------------------------------------------------------- figure options
    width: float = typer.Option(
        10.0, "--width", metavar="FLOAT",
        help="Figure width in inches.",
    ),
    height: float = typer.Option(
        5.0, "--height", metavar="FLOAT",
        help="Figure height in inches.",
    ),
    dpi: int = typer.Option(
        150, "--dpi", metavar="INTEGER", min=72,
        help="Figure DPI (raster formats).",
    ),
    output_format: str = typer.Option(
        "png", "--output-format", metavar="TEXT",
        help="Output format: png | svg | pdf.",
    ),
    plot_scope: str = typer.Option(
        "chromosome", "--plot-scope", metavar="TEXT",
        help="Landscape plot scope: chromosome | genome | both.",
    ),
    contig: str | None = typer.Option(
        None, "--contig", metavar="TEXT",
        help="Landscape contig/chromosome to plot.",
    ),
    contigs: str | None = typer.Option(
        None, "--contigs", metavar="TEXT",
        help="Comma-separated landscape contigs/chromosomes to plot.",
    ),
    # --------------------------------------------------------------- outputs
    outdir: Path | None = typer.Option(
        None, "--outdir", metavar="PATH",
        help="Output directory (overrides global --outdir).",
    ),
) -> None:
    """Generate focused diagnostic plots from existing Privy outputs."""
    state = get_state()
    effective_outdir = outdir or state.outdir
    normalized_plot_set = plot_set.lower()

    if normalized_plot_set not in PLOT_SETS:
        typer.echo(
            "[error] --plot-set must be one of: scan, landscape, pangenome.",
            err=True,
        )
        raise typer.Exit(code=1)
    if output_format not in {"png", "svg", "pdf"}:
        typer.echo("[error] --output-format must be one of: png, svg, pdf.", err=True)
        raise typer.Exit(code=1)

    if normalized_plot_set == "landscape":
        if plot_scope not in {"chromosome", "genome", "both"}:
            typer.echo(
                "[error] --plot-scope must be one of: chromosome, genome, both.",
                err=True,
            )
            raise typer.Exit(code=1)
        source_dir = input_dir or effective_outdir
        plot_outdir = outdir if input_dir is not None and outdir is not None else None
        try:
            from privy.backends.landscape import run_landscape_plots  # noqa: PLC0415

            generated = run_landscape_plots(
                input_dir=source_dir,
                plot_format=output_format,
                outdir=plot_outdir,
                plot_scope=plot_scope,
                contigs=_parse_contigs(contig, contigs),
            )
        except (FileNotFoundError, ValueError) as exc:
            typer.echo(f"[error] {exc}", err=True)
            raise typer.Exit(code=1) from exc
        if not state.quiet:
            for path in generated:
                typer.echo(f"  {path}")
        return

    if contig is not None or contigs is not None or plot_scope != "chromosome":
        typer.echo(
            "[error] --plot-scope, --contig, and --contigs are only used for "
            "landscape plots.",
            err=True,
        )
        raise typer.Exit(code=1)

    if normalized_plot_set == "pangenome":
        source_dir = input_dir or effective_outdir
        plot_outdir = outdir if input_dir is not None and outdir is not None else None
        try:
            from privy.backends.pangenome import run_pangenome_plots  # noqa: PLC0415

            generated = run_pangenome_plots(
                input_dir=source_dir,
                plot_format=output_format,
                outdir=plot_outdir,
            )
        except (FileNotFoundError, ValueError) as exc:
            typer.echo(f"[error] {exc}", err=True)
            raise typer.Exit(code=1) from exc
        if not state.quiet:
            for path in generated:
                typer.echo(f"  {path}")
        return

    if input_dir is not None:
        typer.echo("[error] --input-dir is only used for landscape or pangenome plots.",
                   err=True)
        raise typer.Exit(code=1)
    if hits is None:
        typer.echo("[error] --hits is required for --plot-set scan.", err=True)
        raise typer.Exit(code=1)
    if not hits.exists():
        typer.echo(f"[error] --hits not found: {hits}", err=True)
        raise typer.Exit(code=1)

    from privy.plot.loci import run_plot  # noqa: PLC0415

    cfg: PrivyConfig
    if state.config_path is not None:
        cfg = load_config(state.config_path)
    else:
        from privy.core.config import default_config  # noqa: PLC0415
        cfg = default_config()

    effective_outdir.mkdir(parents=True, exist_ok=True)
    log.info(
        "Starting plot | set=scan | type=%s | outdir=%s",
        plot_type,
        effective_outdir,
    )

    generated = run_plot(
        hits=hits,
        regions=regions,
        evidence=evidence,
        vcf=None,
        bam=None,
        bam_manifest=None,
        gfa=None,
        xmfa=None,
        compare=compare,
        cfg=cfg,
        locus_id=None,
        region_id=None,
        top_n=top_n,
        contig=None,
        region=None,
        plot_type=plot_type,
        width=width,
        height=height,
        dpi=dpi,
        output_format=output_format,
        show_labels=show_labels,
        outdir=effective_outdir,
    )

    if not state.quiet:
        for path in generated:
            typer.echo(f"  {path}")


def _parse_contigs(contig: str | None, contigs: str | None) -> list[str] | None:
    values: list[str] = []
    if contig:
        values.append(contig)
    if contigs:
        values.extend(part.strip() for part in contigs.split(",") if part.strip())
    return values or None
