"""``privy plot`` — focused evidence visualization.

Generates publication-quality diagnostic figures from privy scan and
compare outputs.  Designed to explain findings, not to replace a genome
browser.
"""

from __future__ import annotations

import logging
from pathlib import Path

import typer

from privy.cli.context import get_state
from privy.core.config import PrivyConfig, load_config

log = logging.getLogger("privy.cli.plot")

app = typer.Typer(
    name="plot",
    help=(
        "Generate focused diagnostic plots from privy scan and compare outputs.\n\n"
        "[bold]Plot types:[/bold] "
        "locus_panel | strictness_bar | score_distribution | "
        "support_bar | compare_summary | all\n\n"
        "Pass [bold]--plot-type all[/bold] (default) to generate every applicable "
        "figure for the given inputs."
    ),
    rich_markup_mode="rich",
    no_args_is_help=True,
)


@app.callback(invoke_without_command=True)
def plot(
    # ---------------------------------------------------------------- inputs
    hits: Path = typer.Option(
        ..., "--hits", metavar="PATH",
        help="hits.tsv from privy scan (required).",
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
    # --------------------------------------------------------------- outputs
    outdir: Path | None = typer.Option(
        None, "--outdir", metavar="PATH",
        help="Output directory (overrides global --outdir).",
    ),
) -> None:
    """Generate focused diagnostic plots from privy scan and compare outputs."""
    from privy.plot.loci import run_plot  # noqa: PLC0415

    state = get_state()
    effective_outdir = outdir or state.outdir

    cfg: PrivyConfig
    if state.config_path is not None:
        cfg = load_config(state.config_path)
    else:
        from privy.core.config import default_config  # noqa: PLC0415
        cfg = default_config()

    if not hits.exists():
        typer.echo(f"[error] --hits not found: {hits}", err=True)
        raise typer.Exit(code=1)

    effective_outdir.mkdir(parents=True, exist_ok=True)
    log.info("Starting plot | type=%s outdir=%s", plot_type, effective_outdir)

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
