"""``privy plot`` — focused evidence visualization engine.

Generates publication-quality plots for loci, regions, and summary diagnostics.
Designed to explain and diagnose findings — not to be a genome browser.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import List, Optional

import typer

from privy.cli.context import get_state
from privy.core.config import PrivyConfig, load_config

log = logging.getLogger("privy.cli.plot")

app = typer.Typer(
    name="plot",
    help=(
        "Generate focused plots for loci, regions, and summary diagnostics.\n\n"
        "[bold]Plot types:[/bold] locus_panel | region_summary | "
        "genotype_heatmap | strictness_bar | support_bar | depth_panel\n\n"
        "Focused explanatory plots, not a genome browser."
    ),
    rich_markup_mode="rich",
    no_args_is_help=True,
)


@app.callback(invoke_without_command=True)
def plot(
    # --------------------------------------------------------------- inputs
    hits: Optional[Path] = typer.Option(
        None, "--hits", metavar="PATH",
        help="hits.tsv from privy scan.",
    ),
    regions: Optional[Path] = typer.Option(
        None, "--regions", metavar="PATH",
        help="regions.tsv from privy scan.",
    ),
    evidence: Optional[Path] = typer.Option(
        None, "--evidence", metavar="PATH",
        help="evidence.tsv from privy scan.",
    ),
    vcf: Optional[Path] = typer.Option(
        None, "--vcf", metavar="PATH",
        help="Indexed VCF for genotype visualization.",
    ),
    bam: Optional[List[Path]] = typer.Option(
        None, "--bam", metavar="PATH",
        help="BAM files for depth panels. Repeat flag for multiple files.",
    ),
    bam_manifest: Optional[Path] = typer.Option(
        None, "--bam-manifest", metavar="PATH",
        help="BAM manifest TSV.",
    ),
    gfa: Optional[Path] = typer.Option(
        None, "--gfa", metavar="PATH",
        help="GFA graph file for context strips.",
    ),
    xmfa: Optional[Path] = typer.Option(
        None, "--xmfa", metavar="PATH",
        help="XMFA alignment file for corroboration strips.",
    ),
    # ------------------------------------------------------------- selection
    locus_id: Optional[str] = typer.Option(
        None, "--locus-id", metavar="TEXT",
        help="Plot a specific locus by ID.",
    ),
    region_id: Optional[str] = typer.Option(
        None, "--region-id", metavar="TEXT",
        help="Plot a specific region by ID.",
    ),
    top_n: Optional[int] = typer.Option(
        None, "--top-n", metavar="INTEGER", min=1,
        help="Plot top N loci or regions by final_score.",
    ),
    contig: Optional[str] = typer.Option(
        None, "--contig", metavar="TEXT",
        help="Restrict plots to a contig.",
    ),
    region: Optional[str] = typer.Option(
        None, "--region", metavar="TEXT",
        help="Restrict plots to contig:start-end.",
    ),
    # ----------------------------------------------------------- plot types
    plot_type: str = typer.Option(
        "locus_panel", "--plot-type", metavar="TEXT",
        help=(
            "Plot type: locus_panel | region_summary | genotype_heatmap | "
            "strictness_bar | support_bar | depth_panel."
        ),
    ),
    # --------------------------------------------------------- plot options
    width: float = typer.Option(
        12.0, "--width", metavar="FLOAT",
        help="Figure width in inches.",
    ),
    height: float = typer.Option(
        6.0, "--height", metavar="FLOAT",
        help="Figure height in inches.",
    ),
    dpi: int = typer.Option(
        150, "--dpi", metavar="INTEGER", min=72,
        help="Figure DPI.",
    ),
    output_format: str = typer.Option(
        "png", "--output-format", metavar="TEXT",
        help="Output format: png | pdf | svg.",
    ),
    show_labels: bool = typer.Option(
        True, "--show-labels/--no-show-labels",
        help="Show sample or locus labels where applicable.",
    ),
    outdir: Optional[Path] = typer.Option(
        None, "--outdir", metavar="PATH",
        help="Output directory. Overrides global --outdir.",
    ),
) -> None:
    """Generate focused plots for loci, regions, and summary diagnostics.

    Designed for diagnostics and publication-ready figure generation.
    Not a genome browser — each plot type explains a specific finding.
    """
    from privy.plot.loci import run_plot  # noqa: PLC0415

    state = get_state()
    effective_outdir = outdir or state.outdir

    cfg: PrivyConfig
    if state.config_path is not None:
        cfg = load_config(state.config_path)
    else:
        from privy.core.config import default_config
        cfg = default_config()

    if hits is None:
        typer.echo("[error] --hits is required.", err=True)
        raise typer.Exit(code=1)

    if not hits.exists():
        typer.echo(f"[error] hits.tsv not found: {hits}", err=True)
        raise typer.Exit(code=1)

    effective_outdir.mkdir(parents=True, exist_ok=True)
    log.info("Starting plot | type=%s", plot_type)

    run_plot(
        hits=hits,
        regions=regions,
        evidence=evidence,
        vcf=vcf,
        bam=bam,
        bam_manifest=bam_manifest,
        gfa=gfa,
        xmfa=xmfa,
        cfg=cfg,
        locus_id=locus_id,
        region_id=region_id,
        top_n=top_n,
        contig=contig,
        region=region,
        plot_type=plot_type,
        width=width,
        height=height,
        dpi=dpi,
        output_format=output_format,
        show_labels=show_labels,
        outdir=effective_outdir,
    )
