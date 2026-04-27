"""Main CLI entry point for Panex Privus (``privy``).

Registers all top-level subcommands and handles global options that apply
across every subcommand (config path, output directory, thread count, etc.).

Usage::

    privy --help
    privy scan --vcf cohort.vcf.gz --targets S1 S2 --off-targets S3 S4 --outdir results/
    privy compare   --hits-a results/vcf/hits.tsv --hits-b results/gfa/hits.tsv \
                    --outdir results/compare/
    privy report    --hits results/vcf/hits.tsv --outdir report/
    privy plot      --hits results/vcf/hits.tsv --top-n 10 --outdir plots/
    privy annotate  --hits results/vcf/hits.tsv --gff annotation.gff3.gz --outdir annotated/
    privy export    --hits results/vcf/hits.tsv --regions results/vcf/regions.tsv --outdir exported/
"""

from __future__ import annotations

from pathlib import Path

import typer

from privy import __version__
from privy.cli import annotate, compare, export, pangenome, plot, report, scan
from privy.cli.context import get_state
from privy.utils.logging import configure_logging

app = typer.Typer(
    name="privy",
    help=(
        "[bold]Panex Privus[/bold] — a comparative genomics toolkit for discovering\n"
        "target-private alleles and regions shared within a focal cohort and absent\n"
        "from off-target genomes.\n\n"
        "Primary discovery runs via VCF or GFA. BAM provides a read-depth support layer.\n"
        "Missingness is always reported explicitly via [italic]strictness_class[/italic]."
    ),
    add_completion=True,
    rich_markup_mode="rich",
    no_args_is_help=True,
)

app.add_typer(scan.app, name="scan")
app.add_typer(compare.app, name="compare")
app.add_typer(pangenome.app, name="pangenome")
app.add_typer(report.app, name="report")
app.add_typer(plot.app, name="plot")
app.add_typer(annotate.app, name="annotate")
app.add_typer(export.app, name="export")


def _version_callback(value: bool) -> None:
    if value:
        typer.echo(f"privy {__version__}")
        raise typer.Exit()


@app.callback()
def global_options(
    config: Path | None = typer.Option(
        None,
        "--config",
        metavar="PATH",
        help="Path to YAML configuration file.",
        show_default=False,
    ),
    project_name: str | None = typer.Option(
        None,
        "--project-name",
        metavar="TEXT",
        help="Optional project name written into outputs.",
        show_default=False,
    ),
    outdir: Path = typer.Option(
        Path("."),
        "--outdir",
        metavar="PATH",
        help="Default output directory (may be overridden per subcommand).",
    ),
    threads: int = typer.Option(
        1,
        "--threads",
        metavar="INTEGER",
        min=1,
        help="Worker threads where supported; current scan backends run serially.",
    ),
    log_level: str = typer.Option(
        "info",
        "--log-level",
        metavar="TEXT",
        help="Logging level: debug, info, warning, error.",
    ),
    quiet: bool = typer.Option(
        False,
        "--quiet",
        help="Reduce console output.",
    ),
    version: bool | None = typer.Option(
        None,
        "--version",
        callback=_version_callback,
        is_eager=True,
        help="Show version and exit.",
    ),
) -> None:
    """Panex Privus — target-private genomic signal discovery toolkit."""
    state = get_state()
    state.config_path = config
    state.project_name = project_name
    state.outdir = outdir
    state.threads = threads
    state.log_level = log_level
    state.quiet = quiet
    configure_logging(level=log_level, quiet=quiet)


if __name__ == "__main__":
    app()
