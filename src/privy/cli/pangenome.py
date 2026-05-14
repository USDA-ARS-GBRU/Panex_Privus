"""``privy pangenome`` — whole and sub-pangenome summaries."""

from __future__ import annotations

import logging
from pathlib import Path

import typer

from privy.cli.cohort_args import (
    collect_sample_values,
    load_cohort_file,
    parse_grouped_cohort_args,
)
from privy.cli.context import get_state

log = logging.getLogger("privy.cli.pangenome")

PANGENOME_HELP = (
    "Analyze full, target, and off-target pangenomes from shared feature "
    "matrices. GFA and VCF inputs use the same cohort syntax.\n\n"
    "[bold]Cohort syntax:[/bold] "
    "--targets SAMPLE [SAMPLE ...] --off-targets SAMPLE [SAMPLE ...]"
)

PANGENOME_CONTEXT_SETTINGS = {
    "allow_extra_args": True,
    "ignore_unknown_options": True,
}

app = typer.Typer(
    name="pangenome",
    help=PANGENOME_HELP,
    context_settings=PANGENOME_CONTEXT_SETTINGS,
    rich_markup_mode="rich",
    no_args_is_help=True,
)


@app.callback(
    invoke_without_command=True,
    context_settings=PANGENOME_CONTEXT_SETTINGS,
)
def pangenome(
    ctx: typer.Context,
    gfa: Path | None = typer.Option(
        None, "--gfa", metavar="PATH",
        help="GFA graph file (.gfa or .gfa.gz) to analyze.",
    ),
    vcf: Path | None = typer.Option(
        None, "--vcf", metavar="PATH",
        help="Multisample VCF/BCF file to analyze.",
    ),
    targets_file: Path | None = typer.Option(
        None, "--targets-file", metavar="PATH",
        help="Text file with one target sample name per line.",
    ),
    off_targets_file: Path | None = typer.Option(
        None, "--off-targets-file", metavar="PATH",
        help="Text file with one off-target sample name per line.",
    ),
    ignore_samples_file: Path | None = typer.Option(
        None, "--ignore-samples-file", metavar="PATH",
        help="Text file with one sample name to ignore per line.",
    ),
    cohort_file: Path | None = typer.Option(
        None, "--cohort-file", metavar="PATH",
        help="Optional cohort definition file (TSV or YAML).",
    ),
    permutations: int = typer.Option(
        100, "--permutations", metavar="INTEGER", min=1,
        help="Number of rarefaction permutations.",
    ),
    seed: int = typer.Option(
        42, "--seed", metavar="INTEGER",
        help="Random seed used for deterministic permutations.",
    ),
    write_plots: bool = typer.Option(
        False, "--plots/--no-plots",
        help=(
            "Write pangenome plots immediately. By default, pangenome writes "
            "data tables only; use privy plot --plot-set pangenome afterwards."
        ),
    ),
    plot_format: str = typer.Option(
        "png", "--plot-format", metavar="TEXT",
        help="Plot format: png, svg, or pdf.",
    ),
    outdir: Path | None = typer.Option(
        None, "--outdir", metavar="PATH",
        help="Output directory for pangenome tables.",
    ),
) -> None:
    """Create pangenome summary tables for full, target, and off-target groups."""
    state = get_state()
    effective_outdir = outdir or state.outdir

    if gfa is None and vcf is None:
        typer.echo("[error] Provide at least one input: --gfa or --vcf.", err=True)
        raise typer.Exit(code=1)
    if gfa is not None and not gfa.exists():
        typer.echo(f"[error] GFA file not found: {gfa}", err=True)
        raise typer.Exit(code=1)
    if vcf is not None and not vcf.exists():
        typer.echo(f"[error] VCF file not found: {vcf}", err=True)
        raise typer.Exit(code=1)

    try:
        cohort_cli_args = parse_grouped_cohort_args(ctx.args)
        cohort_from_file = (
            load_cohort_file(cohort_file) if cohort_file is not None else None
        )
        cli_targets = collect_sample_values(
            cohort_cli_args["targets"],
            targets_file,
        )
        cli_off_targets = collect_sample_values(
            cohort_cli_args["off_targets"],
            off_targets_file,
        )
        cli_ignored = collect_sample_values(
            cohort_cli_args["ignore_samples"],
            ignore_samples_file,
        )
    except (FileNotFoundError, ValueError) as exc:
        typer.echo(f"[error] {exc}", err=True)
        raise typer.Exit(code=2) from exc

    target_list = (
        cli_targets
        if cli_targets is not None
        else (list(cohort_from_file.targets) if cohort_from_file is not None else [])
    )
    off_target_list = (
        cli_off_targets
        if cli_off_targets is not None
        else (
            list(cohort_from_file.off_targets)
            if cohort_from_file is not None
            else []
        )
    )
    ignored_list = (
        cli_ignored
        if cli_ignored is not None
        else (
            list(cohort_from_file.ignored_samples)
            if cohort_from_file is not None
            else []
        )
    )

    if not target_list:
        typer.echo(
            "[error] Provide target samples with --targets, --targets-file, "
            "or --cohort-file.",
            err=True,
        )
        raise typer.Exit(code=1)
    if plot_format not in {"png", "svg", "pdf"}:
        typer.echo("[error] --plot-format must be one of: png, svg, pdf.", err=True)
        raise typer.Exit(code=1)

    try:
        from privy.backends.pangenome import (  # noqa: PLC0415
            run_pangenome_gfa,
            run_pangenome_vcf,
        )

        if gfa is not None and vcf is not None:
            gfa_outdir = effective_outdir / "gfa"
            vcf_outdir = effective_outdir / "vcf"
            run_pangenome_gfa(
                gfa=gfa,
                targets=target_list,
                off_targets=off_target_list or None,
                ignored_samples=ignored_list,
                outdir=gfa_outdir,
                permutations=permutations,
                seed=seed,
                write_plots=write_plots,
                plot_format=plot_format,
            )
            run_pangenome_vcf(
                vcf=vcf,
                targets=target_list,
                off_targets=off_target_list or None,
                ignored_samples=ignored_list,
                outdir=vcf_outdir,
                permutations=permutations,
                seed=seed,
                write_plots=write_plots,
                plot_format=plot_format,
            )
        elif gfa is not None:
            run_pangenome_gfa(
                gfa=gfa,
                targets=target_list,
                off_targets=off_target_list or None,
                ignored_samples=ignored_list,
                outdir=effective_outdir,
                permutations=permutations,
                seed=seed,
                write_plots=write_plots,
                plot_format=plot_format,
            )
        elif vcf is not None:
            run_pangenome_vcf(
                vcf=vcf,
                targets=target_list,
                off_targets=off_target_list or None,
                ignored_samples=ignored_list,
                outdir=effective_outdir,
                permutations=permutations,
                seed=seed,
                write_plots=write_plots,
                plot_format=plot_format,
            )
    except (FileNotFoundError, ValueError) as exc:
        typer.echo(f"[error] {exc}", err=True)
        raise typer.Exit(code=1) from exc
    except NotImplementedError as exc:
        typer.echo(f"[error] Not implemented: {exc}", err=True)
        raise typer.Exit(code=2) from exc

    typer.echo(f"privy pangenome complete. Outputs in: {effective_outdir}")
