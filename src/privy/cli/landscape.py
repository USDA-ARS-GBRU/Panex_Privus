"""``privy landscape`` — windowed VCF landscape and local background maps."""

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

log = logging.getLogger("privy.cli.landscape")

LANDSCAPE_HELP = (
    "Create target/off-target-aware VCF sliding-window landscapes, local "
    "background maps, and sample-similarity summaries.\n\n"
    "[bold]Cohort syntax:[/bold] "
    "--targets SAMPLE [SAMPLE ...] --off-targets SAMPLE [SAMPLE ...]"
)

LANDSCAPE_CONTEXT_SETTINGS = {
    "allow_extra_args": True,
    "ignore_unknown_options": True,
}

app = typer.Typer(
    name="landscape",
    help=LANDSCAPE_HELP,
    context_settings=LANDSCAPE_CONTEXT_SETTINGS,
    rich_markup_mode="rich",
    no_args_is_help=True,
)


@app.callback(
    invoke_without_command=True,
    context_settings=LANDSCAPE_CONTEXT_SETTINGS,
)
def landscape(
    ctx: typer.Context,
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
    window_records: int = typer.Option(
        200, "--window-records", metavar="INTEGER", min=1,
        help="Number of VCF records per record-based window.",
    ),
    step_records: int = typer.Option(
        50, "--step-records", metavar="INTEGER", min=1,
        help="Number of VCF records to advance between record-based windows.",
    ),
    window_bp: int | None = typer.Option(
        None, "--window-bp", metavar="INTEGER", min=1,
        help="Use base-pair windows of this size instead of record-count windows.",
    ),
    step_bp: int | None = typer.Option(
        None, "--step-bp", metavar="INTEGER", min=1,
        help="Base-pair step size. Defaults to --window-bp when omitted.",
    ),
    pass_only: bool = typer.Option(
        True, "--pass-only/--no-pass-only",
        help="Use only PASS VCF records.",
    ),
    min_qual: float | None = typer.Option(
        None, "--min-qual", metavar="FLOAT",
        help="Minimum VCF QUAL to include.",
    ),
    rare_max_count: int = typer.Option(
        1, "--rare-max-count", metavar="INTEGER", min=0,
        help="ALT allele carrier-count threshold used for rare-burden windows.",
    ),
    rare_max_freq: float = typer.Option(
        0.05, "--rare-max-freq", metavar="FLOAT", min=0.0, max=1.0,
        help="ALT allele carrier-frequency threshold used for rare-burden windows.",
    ),
    min_called_for_freq: int = typer.Option(
        10, "--min-called-for-freq", metavar="INTEGER", min=0,
        help="Minimum called samples at a variant before genotype-frequency values are used.",
    ),
    min_freq_values: int = typer.Option(
        10, "--min-freq-values", metavar="INTEGER", min=0,
        help="Minimum per-sample values in a window before median_call_freq is reported.",
    ),
    min_background_similarity: float = typer.Option(
        0.65, "--min-background-similarity", metavar="FLOAT", min=0.0, max=1.0,
        help="Minimum nearest-sample similarity for assigning local background blocks.",
    ),
    min_introgression_similarity: float | None = typer.Option(
        None, "--min-introgression-similarity", metavar="FLOAT", min=0.0, max=1.0,
        help=(
            "Minimum off-target similarity for candidate introgression blocks. "
            "Defaults to --min-background-similarity."
        ),
    ),
    min_introgression_delta: float = typer.Option(
        0.0, "--min-introgression-delta", metavar="FLOAT", min=0.0, max=1.0,
        help=(
            "Minimum similarity advantage of the nearest off-target over the "
            "nearest target sample."
        ),
    ),
    max_introgression_missing_rate: float = typer.Option(
        0.5, "--max-introgression-missing-rate", metavar="FLOAT", min=0.0, max=1.0,
        help="Maximum target-sample missingness allowed in candidate introgression windows.",
    ),
    min_introgression_windows: int = typer.Option(
        1, "--min-introgression-windows", metavar="INTEGER", min=1,
        help="Minimum adjacent windows required to emit a candidate introgression block.",
    ),
    similarity_output: str = typer.Option(
        "summary", "--similarity-output", metavar="TEXT",
        help="Pairwise similarity output mode: full, summary, or none.",
    ),
    vcf_engine: str = typer.Option(
        "auto", "--vcf-engine", metavar="TEXT",
        help="VCF parser engine: auto, pysam, or cyvcf2.",
    ),
    local_pca: bool = typer.Option(
        False, "--local-pca/--no-local-pca",
        help="Write exploratory local PCA coordinates from each window similarity matrix.",
    ),
    write_plots: bool = typer.Option(
        True, "--plots/--no-plots",
        help="Write landscape figures alongside TSV outputs.",
    ),
    plot_format: str = typer.Option(
        "png", "--plot-format", metavar="TEXT",
        help="Plot format: png, svg, or pdf.",
    ),
    outdir: Path | None = typer.Option(
        None, "--outdir", metavar="PATH",
        help="Output directory for landscape tables and plots.",
    ),
) -> None:
    """Create VCF window metrics and local background maps."""
    state = get_state()
    effective_outdir = outdir or state.outdir

    if vcf is None:
        typer.echo("[error] Provide a VCF input with --vcf.", err=True)
        raise typer.Exit(code=1)
    if not vcf.exists():
        typer.echo(f"[error] VCF file not found: {vcf}", err=True)
        raise typer.Exit(code=1)
    if plot_format not in {"png", "svg", "pdf"}:
        typer.echo("[error] --plot-format must be one of: png, svg, pdf.", err=True)
        raise typer.Exit(code=1)
    if similarity_output not in {"full", "summary", "none"}:
        typer.echo(
            "[error] --similarity-output must be one of: full, summary, none.",
            err=True,
        )
        raise typer.Exit(code=1)
    if vcf_engine not in {"auto", "pysam", "cyvcf2"}:
        typer.echo("[error] --vcf-engine must be one of: auto, pysam, cyvcf2.", err=True)
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

    try:
        from privy.backends.landscape import run_landscape_vcf  # noqa: PLC0415

        run_landscape_vcf(
            vcf=vcf,
            targets=target_list,
            off_targets=off_target_list or None,
            ignored_samples=ignored_list,
            outdir=effective_outdir,
            window_records=window_records,
            step_records=step_records,
            window_bp=window_bp,
            step_bp=step_bp,
            pass_only=pass_only,
            min_qual=min_qual,
            rare_max_count=rare_max_count,
            rare_max_freq=rare_max_freq,
            min_called_for_freq=min_called_for_freq,
            min_freq_values=min_freq_values,
            min_background_similarity=min_background_similarity,
            min_introgression_similarity=min_introgression_similarity,
            min_introgression_delta=min_introgression_delta,
            max_introgression_missing_rate=max_introgression_missing_rate,
            min_introgression_windows=min_introgression_windows,
            similarity_output=similarity_output,
            vcf_engine=vcf_engine,
            local_pca=local_pca,
            write_plots=write_plots,
            plot_format=plot_format,
        )
    except (FileNotFoundError, ValueError) as exc:
        typer.echo(f"[error] {exc}", err=True)
        raise typer.Exit(code=1) from exc

    typer.echo(f"privy landscape complete. Outputs in: {effective_outdir}")
