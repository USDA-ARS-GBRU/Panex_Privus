"""``privy landscape`` — windowed VCF landscape and local background maps."""

from __future__ import annotations

import logging
from pathlib import Path

import typer

from privy.cli.context import get_state

log = logging.getLogger("privy.cli.landscape")

app = typer.Typer(
    name="landscape",
    help=(
        "Create target/off-target-aware VCF sliding-window landscapes, local "
        "background maps, and sample-similarity summaries."
    ),
    rich_markup_mode="rich",
    no_args_is_help=True,
)


@app.callback(invoke_without_command=True)
def landscape(
    vcf: Path | None = typer.Option(
        None, "--vcf", metavar="PATH",
        help="Multisample VCF/BCF file to analyze.",
    ),
    targets: list[str] | None = typer.Option(
        None, "--targets", metavar="TEXT",
        help="Target sample name. Repeat flag for multiple samples.",
    ),
    targets_file: Path | None = typer.Option(
        None, "--targets-file", metavar="PATH",
        help="Text file with one target sample name per line.",
    ),
    off_targets: list[str] | None = typer.Option(
        None, "--off-targets", metavar="TEXT",
        help="Off-target sample name. Repeat flag for multiple samples.",
    ),
    off_targets_file: Path | None = typer.Option(
        None, "--off-targets-file", metavar="PATH",
        help="Text file with one off-target sample name per line.",
    ),
    ignore_samples: list[str] | None = typer.Option(
        None, "--ignore-samples", metavar="TEXT",
        help="Sample to exclude. Repeat flag for multiple samples.",
    ),
    ignore_samples_file: Path | None = typer.Option(
        None, "--ignore-samples-file", metavar="PATH",
        help="Text file with one sample name to ignore per line.",
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

    target_list = _merge_sample_flags(targets, targets_file)
    off_target_list = _merge_sample_flags(off_targets, off_targets_file)
    ignored_list = _merge_sample_flags(ignore_samples, ignore_samples_file)
    if not target_list:
        typer.echo(
            "[error] Provide target samples with --targets and/or --targets-file.",
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
            write_plots=write_plots,
            plot_format=plot_format,
        )
    except (FileNotFoundError, ValueError) as exc:
        typer.echo(f"[error] {exc}", err=True)
        raise typer.Exit(code=1) from exc

    typer.echo(f"privy landscape complete. Outputs in: {effective_outdir}")


def _merge_sample_flags(values: list[str] | None, path: Path | None) -> list[str]:
    merged: list[str] = []
    for value in values or []:
        if value.strip():
            merged.append(value.strip())
    if path is not None:
        if not path.exists():
            raise FileNotFoundError(f"Sample list file not found: {path}")
        merged.extend(_read_sample_list(path))
    return list(dict.fromkeys(merged))


def _read_sample_list(path: Path) -> list[str]:
    samples: list[str] = []
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        samples.append(line.split()[0])
    return samples
