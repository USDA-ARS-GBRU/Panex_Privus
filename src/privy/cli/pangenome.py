"""``privy pangenome`` — whole and sub-pangenome summaries."""

from __future__ import annotations

import logging
from pathlib import Path

import typer

from privy.cli.context import get_state

log = logging.getLogger("privy.cli.pangenome")

app = typer.Typer(
    name="pangenome",
    help=(
        "Analyze full, target, and off-target pangenomes from shared feature "
        "matrices. GFA segment-level analysis is available first."
    ),
    rich_markup_mode="rich",
    no_args_is_help=True,
)


@app.callback(invoke_without_command=True)
def pangenome(
    gfa: Path | None = typer.Option(
        None, "--gfa", metavar="PATH",
        help="GFA graph file to analyze.",
    ),
    vcf: Path | None = typer.Option(
        None, "--vcf", metavar="PATH",
        help="Reserved for the upcoming VCF pangenome adapter.",
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
    permutations: int = typer.Option(
        100, "--permutations", metavar="INTEGER", min=1,
        help="Number of rarefaction permutations.",
    ),
    seed: int = typer.Option(
        42, "--seed", metavar="INTEGER",
        help="Random seed used for deterministic permutations.",
    ),
    write_plots: bool = typer.Option(
        True, "--plots/--no-plots",
        help="Write first-pass pangenome plots alongside TSV outputs.",
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

    target_list = _merge_sample_flags(targets, targets_file)
    off_target_list = _merge_sample_flags(off_targets, off_targets_file)
    ignored_list = _merge_sample_flags(ignore_samples, ignore_samples_file)

    if not target_list:
        typer.echo(
            "[error] Provide target samples with --targets and/or --targets-file.",
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
