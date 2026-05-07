"""``privy compare`` — cross-source scan reconciliation.

Compares two ``hits.tsv`` files produced by separate ``privy scan`` runs
(typically one VCF scan and one GFA scan) and classifies each locus pair
as supported, partially_supported, contradicted, source_specific,
uninformative, or missing_data.
"""

from __future__ import annotations

import logging
from pathlib import Path

import typer

from privy.cli.context import get_state
from privy.core.config import PrivyConfig, load_config

log = logging.getLogger("privy.cli.compare")

app = typer.Typer(
    name="compare",
    help=(
        "Reconcile two [bold]privy scan[/bold] result sets by coordinate overlap.\n\n"
        "Provide one [italic]hits.tsv[/italic] from a VCF scan (--hits-a) and one "
        "from a GFA scan (--hits-b), or any two scan runs you want to compare.\n\n"
        "[bold]Match classes:[/bold] "
        "supported | partially_supported | contradicted | "
        "source_specific | uninformative | missing_data\n\n"
        "[bold]Outputs:[/bold] compare.tsv, compare_summary.tsv, compare.json"
    ),
    rich_markup_mode="rich",
    no_args_is_help=True,
)


@app.callback(invoke_without_command=True)
def compare(
    # ---------------------------------------------------------------- inputs
    hits_a: Path = typer.Option(
        ...,
        "--hits-a",
        metavar="PATH",
        help="hits.tsv from the first privy scan run (source A).",
    ),
    hits_b: Path = typer.Option(
        ...,
        "--hits-b",
        metavar="PATH",
        help="hits.tsv from the second privy scan run (source B).",
    ),
    source_label_a: str | None = typer.Option(
        None,
        "--source-a",
        metavar="TEXT",
        help="Display label for source A (default: inferred from locus_id prefix).",
    ),
    source_label_b: str | None = typer.Option(
        None,
        "--source-b",
        metavar="TEXT",
        help="Display label for source B (default: inferred from locus_id prefix).",
    ),
    # --------------------------------------------------- comparison options
    overlap_mode: str | None = typer.Option(
        None,
        "--overlap-mode",
        metavar="TEXT",
        help=("Coordinate matching mode: contained, reciprocal, or any [default: contained]."),
    ),
    min_reciprocal_overlap: float | None = typer.Option(
        None,
        "--min-reciprocal-overlap",
        metavar="FLOAT",
        min=0.0,
        max=1.0,
        help=("Minimum overlap fraction for contained/reciprocal matching [default: 0.5]."),
    ),
    breakpoint_tolerance_bp: int | None = typer.Option(
        None,
        "--breakpoint-tolerance-bp",
        metavar="INT",
        min=0,
        help="Gap tolerance (bp) for near-miss breakpoint matching [default: 200].",
    ),
    require_state_compatibility: bool | None = typer.Option(
        None,
        "--require-state-compatibility/--no-require-state-compatibility",
        help=(
            "Also require strictness-class compatibility (strict_* vs relaxed_threshold) "
            "in addition to coordinate overlap."
        ),
    ),
    normalize_contigs: bool | None = typer.Option(
        None,
        "--normalize-contigs/--no-normalize-contigs",
        help=(
            "Canonicalize minigraph-cactus contigs such as SAMPLE#HAP#CONTIG "
            "to CONTIG before comparing [default: enabled]."
        ),
    ),
    # --------------------------------------------------------------- outputs
    write_compare_tsv: bool = typer.Option(
        True,
        "--write-compare-tsv/--no-write-compare-tsv",
        help="Write compare.tsv.",
    ),
    write_summary_tsv: bool = typer.Option(
        True,
        "--write-summary-tsv/--no-write-summary-tsv",
        help="Write compare_summary.tsv.",
    ),
    write_json: bool = typer.Option(
        True,
        "--write-json/--no-write-json",
        help="Write compare.json.",
    ),
    outdir: Path | None = typer.Option(
        None,
        "--outdir",
        metavar="PATH",
        help="Output directory (overrides global --outdir).",
    ),
) -> None:
    """Reconcile two privy scan result sets by coordinate overlap."""
    from privy.backends.compare import run_compare  # noqa: PLC0415

    state = get_state()
    effective_outdir = outdir or state.outdir

    cfg: PrivyConfig
    if state.config_path is not None:
        cfg = load_config(state.config_path)
    else:
        from privy.core.config import default_config  # noqa: PLC0415

        cfg = default_config()

    compare_updates: dict[str, object] = {}
    if overlap_mode is not None:
        normalized_mode = overlap_mode.lower()
        if normalized_mode not in {"any", "reciprocal", "contained"}:
            typer.echo(
                "[error] --overlap-mode must be one of: any, reciprocal, contained",
                err=True,
            )
            raise typer.Exit(code=1)
        compare_updates["overlap_mode"] = normalized_mode
    if min_reciprocal_overlap is not None:
        compare_updates["min_reciprocal_overlap"] = min_reciprocal_overlap
    if breakpoint_tolerance_bp is not None:
        compare_updates["breakpoint_tolerance_bp"] = breakpoint_tolerance_bp
    if require_state_compatibility is not None:
        compare_updates["require_state_compatibility"] = require_state_compatibility
    if normalize_contigs is not None:
        compare_updates["normalize_contigs"] = normalize_contigs
    if compare_updates:
        cfg = cfg.model_copy(update={"compare": cfg.compare.model_copy(update=compare_updates)})

    if not hits_a.exists():
        typer.echo(f"[error] --hits-a not found: {hits_a}", err=True)
        raise typer.Exit(code=1)
    if not hits_b.exists():
        typer.echo(f"[error] --hits-b not found: {hits_b}", err=True)
        raise typer.Exit(code=1)

    effective_outdir.mkdir(parents=True, exist_ok=True)
    log.info("Starting compare | hits_a=%s hits_b=%s", hits_a, hits_b)

    run_compare(
        hits_a=hits_a,
        hits_b=hits_b,
        outdir=effective_outdir,
        cfg=cfg,
        source_label_a=source_label_a,
        source_label_b=source_label_b,
        write_compare_tsv=write_compare_tsv,
        write_summary_tsv=write_summary_tsv,
        write_json=write_json,
    )
