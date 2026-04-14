"""``privy compare`` — cross-evidence reconciliation engine.

Compares loci or regions across VCF, BAM, GFA, and XMFA evidence sources
and classifies each locus as supported, contradicted, source-specific,
uninformative, or missing_data.

This is the intellectual differentiator of Panex Privus: it treats
cross-evidence agreement as a first-class analytical feature rather than
an afterthought.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import List, Optional

import typer

from privy.cli.context import get_state
from privy.core.config import PrivyConfig, load_config

log = logging.getLogger("privy.cli.compare")

app = typer.Typer(
    name="compare",
    help=(
        "Compare loci or regions across evidence sources.\n\n"
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
    # --------------------------------------------------------------- inputs
    hits: Optional[Path] = typer.Option(
        None, "--hits", metavar="PATH",
        help="hits.tsv from privy scan.",
    ),
    regions: Optional[Path] = typer.Option(
        None, "--regions", metavar="PATH",
        help="regions.tsv from privy scan.",
    ),
    vcf: Optional[Path] = typer.Option(
        None, "--vcf", metavar="PATH",
        help="Indexed VCF for comparison.",
    ),
    bam: Optional[List[Path]] = typer.Option(
        None, "--bam", metavar="PATH",
        help="BAM files for comparison. Repeat flag for multiple files.",
    ),
    bam_manifest: Optional[Path] = typer.Option(
        None, "--bam-manifest", metavar="PATH",
        help="BAM manifest TSV.",
    ),
    gfa: Optional[Path] = typer.Option(
        None, "--gfa", metavar="PATH",
        help="GFA graph file.",
    ),
    xmfa: Optional[Path] = typer.Option(
        None, "--xmfa", metavar="PATH",
        help="XMFA alignment file.",
    ),
    source_a: Optional[Path] = typer.Option(
        None, "--a", metavar="PATH",
        help="First comparison input (for scan_vs_scan mode).",
    ),
    source_b: Optional[Path] = typer.Option(
        None, "--b", metavar="PATH",
        help="Second comparison input (for scan_vs_scan mode).",
    ),
    # ----------------------------------------------------------- compare options
    mode: str = typer.Option(
        "multi_evidence", "--mode", metavar="TEXT",
        help=(
            "Comparison mode: vcf_vs_bam | vcf_vs_gfa | vcf_vs_xmfa | "
            "scan_vs_scan | multi_evidence."
        ),
    ),
    overlap_mode: str = typer.Option(
        "reciprocal", "--overlap-mode", metavar="TEXT",
        help="Overlap mode: any | reciprocal | contained.",
    ),
    min_reciprocal_overlap: Optional[float] = typer.Option(
        None, "--min-reciprocal-overlap", metavar="FLOAT", min=0.0, max=1.0,
        help="Minimum reciprocal overlap for interval matching.",
    ),
    breakpoint_tolerance_bp: Optional[int] = typer.Option(
        None, "--breakpoint-tolerance-bp", metavar="INTEGER", min=0,
        help="Tolerance (bp) for breakpoint-aware comparisons.",
    ),
    require_state_compatibility: bool = typer.Option(
        False, "--require-state-compatibility/--no-require-state-compatibility",
        help="Require allele/state compatibility in addition to coordinate overlap.",
    ),
    # --------------------------------------------------------------- outputs
    write_compare_tsv: bool = typer.Option(
        True, "--write-compare-tsv/--no-write-compare-tsv",
        help="Write compare.tsv.",
    ),
    write_summary_tsv: bool = typer.Option(
        True, "--write-summary-tsv/--no-write-summary-tsv",
        help="Write compare_summary.tsv.",
    ),
    write_json: bool = typer.Option(
        True, "--write-json/--no-write-json",
        help="Write compare.json.",
    ),
    outdir: Optional[Path] = typer.Option(
        None, "--outdir", metavar="PATH",
        help="Output directory. Overrides global --outdir.",
    ),
) -> None:
    """Compare loci or regions across evidence sources.

    Comparison is based on interval overlap plus evidence compatibility.
    Different evidence layers may support, contradict, or fail to inform a locus.
    """
    from privy.compare.engine import run_compare  # noqa: PLC0415

    state = get_state()
    effective_outdir = outdir or state.outdir

    cfg: PrivyConfig
    if state.config_path is not None:
        cfg = load_config(state.config_path)
    else:
        from privy.core.config import default_config
        cfg = default_config()

    if min_reciprocal_overlap is not None:
        cfg = cfg.model_copy(
            update={"compare": cfg.compare.model_copy(
                update={"min_reciprocal_overlap": min_reciprocal_overlap}
            )}
        )
    if breakpoint_tolerance_bp is not None:
        cfg = cfg.model_copy(
            update={"compare": cfg.compare.model_copy(
                update={"breakpoint_tolerance_bp": breakpoint_tolerance_bp}
            )}
        )

    if hits is None and source_a is None:
        typer.echo(
            "[error] Provide --hits (from privy scan) or --a/--b for scan_vs_scan mode.",
            err=True,
        )
        raise typer.Exit(code=1)

    effective_outdir.mkdir(parents=True, exist_ok=True)
    log.info("Starting compare | mode=%s", mode)

    run_compare(
        hits=hits,
        regions=regions,
        vcf=vcf,
        bam=bam,
        bam_manifest=bam_manifest,
        gfa=gfa,
        xmfa=xmfa,
        source_a=source_a,
        source_b=source_b,
        cfg=cfg,
        mode=mode,
        outdir=effective_outdir,
        write_compare_tsv=write_compare_tsv,
        write_summary_tsv=write_summary_tsv,
        write_json=write_json,
    )
