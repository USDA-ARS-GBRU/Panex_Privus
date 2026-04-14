"""``privy scan`` — primary discovery engine.

Discovers target-private alleles and candidate private regions from
VCF-first workflows with optional BAM, GFA, and XMFA support layers.

This module wires all CLI options and resolves them against the YAML config
before dispatching to :mod:`privy.backends.vcf_scan`.
"""

from __future__ import annotations

import logging
from pathlib import Path

import typer

from privy.cli.context import get_state
from privy.core.cohort import CohortDefinition
from privy.core.config import PrivyConfig, load_config

log = logging.getLogger("privy.cli.scan")

app = typer.Typer(
    name="scan",
    help=(
        "Discover target-private alleles and candidate private regions.\n\n"
        "VCF is the primary discovery backend. BAM, GFA, and XMFA are support layers.\n"
        "Missingness is reported via [italic]strictness_class[/italic] in all outputs.\n\n"
        "[bold]Outputs:[/bold] hits.tsv, regions.tsv, evidence.tsv, "
        "sample_support.tsv, qc.tsv, run.json"
    ),
    rich_markup_mode="rich",
    no_args_is_help=True,
)


@app.callback(invoke_without_command=True)
def scan(
    # ------------------------------------------------------------------ inputs
    vcf: Path | None = typer.Option(
        None, "--vcf", metavar="PATH",
        help="Indexed multisample VCF (.vcf.gz + .tbi); primary v1 backend.",
    ),
    xmfa: Path | None = typer.Option(
        None, "--xmfa", metavar="PATH",
        help="XMFA alignment file; optional secondary or alternate input.",
    ),
    gfa: Path | None = typer.Option(
        None, "--gfa", metavar="PATH",
        help="GFA graph file; optional graph-context support layer.",
    ),
    bam: list[Path] | None = typer.Option(
        None, "--bam", metavar="PATH",
        help="One or more BAM files. Repeat flag for multiple files.",
    ),
    bam_manifest: Path | None = typer.Option(
        None, "--bam-manifest", metavar="PATH",
        help="TSV manifest mapping BAM files to sample names/groups.",
    ),
    # --------------------------------------------------------------- cohort
    targets: list[str] | None = typer.Option(
        None, "--targets", metavar="TEXT",
        help="Target sample names. Repeat flag for multiple samples.",
    ),
    off_targets: list[str] | None = typer.Option(
        None, "--off-targets", metavar="TEXT",
        help="Off-target sample names. Repeat flag for multiple samples.",
    ),
    ignore_samples: list[str] | None = typer.Option(
        None, "--ignore-samples", metavar="TEXT",
        help="Samples to ignore during discovery.",
    ),
    cohort_file: Path | None = typer.Option(
        None, "--cohort-file", metavar="PATH",
        help="Optional cohort definition file (TSV or YAML).",
    ),
    # ------------------------------------------------------------- discovery
    mode: str = typer.Option(
        "private_allele", "--mode", metavar="TEXT",
        help="Discovery mode: private_allele | private_genotype | private_sv_state.",
    ),
    min_target_support: float | None = typer.Option(
        None, "--min-target-support", metavar="FLOAT", min=0.0, max=1.0,
        help="Minimum fraction of target samples supporting allele.",
    ),
    max_off_target_support: float | None = typer.Option(
        None, "--max-off-target-support", metavar="FLOAT", min=0.0, max=1.0,
        help="Maximum fraction of off-target samples supporting allele.",
    ),
    allow_multiallelic: bool = typer.Option(
        True, "--allow-multiallelic/--no-allow-multiallelic",
        help="Whether to evaluate multiallelic records.",
    ),
    pass_only: bool = typer.Option(
        True, "--pass-only/--no-pass-only",
        help="Require VCF FILTER=PASS.",
    ),
    min_qual: float | None = typer.Option(
        None, "--min-qual", metavar="FLOAT",
        help="Minimum VCF QUAL score.",
    ),
    region: str | None = typer.Option(
        None, "--region", metavar="TEXT",
        help="Restrict scan to region: contig:start-end.",
    ),
    contig: str | None = typer.Option(
        None, "--contig", metavar="TEXT",
        help="Restrict scan to a single contig.",
    ),
    chunk_size: int | None = typer.Option(
        None, "--chunk-size", metavar="INTEGER", min=1000,
        help="Chunk size (bp) for streaming large contigs.",
    ),
    merge_distance: int | None = typer.Option(
        None, "--merge-distance", metavar="INTEGER", min=0,
        help="Merge nearby passing loci into candidate regions within this bp distance.",
    ),
    same_variant_class_only: bool = typer.Option(
        False, "--same-variant-class-only/--no-same-variant-class-only",
        help="Only merge loci of the same variant class.",
    ),
    # ----------------------------------------------------------- strictness
    strictness_report: bool = typer.Option(
        True, "--strictness-report/--no-strictness-report",
        help="Report strictness classes explicitly in outputs.",
    ),
    relaxed_target_missing: float | None = typer.Option(
        None, "--relaxed-target-missing", metavar="FLOAT", min=0.0, max=1.0,
        help="Tolerated target missingness fraction for relaxed_threshold class.",
    ),
    relaxed_offtarget_missing: float | None = typer.Option(
        None, "--relaxed-offtarget-missing", metavar="FLOAT", min=0.0, max=1.0,
        help="Tolerated off-target missingness fraction for relaxed_threshold class.",
    ),
    # -------------------------------------------------------- BAM support
    bam_min_depth: int | None = typer.Option(
        None, "--bam-min-depth", metavar="INTEGER", min=0,
        help="Minimum depth for BAM evidence evaluation.",
    ),
    bam_min_alt_count: int | None = typer.Option(
        None, "--bam-min-alt-count", metavar="INTEGER", min=0,
        help="Minimum alternate-supporting read count.",
    ),
    bam_min_alt_fraction: float | None = typer.Option(
        None, "--bam-min-alt-fraction", metavar="FLOAT", min=0.0, max=1.0,
        help="Minimum alternate allele fraction.",
    ),
    summarize_softclips: bool = typer.Option(
        False, "--summarize-softclips/--no-summarize-softclips",
        help="Summarize soft-clipped reads near candidate loci.",
    ),
    summarize_splitreads: bool = typer.Option(
        False, "--summarize-splitreads/--no-summarize-splitreads",
        help="Summarize split-read support near candidate loci.",
    ),
    # --------------------------------------------------------- GFA support
    junction_window_bp: int | None = typer.Option(
        None, "--junction-window-bp", metavar="INTEGER", min=0,
        help="Window around locus for branch-junction annotation.",
    ),
    report_path_membership: bool = typer.Option(
        True, "--report-path-membership/--no-report-path-membership",
        help="Report GFA path membership where available.",
    ),
    report_graph_complexity: bool = typer.Option(
        True, "--report-graph-complexity/--no-report-graph-complexity",
        help="Summarize local graph complexity near candidate loci.",
    ),
    # -------------------------------------------------------- XMFA support
    gap_aware: bool = typer.Option(
        True, "--gap-aware/--no-gap-aware",
        help="Use gap-aware alignment corroboration.",
    ),
    xmfa_window_bp: int | None = typer.Option(
        None, "--xmfa-window-bp", metavar="INTEGER", min=0,
        help="Window (bp) for local XMFA corroboration.",
    ),
    # --------------------------------------------------------------- scoring
    discovery_weight: float | None = typer.Option(
        None, "--discovery-weight", metavar="FLOAT", min=0.0,
        help="Weight for discovery score.",
    ),
    support_weight: float | None = typer.Option(
        None, "--support-weight", metavar="FLOAT", min=0.0,
        help="Weight for support score.",
    ),
    penalty_weight: float | None = typer.Option(
        None, "--penalty-weight", metavar="FLOAT", min=0.0,
        help="Weight for penalty score.",
    ),
    # --------------------------------------------------------------- outputs
    write_hits: bool = typer.Option(True, "--write-hits/--no-write-hits",
                                    help="Write hits.tsv."),
    write_regions: bool = typer.Option(True, "--write-regions/--no-write-regions",
                                       help="Write regions.tsv."),
    write_evidence: bool = typer.Option(True, "--write-evidence/--no-write-evidence",
                                        help="Write evidence.tsv."),
    write_sample_support: bool = typer.Option(
        True, "--write-sample-support/--no-write-sample-support",
        help="Write sample_support.tsv.",
    ),
    write_qc: bool = typer.Option(True, "--write-qc/--no-write-qc",
                                  help="Write qc.tsv."),
    write_run_json: bool = typer.Option(True, "--write-run-json/--no-write-run-json",
                                        help="Write run.json."),
    outdir: Path | None = typer.Option(
        None, "--outdir", metavar="PATH",
        help="Output directory. Overrides global --outdir.",
    ),
) -> None:
    """Discover target-private alleles and candidate private regions.

    VCF and GFA are primary discovery backends.  Missingness is reported
    via strictness_class and never silently folded into pass/fail logic.
    """
    state = get_state()
    effective_outdir = outdir or state.outdir

    # ------------------------------------------------------------------ config
    cfg: PrivyConfig
    if state.config_path is not None:
        log.info("Loading config from %s", state.config_path)
        cfg = load_config(state.config_path)
    else:
        from privy.core.config import default_config
        cfg = default_config()

    # Apply CLI overrides to config
    if state.project_name:
        cfg = cfg.model_copy(update={"project_name": state.project_name})
    if min_target_support is not None:
        cfg = cfg.model_copy(
            update={"scan": cfg.scan.model_copy(update={"min_target_support": min_target_support})}
        )
    if max_off_target_support is not None:
        cfg = cfg.model_copy(
            update={"scan": cfg.scan.model_copy(
                update={"max_off_target_support": max_off_target_support}
            )}
        )
    if merge_distance is not None:
        cfg = cfg.model_copy(
            update={"scan": cfg.scan.model_copy(update={"merge_distance": merge_distance})}
        )
    if chunk_size is not None:
        cfg = cfg.model_copy(
            update={"scan": cfg.scan.model_copy(update={"chunk_size": chunk_size})}
        )

    # ----------------------------------------------------------------- cohort
    # Cohort from CLI takes precedence over YAML cohort
    effective_targets: list[str] = targets or list(cfg.cohorts.targets)
    effective_off_targets: list[str] = off_targets or list(cfg.cohorts.off_targets)
    effective_ignored: list[str] = ignore_samples or list(cfg.cohorts.ignored_samples)

    if not effective_targets or not effective_off_targets:
        typer.echo(
            "[error] Cohort is incomplete. Provide --targets and --off-targets "
            "(or define cohorts in your YAML config).",
            err=True,
        )
        raise typer.Exit(code=1)

    try:
        cohort = CohortDefinition.from_lists(
            targets=effective_targets,
            off_targets=effective_off_targets,
            ignored_samples=effective_ignored,
        )
    except ValueError as exc:
        typer.echo(f"[error] Cohort definition error: {exc}", err=True)
        raise typer.Exit(code=1) from exc

    # ------------------------------------------------------------------ input
    if vcf is None and gfa is None and xmfa is None:
        typer.echo(
            "[error] At least one primary input is required: --vcf, --gfa, or --xmfa",
            err=True,
        )
        raise typer.Exit(code=1)

    if vcf is not None and not vcf.exists():
        typer.echo(f"[error] VCF file not found: {vcf}", err=True)
        raise typer.Exit(code=1)

    if gfa is not None and not gfa.exists():
        typer.echo(f"[error] GFA file not found: {gfa}", err=True)
        raise typer.Exit(code=1)

    # --------------------------------------------------------------- outdir
    effective_outdir.mkdir(parents=True, exist_ok=True)
    log.info("Output directory: %s", effective_outdir)

    # ------------------------------------------------------------------- run
    log.info(
        "Starting scan | targets=%d off_targets=%d mode=%s",
        cohort.n_targets, cohort.n_off_targets, mode,
    )

    try:
        if vcf is not None:
            # VCF-first: use the VCF backend; GFA/XMFA/BAM are optional layers
            from privy.backends.vcf_scan import run_vcf_scan  # noqa: PLC0415

            run_vcf_scan(
                vcf=vcf,
                cohort=cohort,
                cfg=cfg,
                outdir=effective_outdir,
                mode=mode,
                bam=bam,
                bam_manifest=bam_manifest,
                gfa=gfa,
                xmfa=xmfa,
                region=region,
                contig=contig,
                write_hits=write_hits,
                write_regions=write_regions,
                write_evidence=write_evidence,
                write_sample_support=write_sample_support,
                write_qc=write_qc,
                write_run_json=write_run_json,
                threads=state.threads,
            )
        elif gfa is not None:
            # GFA-only primary scan
            from privy.backends.gfa_scan import run_gfa_scan  # noqa: PLC0415

            run_gfa_scan(
                gfa=gfa,
                cohort=cohort,
                cfg=cfg,
                outdir=effective_outdir,
                mode=mode,
                region=region,
                contig=contig,
                write_hits=write_hits,
                write_regions=write_regions,
                write_evidence=write_evidence,
                write_sample_support=write_sample_support,
                write_qc=write_qc,
                write_run_json=write_run_json,
                threads=state.threads,
            )
        else:
            raise NotImplementedError("XMFA-only scan is not yet implemented.")
    except (FileNotFoundError, ValueError) as exc:
        typer.echo(f"[error] {exc}", err=True)
        raise typer.Exit(code=1) from exc
    except NotImplementedError as exc:
        typer.echo(f"[error] Not implemented: {exc}", err=True)
        raise typer.Exit(code=2) from exc
