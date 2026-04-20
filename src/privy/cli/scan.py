"""``privy scan`` — primary discovery engine.

Discovers target-private alleles and candidate private regions from
VCF-first workflows with optional BAM, GFA, and XMFA support layers.

This module wires all CLI options and resolves them against the YAML config
before dispatching to :mod:`privy.backends.vcf_scan`.
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path

import typer
import yaml
from click.core import ParameterSource

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
    ctx: typer.Context,
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
        help="GFA graph file; primary backend when used without --vcf.",
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
    bam_min_mapq: int | None = typer.Option(
        None, "--bam-min-mapq", metavar="INTEGER", min=0,
        help="Minimum mapping quality for BAM reads.",
    ),
    bam_min_baseq: int | None = typer.Option(
        None, "--bam-min-baseq", metavar="INTEGER", min=0,
        help="Minimum base quality for BAM pileup.",
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
    min_segment_length: int | None = typer.Option(
        None, "--min-segment-length", metavar="INTEGER", min=1,
        help="Minimum GFA segment length (bp) to evaluate.",
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

    cfg = _apply_cli_overrides(
        ctx=ctx,
        cfg=cfg,
        project_name=state.project_name,
        min_target_support=min_target_support,
        max_off_target_support=max_off_target_support,
        allow_multiallelic=allow_multiallelic,
        pass_only=pass_only,
        min_qual=min_qual,
        chunk_size=chunk_size,
        merge_distance=merge_distance,
        same_variant_class_only=same_variant_class_only,
        strictness_report=strictness_report,
        relaxed_target_missing=relaxed_target_missing,
        relaxed_offtarget_missing=relaxed_offtarget_missing,
        bam_min_depth=bam_min_depth,
        bam_min_alt_count=bam_min_alt_count,
        bam_min_alt_fraction=bam_min_alt_fraction,
        bam_min_mapq=bam_min_mapq,
        bam_min_baseq=bam_min_baseq,
        summarize_softclips=summarize_softclips,
        summarize_splitreads=summarize_splitreads,
        junction_window_bp=junction_window_bp,
        report_path_membership=report_path_membership,
        report_graph_complexity=report_graph_complexity,
        min_segment_length=min_segment_length,
        gap_aware=gap_aware,
        xmfa_window_bp=xmfa_window_bp,
        discovery_weight=discovery_weight,
        support_weight=support_weight,
        penalty_weight=penalty_weight,
    )

    # ----------------------------------------------------------------- cohort
    cohort_from_file = _load_cohort_file(cohort_file) if cohort_file is not None else None

    # CLI cohort flags take precedence over cohort file, which takes precedence over config
    effective_targets: list[str] = (
        targets
        or (list(cohort_from_file.targets) if cohort_from_file is not None else None)
        or list(cfg.cohorts.targets)
    )
    effective_off_targets: list[str] = (
        off_targets
        or (list(cohort_from_file.off_targets) if cohort_from_file is not None else None)
        or list(cfg.cohorts.off_targets)
    )
    effective_ignored: list[str] = (
        ignore_samples
        or (list(cohort_from_file.ignored_samples) if cohort_from_file is not None else None)
        or list(cfg.cohorts.ignored_samples)
    )

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


def _apply_cli_overrides(
    ctx: typer.Context,
    cfg: PrivyConfig,
    **values: object,
) -> PrivyConfig:
    """Apply only explicitly provided CLI options to the config model."""
    if values.get("project_name"):
        cfg = cfg.model_copy(update={"project_name": values["project_name"]})

    scan_updates = _provided_updates(
        ctx,
        values,
        names=(
            "min_target_support",
            "max_off_target_support",
            "allow_multiallelic",
            "pass_only",
            "min_qual",
            "chunk_size",
            "merge_distance",
            "same_variant_class_only",
            "strictness_report",
            "relaxed_target_missing",
            "relaxed_offtarget_missing",
        ),
    )
    if scan_updates:
        cfg = cfg.model_copy(update={"scan": cfg.scan.model_copy(update=scan_updates)})

    bam_updates = _provided_updates(
        ctx,
        values,
        names=(
            "bam_min_depth",
            "bam_min_alt_count",
            "bam_min_alt_fraction",
            "bam_min_mapq",
            "bam_min_baseq",
            "summarize_softclips",
            "summarize_splitreads",
        ),
        rename={
            "bam_min_depth": "min_depth",
            "bam_min_alt_count": "min_alt_count",
            "bam_min_alt_fraction": "allele_fraction_min",
            "bam_min_mapq": "min_mapq",
            "bam_min_baseq": "min_baseq",
        },
    )
    if bam_updates:
        cfg = cfg.model_copy(update={"bam": cfg.bam.model_copy(update=bam_updates)})

    gfa_updates = _provided_updates(
        ctx,
        values,
        names=(
            "junction_window_bp",
            "report_path_membership",
            "report_graph_complexity",
            "min_segment_length",
        ),
    )
    if gfa_updates:
        cfg = cfg.model_copy(update={"gfa": cfg.gfa.model_copy(update=gfa_updates)})

    xmfa_updates = _provided_updates(
        ctx,
        values,
        names=("gap_aware", "xmfa_window_bp"),
        rename={"xmfa_window_bp": "window_bp"},
    )
    if xmfa_updates:
        cfg = cfg.model_copy(update={"xmfa": cfg.xmfa.model_copy(update=xmfa_updates)})

    scoring_updates = _provided_updates(
        ctx,
        values,
        names=("discovery_weight", "support_weight", "penalty_weight"),
    )
    if scoring_updates:
        cfg = cfg.model_copy(
            update={"scoring": cfg.scoring.model_copy(update=scoring_updates)}
        )

    return cfg


def _provided_updates(
    ctx: typer.Context,
    values: dict[str, object],
    names: tuple[str, ...],
    rename: dict[str, str] | None = None,
) -> dict[str, object]:
    """Return a model-update dict for options explicitly set on the command line."""
    updates: dict[str, object] = {}
    rename = rename or {}
    for name in names:
        if _was_provided(ctx, name):
            updates[rename.get(name, name)] = values[name]
    return updates


def _was_provided(ctx: typer.Context, param_name: str) -> bool:
    """Return True when a parameter came from an explicit CLI flag."""
    return ctx.get_parameter_source(param_name) is ParameterSource.COMMANDLINE


def _load_cohort_file(path: Path) -> CohortDefinition:
    """Load a cohort definition from YAML or TSV."""
    if not path.exists():
        raise FileNotFoundError(f"Cohort file not found: {path}")

    suffix = path.suffix.lower()
    if suffix in {".yaml", ".yml"}:
        return _load_cohort_yaml(path)
    if suffix == ".tsv":
        return _load_cohort_tsv(path)
    raise ValueError(
        f"Unsupported cohort file format: {path.suffix!r}. Use .yaml, .yml, or .tsv."
    )


def _load_cohort_yaml(path: Path) -> CohortDefinition:
    """Load a cohort definition from YAML."""
    with open(path, encoding="utf-8") as fh:
        raw = yaml.safe_load(fh) or {}

    if "cohorts" in raw and isinstance(raw["cohorts"], dict):
        raw = raw["cohorts"]

    return CohortDefinition.from_lists(
        targets=list(raw.get("targets", [])),
        off_targets=list(raw.get("off_targets", [])),
        ignored_samples=list(raw.get("ignored_samples", [])),
    )


def _load_cohort_tsv(path: Path) -> CohortDefinition:
    """Load a cohort definition from a TSV with sample/role columns."""
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Cohort TSV is missing a header row.")

        fieldnames = set(reader.fieldnames)
        sample_field = "sample_id" if "sample_id" in fieldnames else "sample"
        role_field = "cohort_role" if "cohort_role" in fieldnames else "role"

        if sample_field not in fieldnames or role_field not in fieldnames:
            raise ValueError(
                "Cohort TSV must contain sample_id/sample and cohort_role/role columns."
            )

        targets: list[str] = []
        off_targets: list[str] = []
        ignored: list[str] = []

        for row in reader:
            sample = (row.get(sample_field) or "").strip()
            role = (row.get(role_field) or "").strip().lower()
            if not sample:
                continue
            if role in {"target", "targets"}:
                targets.append(sample)
            elif role in {"off_target", "off-target", "offtarget", "background"}:
                off_targets.append(sample)
            elif role in {"ignored", "ignore"}:
                ignored.append(sample)
            else:
                raise ValueError(
                    f"Unsupported cohort role {role!r} for sample {sample!r} in {path}."
                )

    return CohortDefinition.from_lists(
        targets=targets,
        off_targets=off_targets,
        ignored_samples=ignored,
    )
