"""``privy interactive`` — self-contained HTML dashboards."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import cast

import typer

from privy.cli.context import get_state
from privy.interactive.focus import RECOMMENDED_FOCUS_BP, run_focus_dashboards
from privy.interactive.genotypes import VariantFilter
from privy.interactive.models import FocusRegion, parse_focus_region

log = logging.getLogger("privy.cli.interactive")

app = typer.Typer(
    name="interactive",
    help=(
        "Build self-contained interactive HTML dashboards.\n\n"
        "The first supported mode is one static dashboard per --focus region. "
        "Start with regions around 4 Mbp or smaller for responsive review."
    ),
    rich_markup_mode="rich",
    no_args_is_help=True,
)


@app.callback(invoke_without_command=True)
def interactive(
    focus: list[str] | None = typer.Option(
        None,
        "--focus",
        metavar="REGION",
        help=(
            "Focus region to render as CONTIG:START-END. Repeat for multiple "
            "regions; one HTML file is written per region."
        ),
    ),
    sites_tsv: Path | None = typer.Option(
        None,
        "--sites-tsv",
        metavar="PATH",
        help=(
            "Precomputed focal genotype sites TSV. When omitted, Privy extracts "
            "one sites TSV per --focus region from --vcf."
        ),
    ),
    vcf: Path | None = typer.Option(
        None,
        "--vcf",
        metavar="PATH",
        help="Multisample VCF/BCF to extract focus-region sites TSVs from.",
    ),
    gff3: Path | None = typer.Option(
        None,
        "--gff3",
        metavar="PATH",
        help="Reference gene GFF3/GFF3.gz used for genes, exons, CDS, introns, and promoters.",
    ),
    samples: tuple[str, str, str] | None = typer.Option(
        None,
        "--samples",
        metavar="OFFTARGET DERIVED DONOR",
        help="Three sample names in the sites TSV: off-target, derived line, donor.",
    ),
    track_gff: list[str] | None = typer.Option(
        None,
        "--track-gff",
        metavar="LABEL=PATH",
        help="Additional GFF3/GFF3.gz feature track to embed. Repeatable.",
    ),
    functional_tsv: Path | None = typer.Option(
        None,
        "--functional-tsv",
        metavar="PATH",
        help="Optional gene functional annotation TSV to join by gene or gene_id.",
    ),
    sample_abbrev: list[str] | None = typer.Option(
        None,
        "--sample-abbrev",
        metavar="ABBR=SAMPLE",
        help="Display abbreviation for sample genotypes. Repeatable.",
    ),
    keyword_group: list[str] | None = typer.Option(
        None,
        "--keyword-group",
        metavar="NAME=term1,term2",
        help=(
            "Phenotype-oriented candidate group selected from feature and "
            "functional text. Repeatable."
        ),
    ),
    promoter_bp: int = typer.Option(
        2000,
        "--promoter-bp",
        metavar="INTEGER",
        min=0,
        help="Strand-aware upstream promoter length to display and rank.",
    ),
    sv_size_threshold: int = typer.Option(
        50,
        "--sv-size-threshold",
        metavar="INTEGER",
        min=1,
        help="Minimum allele length for the size-based SV-like display layer.",
    ),
    candidate_limit: int = typer.Option(
        60,
        "--candidate-limit",
        metavar="INTEGER",
        min=1,
        help="Maximum variant-supported feature rows to rank per focus region.",
    ),
    pass_only: bool = typer.Option(
        True,
        "--pass-only/--no-pass-only",
        help="When extracting from --vcf, keep only records with FILTER=PASS.",
    ),
    require_all_called: bool = typer.Option(
        True,
        "--require-all-called/--allow-missing-genotypes",
        help="When extracting from --vcf, require all three focus samples to have called GT.",
    ),
    variant_type: str = typer.Option(
        "all",
        "--variant-type",
        metavar="TEXT",
        help="When extracting from --vcf, keep variant type: all, snp, indel, or sv.",
    ),
    biallelic_only: bool = typer.Option(
        False,
        "--biallelic-only/--allow-multiallelic",
        help="When extracting from --vcf, keep only records with one ALT allele.",
    ),
    title: str | None = typer.Option(
        None,
        "--title",
        metavar="TEXT",
        help="Optional dashboard title.",
    ),
    subtitle: str | None = typer.Option(
        None,
        "--subtitle",
        metavar="TEXT",
        help="Optional dashboard subtitle.",
    ),
    outdir: Path | None = typer.Option(
        None,
        "--outdir",
        metavar="PATH",
        help="Output directory. Overrides global --outdir.",
    ),
) -> None:
    """Build one self-contained interactive HTML file per focus region."""
    state = get_state()
    effective_outdir = outdir or state.outdir

    try:
        focus_regions = _parse_focuses(focus)
        parsed_tracks = _parse_label_paths(track_gff or [], "--track-gff")
        parsed_abbrev = _parse_sample_abbrev(sample_abbrev or [])
        parsed_keyword_groups = _parse_keyword_groups(keyword_group or [])
    except ValueError as exc:
        typer.echo(f"[error] {exc}", err=True)
        raise typer.Exit(code=2) from exc

    if sites_tsv is None and vcf is None:
        typer.echo("[error] Provide either --sites-tsv or --vcf.", err=True)
        raise typer.Exit(code=1)
    if variant_type not in {"all", "snp", "indel", "sv"}:
        typer.echo("[error] --variant-type must be one of: all, snp, indel, sv.", err=True)
        raise typer.Exit(code=1)
    if gff3 is None:
        typer.echo("[error] Provide --gff3 for gene and feature tracks.", err=True)
        raise typer.Exit(code=1)
    if samples is None:
        typer.echo(
            "[error] Provide --samples OFFTARGET DERIVED DONOR matching the sites TSV.",
            err=True,
        )
        raise typer.Exit(code=1)
    if vcf is not None and not state.quiet:
        if sites_tsv is not None:
            typer.echo("[warn] --sites-tsv supplied; --vcf extraction will be skipped.")
    for region in focus_regions:
        if region.length > RECOMMENDED_FOCUS_BP and not state.quiet:
            typer.echo(
                "[warn] Focus region "
                f"{region.label} is {region.length / 1_000_000:.2f} Mbp. "
                "For novice use, start with regions around 4 Mbp or smaller.",
            )

    try:
        generated = run_focus_dashboards(
            focuses=focus_regions,
            sites_tsv=sites_tsv,
            vcf=vcf,
            gff3=gff3,
            samples=samples,
            outdir=effective_outdir,
            title=title,
            subtitle=subtitle,
            functional_tsv=functional_tsv,
            track_gff=parsed_tracks,
            sample_abbrev=parsed_abbrev,
            promoter_bp=promoter_bp,
            sv_size_threshold=sv_size_threshold,
            candidate_limit=candidate_limit,
            pass_only=pass_only,
            require_all_called=require_all_called,
            variant_filter=cast(VariantFilter, variant_type),
            biallelic_only=biallelic_only,
            keyword_groups=parsed_keyword_groups,
        )
    except (FileNotFoundError, ValueError) as exc:
        typer.echo(f"[error] {exc}", err=True)
        raise typer.Exit(code=1) from exc

    if not state.quiet:
        for path in generated:
            typer.echo(f"  {path}")


def _parse_focuses(values: list[str] | None) -> list[FocusRegion]:
    if not values:
        raise ValueError("At least one --focus REGION is required.")
    return [parse_focus_region(value) for value in values]


def _parse_label_paths(values: list[str], option_name: str) -> list[tuple[str, Path]]:
    parsed: list[tuple[str, Path]] = []
    for value in values:
        if "=" not in value:
            raise ValueError(f"{option_name} values must be LABEL=PATH.")
        label, raw_path = value.split("=", 1)
        label = label.strip()
        if not label:
            raise ValueError(f"{option_name} label cannot be empty.")
        parsed.append((label, Path(raw_path)))
    return parsed


def _parse_sample_abbrev(values: list[str]) -> dict[str, str]:
    parsed: dict[str, str] = {}
    for value in values:
        if "=" not in value:
            raise ValueError("--sample-abbrev values must be ABBR=SAMPLE.")
        abbr, sample = value.split("=", 1)
        if not abbr or not sample:
            raise ValueError("--sample-abbrev values must be ABBR=SAMPLE.")
        parsed[sample] = abbr
    return parsed


def _parse_keyword_groups(values: list[str]) -> list[tuple[str, list[str]]]:
    parsed: list[tuple[str, list[str]]] = []
    seen: set[str] = set()
    for value in values:
        if "=" not in value:
            raise ValueError("--keyword-group values must be NAME=term1,term2.")
        name, raw_terms = value.split("=", 1)
        name = name.strip()
        terms = [term.strip().lower() for term in raw_terms.split(",") if term.strip()]
        if not name or not terms:
            raise ValueError("--keyword-group values must be NAME=term1,term2.")
        key = name.lower()
        if key in seen:
            raise ValueError(f"Duplicate --keyword-group name: {name}")
        seen.add(key)
        parsed.append((name, terms))
    return parsed
