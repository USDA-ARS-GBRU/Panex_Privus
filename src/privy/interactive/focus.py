"""Build focus-region interactive HTML dashboards."""

from __future__ import annotations

import csv
import gzip
import json
import logging
import re
from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, TypeAlias
from urllib.parse import unquote

from privy.interactive.genotypes import (
    ExtractionSummary,
    VariantFilter,
    extract_focus_sites_from_vcf,
)
from privy.interactive.models import FocusRegion
from privy.interactive.render import render_focus_html, render_index_html

log = logging.getLogger("privy.interactive.focus")

RECOMMENDED_FOCUS_BP = 4_000_000
GeneModelReadResult: TypeAlias = tuple[
    list["GeneModel"],
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[dict[str, Any]],
]


@dataclass
class Transcript:
    transcript_id: str
    gene_id: str
    start: int
    end: int
    strand: str
    longest: bool


@dataclass
class GeneModel:
    gene_id: str
    gene: str
    start: int
    end: int
    strand: str
    transcript_id: str = ""
    exons: list[tuple[int, int]] = field(default_factory=list)
    cds: list[tuple[int, int]] = field(default_factory=list)
    introns: list[tuple[int, int]] = field(default_factory=list)
    promoter: tuple[int, int] | None = None


@dataclass(frozen=True)
class FocusOutput:
    region: FocusRegion
    html: Path
    features_tsv: Path
    metadata_json: Path
    sites_tsv: Path | None
    extraction_summary: ExtractionSummary | None = None


def run_focus_dashboards(
    *,
    focuses: list[FocusRegion],
    sites_tsv: Path | None = None,
    vcf: Path | None = None,
    gff3: Path,
    samples: tuple[str, str, str],
    outdir: Path,
    title: str | None = None,
    subtitle: str | None = None,
    functional_tsv: Path | None = None,
    track_gff: list[tuple[str, Path]] | None = None,
    sample_abbrev: dict[str, str] | None = None,
    promoter_bp: int = 2000,
    sv_size_threshold: int = 50,
    candidate_limit: int = 60,
    pass_only: bool = True,
    require_all_called: bool = True,
    variant_filter: VariantFilter = "all",
    biallelic_only: bool = False,
    keyword_groups: list[tuple[str, list[str]]] | None = None,
) -> list[Path]:
    """Write one static focus HTML per region and an index for multi-focus runs."""
    if not focuses:
        raise ValueError("At least one --focus region is required.")
    if sites_tsv is None and vcf is None:
        raise ValueError("Provide either --sites-tsv or --vcf.")
    if sites_tsv is not None and not sites_tsv.exists():
        raise FileNotFoundError(f"--sites-tsv not found: {sites_tsv}")
    if vcf is not None and not vcf.exists():
        raise FileNotFoundError(f"--vcf not found: {vcf}")
    if not gff3.exists():
        raise FileNotFoundError(f"--gff3 not found: {gff3}")
    if functional_tsv is not None and not functional_tsv.exists():
        raise FileNotFoundError(f"--functional-tsv not found: {functional_tsv}")
    for label, path in track_gff or []:
        if not path.exists():
            raise FileNotFoundError(f"--track-gff {label} path not found: {path}")

    outdir.mkdir(parents=True, exist_ok=True)
    annotations = _read_function_annotations(functional_tsv)
    outputs: list[FocusOutput] = []
    for focus in focuses:
        if focus.length > RECOMMENDED_FOCUS_BP:
            log.warning(
                "Focus region %s is %.2f Mbp; start with <=4 Mbp for responsive "
                "novice-friendly dashboards unless variant density is known to be low.",
                focus.label,
                focus.length / 1_000_000,
            )
        focus_sites_tsv = sites_tsv
        extraction_summary = None
        if focus_sites_tsv is None:
            focus_sites_tsv = outdir / f"{focus.slug}.sites.tsv"
            if vcf is None:
                raise ValueError("Internal error: --vcf is required to extract sites.")
            extraction_summary = extract_focus_sites_from_vcf(
                vcf=vcf,
                gff3=gff3,
                focus=focus,
                samples=samples,
                out_tsv=focus_sites_tsv,
                pass_only=pass_only,
                require_all_called=require_all_called,
                variant_filter=variant_filter,
                biallelic_only=biallelic_only,
            )
        output = _write_one_focus(
            focus=focus,
            sites_tsv=focus_sites_tsv,
            gff3=gff3,
            samples=samples,
            outdir=outdir,
            title=title,
            subtitle=subtitle,
            annotations=annotations,
            functional_tsv=functional_tsv,
            track_gff=track_gff or [],
            sample_abbrev=sample_abbrev or {},
            promoter_bp=promoter_bp,
            sv_size_threshold=sv_size_threshold,
            candidate_limit=candidate_limit,
            extraction_summary=extraction_summary,
            extraction_filters={
                "pass_only": pass_only,
                "require_all_called": require_all_called,
                "variant_filter": variant_filter,
                "biallelic_only": biallelic_only,
            },
            keyword_groups=keyword_groups or [],
        )
        outputs.append(output)

    paths = [output.html for output in outputs]
    if len(outputs) > 1:
        index_path = outdir / "index.html"
        index_path.write_text(
            render_index_html(
                outputs=outputs,
                title=title or "Privy Interactive Focus Regions",
                subtitle=(
                    subtitle
                    or "One self-contained dashboard was written for each --focus region."
                ),
            ),
            encoding="utf-8",
        )
        paths.insert(0, index_path)
    _write_run_metadata(
        outdir / "interactive.json",
        outputs,
        sites_tsv,
        vcf,
        gff3,
        track_gff or [],
    )
    return paths


def _write_one_focus(
    *,
    focus: FocusRegion,
    sites_tsv: Path,
    gff3: Path,
    samples: tuple[str, str, str],
    outdir: Path,
    title: str | None,
    subtitle: str | None,
    annotations: dict[str, dict[str, str]],
    functional_tsv: Path | None,
    track_gff: list[tuple[str, Path]],
    sample_abbrev: dict[str, str],
    promoter_bp: int,
    sv_size_threshold: int,
    candidate_limit: int,
    extraction_summary: ExtractionSummary | None,
    extraction_filters: dict[str, object],
    keyword_groups: list[tuple[str, list[str]]],
) -> FocusOutput:
    variants = _read_variants(
        sites_tsv,
        focus,
        samples,
        sample_abbrev,
        sv_size_threshold,
    )
    genes, exons, introns, cds, promoters = _read_gene_models(gff3, focus, promoter_bp)
    tracks = [_read_gff_track(path, label, focus) for label, path in track_gff]
    gene_data, feature_rows = _annotate_counts(
        genes=genes,
        exons=exons,
        introns=introns,
        cds=cds,
        promoters=promoters,
        tracks=tracks,
        variants=variants,
        annotations=annotations,
        focus=focus,
    )
    candidates = _select_candidates(feature_rows, focus, candidate_limit)
    _add_annotation_to_gene_data(gene_data, annotations)

    feature_tsv = outdir / f"{focus.slug}.features.tsv"
    _write_feature_tsv(feature_tsv, candidates)

    data = _build_browser_data(
        focus=focus,
        variants=variants,
        genes=gene_data,
        tracks=tracks,
        candidates=candidates,
        samples=samples,
        sample_abbrev=sample_abbrev,
        sites_tsv=sites_tsv,
        gff3=gff3,
        functional_tsv=functional_tsv,
        track_gff=track_gff,
        promoter_bp=promoter_bp,
        sv_size_threshold=sv_size_threshold,
        title=title,
        subtitle=subtitle,
        extraction_summary=extraction_summary,
        extraction_filters=extraction_filters,
        keyword_groups=keyword_groups,
    )
    html_path = outdir / f"{focus.slug}.html"
    html_path.write_text(render_focus_html(data), encoding="utf-8")

    metadata_path = outdir / f"{focus.slug}.json"
    metadata_path.write_text(json.dumps(data["summary"], indent=2) + "\n", encoding="utf-8")
    return FocusOutput(
        region=focus,
        html=html_path,
        features_tsv=feature_tsv,
        metadata_json=metadata_path,
        sites_tsv=sites_tsv,
        extraction_summary=extraction_summary,
    )


def _open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def _parse_attrs(raw: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for item in raw.split(";"):
        if not item:
            continue
        if "=" not in item:
            attrs[unquote(item)] = ""
            continue
        key, value = item.split("=", 1)
        attrs[unquote(key)] = unquote(value)
    return attrs


def _short_gene_id(gene_id: str) -> str:
    return re.sub(r"\.Wm82\.a6\.v1$", "", gene_id)


def _clip_interval(region: FocusRegion, start: int, end: int) -> tuple[int, int] | None:
    clipped_start = max(region.start, start)
    clipped_end = min(region.end, end)
    if clipped_start > clipped_end:
        return None
    return clipped_start, clipped_end


def _merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not intervals:
        return []
    merged = [sorted(intervals)[0]]
    for start, end in sorted(intervals)[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end + 1:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))
    return merged


def _read_function_annotations(path: Path | None) -> dict[str, dict[str, str]]:
    annotations: dict[str, dict[str, str]] = {}
    if path is None:
        return annotations
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            for key in (row.get("gene"), row.get("gene_id"), row.get("locus"), row.get("Name")):
                if key:
                    annotations[key] = row
    return annotations


def _read_gene_models(
    path: Path,
    region: FocusRegion,
    promoter_bp: int,
) -> GeneModelReadResult:
    genes_raw: dict[str, dict[str, Any]] = {}
    transcripts: dict[str, Transcript] = {}
    tx_by_gene: dict[str, list[str]] = defaultdict(list)
    exons_by_parent: dict[str, list[tuple[int, int]]] = defaultdict(list)
    cds_by_parent: dict[str, list[tuple[int, int]]] = defaultdict(list)

    with _open_text(path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            seqid, _source, feature_type, start_s, end_s, _score, strand, _phase, attrs_s = fields
            if seqid != region.contig:
                continue
            start, end = int(start_s), int(end_s)
            if end < region.start or start > region.end:
                continue
            attrs = _parse_attrs(attrs_s)
            if feature_type == "gene":
                gene_id = attrs.get("ID", "")
                if not gene_id:
                    continue
                genes_raw[gene_id] = {
                    "gene_id": gene_id,
                    "gene": attrs.get("Name", _short_gene_id(gene_id)),
                    "start": start,
                    "end": end,
                    "strand": strand,
                }
            elif feature_type in {"mRNA", "transcript"}:
                tx_id = attrs.get("ID", "")
                gene_id = attrs.get("Parent", "")
                if not tx_id or not gene_id:
                    continue
                transcripts[tx_id] = Transcript(
                    transcript_id=tx_id,
                    gene_id=gene_id,
                    start=start,
                    end=end,
                    strand=strand,
                    longest=attrs.get("longest") == "1",
                )
                tx_by_gene[gene_id].append(tx_id)
            elif feature_type in {"exon", "CDS"}:
                parents = [parent for parent in attrs.get("Parent", "").split(",") if parent]
                for parent in parents:
                    if feature_type == "exon":
                        exons_by_parent[parent].append((start, end))
                    else:
                        cds_by_parent[parent].append((start, end))

    models: list[GeneModel] = []
    exons: list[dict[str, Any]] = []
    introns: list[dict[str, Any]] = []
    cds: list[dict[str, Any]] = []
    promoters: list[dict[str, Any]] = []
    for gene_id, raw in sorted(
        genes_raw.items(),
        key=lambda item: (item[1]["start"], item[1]["end"]),
    ):
        chosen = _choose_transcript(gene_id, tx_by_gene, transcripts)
        parent_id = chosen.transcript_id if chosen else gene_id
        exon_intervals = _merge_intervals(exons_by_parent.get(parent_id, []))
        cds_intervals = _merge_intervals(cds_by_parent.get(parent_id, []))
        if not exon_intervals and chosen:
            exon_intervals = _merge_intervals(exons_by_parent.get(gene_id, []))
        if not cds_intervals and chosen:
            cds_intervals = _merge_intervals(cds_by_parent.get(gene_id, []))

        intron_intervals: list[tuple[int, int]] = []
        for left, right in zip(exon_intervals, exon_intervals[1:], strict=False):
            intron = _clip_interval(region, left[1] + 1, right[0] - 1)
            if intron:
                intron_intervals.append(intron)
        promoter = (
            _clip_interval(region, raw["start"] - promoter_bp, raw["start"] - 1)
            if raw["strand"] == "+"
            else _clip_interval(region, raw["end"] + 1, raw["end"] + promoter_bp)
        )
        model = GeneModel(
            gene_id=gene_id,
            gene=raw["gene"],
            start=raw["start"],
            end=raw["end"],
            strand=raw["strand"],
            transcript_id=chosen.transcript_id if chosen else "",
            exons=exon_intervals,
            cds=cds_intervals,
            introns=intron_intervals,
            promoter=promoter,
        )
        models.append(model)
        for idx, (start, end) in enumerate(exon_intervals, start=1):
            exons.append(
                {
                    "id": f"{model.gene}:exon:{idx}",
                    "gene": model.gene,
                    "start": start,
                    "end": end,
                }
            )
        for idx, (start, end) in enumerate(intron_intervals, start=1):
            introns.append(
                {
                    "id": f"{model.gene}:intron:{idx}",
                    "gene": model.gene,
                    "start": start,
                    "end": end,
                }
            )
        for idx, (start, end) in enumerate(cds_intervals, start=1):
            cds.append(
                {
                    "id": f"{model.gene}:CDS:{idx}",
                    "gene": model.gene,
                    "start": start,
                    "end": end,
                }
            )
        if promoter:
            promoters.append({
                "id": f"{model.gene}:promoter:{promoter_bp}bp",
                "gene": model.gene,
                "start": promoter[0],
                "end": promoter[1],
            })
    return models, exons, introns, cds, promoters


def _choose_transcript(
    gene_id: str,
    tx_by_gene: dict[str, list[str]],
    transcripts: dict[str, Transcript],
) -> Transcript | None:
    records = [transcripts[tx_id] for tx_id in tx_by_gene.get(gene_id, []) if tx_id in transcripts]
    if not records:
        return None
    records.sort(key=lambda tx: (not tx.longest, -(tx.end - tx.start), tx.transcript_id))
    return records[0]


def _read_gff_track(path: Path, label: str, region: FocusRegion) -> dict[str, Any]:
    features: list[dict[str, Any]] = []
    with _open_text(path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            seqid, source, feature_type, start_s, end_s, score, strand, _phase, attrs_s = fields
            if seqid != region.contig:
                continue
            start, end = int(start_s), int(end_s)
            if end < region.start or start > region.end:
                continue
            clipped = _clip_interval(region, start, end)
            if clipped is None:
                continue
            attrs = _parse_attrs(attrs_s)
            features.append({
                "id": attrs.get("ID", f"{label}:{len(features) + 1}"),
                "track": label,
                "source": source,
                "type": feature_type,
                "start": clipped[0],
                "end": clipped[1],
                "score": score,
                "strand": strand,
                "name": attrs.get("Name", attrs.get("ID", "")),
                "class": attrs.get("class", attrs.get("Class", feature_type)),
            })
    return {"label": label, "features": features}


def _read_variants(
    path: Path,
    region: FocusRegion,
    samples: tuple[str, str, str],
    sample_abbrev: dict[str, str],
    sv_size_threshold: int,
) -> list[dict[str, Any]]:
    variants: list[dict[str, Any]] = []
    display = {sample: sample_abbrev.get(sample, sample) for sample in samples}
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path} has no header.")
        missing = [
            column
            for column in ("contig", "pos", "ref", "alt")
            if column not in reader.fieldnames
        ]
        missing.extend(
            f"{sample}_gt"
            for sample in samples
            if f"{sample}_gt" not in reader.fieldnames
        )
        if missing:
            raise ValueError(f"{path} missing required column(s): {', '.join(missing)}")
        for row in reader:
            if row["contig"] != region.contig:
                continue
            pos = int(row["pos"])
            if pos < region.start or pos > region.end:
                continue
            ref = row.get("ref", "")
            alt = row.get("alt", "")
            variant_type = _classify_variant(row, sv_size_threshold)
            pattern = _classify_pattern(row, samples)
            overlapping_gene = row.get("overlapping_gene_names", "")
            nearest_gene = row.get("nearest_gene_name", "")
            variants.append({
                "p": pos,
                "e": pos + max(1, len(ref)) - 1,
                "type": variant_type,
                "pat": pattern,
                "target": pattern == "target_private",
                "ref": _truncate(ref),
                "alt": _truncate(alt),
                "gts": {display[sample]: row.get(f"{sample}_gt", "") for sample in samples},
                "gene": "" if overlapping_gene == "NA" else overlapping_gene,
                "nearest": "" if nearest_gene == "NA" else nearest_gene,
            })
    return sorted(variants, key=lambda item: item["p"])


def _classify_variant(row: dict[str, str], sv_size_threshold: int) -> str:
    ref = row.get("ref", "")
    alt = row.get("alt", "")
    alts = alt.split(",") if alt else []
    if row.get("variant_type", "").lower() == "snp":
        return "SNP"
    if any(item.startswith("<") or "[" in item or "]" in item for item in alts):
        return "SV-like"
    max_len = max([len(ref)] + [len(item) for item in alts] + [1])
    if max_len >= sv_size_threshold:
        return "SV-like"
    if any(len(item) != len(ref) for item in alts):
        return "INDEL"
    return "complex"


def _classify_pattern(row: dict[str, str], samples: tuple[str, str, str]) -> str:
    if row.get("target_private_alt_pattern", "").lower() == "true":
        return "target_private"
    off_target, derived, donor = samples
    off_gt = row.get(f"{off_target}_gt", "")
    derived_gt = row.get(f"{derived}_gt", "")
    donor_gt = row.get(f"{donor}_gt", "")
    if derived_gt and donor_gt and off_gt and derived_gt == donor_gt and derived_gt != off_gt:
        return "target_private"
    if derived_gt and off_gt and donor_gt and derived_gt == off_gt and derived_gt != donor_gt:
        return "background_like"
    if derived_gt and derived_gt == donor_gt == off_gt:
        return "all_same"
    return "other"


def _truncate(text: str, limit: int = 80) -> str:
    return text if len(text) <= limit else text[:limit] + "..."


def _annotate_counts(
    *,
    genes: list[GeneModel],
    exons: list[dict[str, Any]],
    introns: list[dict[str, Any]],
    cds: list[dict[str, Any]],
    promoters: list[dict[str, Any]],
    tracks: list[dict[str, Any]],
    variants: list[dict[str, Any]],
    annotations: dict[str, dict[str, str]],
    focus: FocusRegion,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    positions = [variant["p"] for variant in variants]
    rows: list[dict[str, Any]] = []

    def add_feature(
        feature_type: str,
        feature_id: str,
        label: str,
        start: int,
        end: int,
        gene: str = "",
        extra: str = "",
    ) -> dict[str, Any]:
        counts = _count_variants(positions, variants, start, end)
        annot = annotations.get(gene, {})
        row: dict[str, Any] = {
            "feature_type": feature_type,
            "feature_id": feature_id,
            "label": label,
            "gene": gene,
            "contig": focus.contig,
            "start": start,
            "end": end,
            "target_private_total": counts.get("target_total", 0),
            "target_private_snp": counts.get("target_SNP", 0),
            "target_private_indel_or_complex": (
                counts.get("target_INDEL", 0) + counts.get("target_complex", 0)
            ),
            "target_private_sv_like": counts.get("target_SV-like", 0),
            "all_variant_total": counts.get("all_total", 0),
            "extra": extra,
            "functional_category": annot.get("functional_category_keywords", ""),
            "screening_priority": annot.get("screening_priority", ""),
            "representative_function": (
                annot.get("representative_predicted_function")
                or annot.get("product")
                or annot.get("description")
                or ""
            ),
            "screening_note": annot.get("screening_note", ""),
        }
        row["score"] = (
            row["target_private_snp"]
            + 2 * row["target_private_indel_or_complex"]
            + 4 * row["target_private_sv_like"]
        )
        rows.append(row)
        return row

    gene_count_rows: dict[str, dict[str, Any]] = {}
    for gene in genes:
        gene_count_rows[gene.gene] = add_feature(
            "gene",
            gene.gene_id,
            gene.gene,
            gene.start,
            gene.end,
            gene.gene,
        )
    for feature in promoters:
        add_feature(
            "promoter",
            feature["id"],
            feature["id"],
            feature["start"],
            feature["end"],
            feature["gene"],
            "strand-aware upstream",
        )
    for feature in exons:
        add_feature(
            "exon",
            feature["id"],
            feature["id"],
            feature["start"],
            feature["end"],
            feature["gene"],
        )
    for feature in introns:
        add_feature(
            "intron",
            feature["id"],
            feature["id"],
            feature["start"],
            feature["end"],
            feature["gene"],
        )
    for feature in cds:
        add_feature(
            "CDS",
            feature["id"],
            feature["id"],
            feature["start"],
            feature["end"],
            feature["gene"],
        )
    for track in tracks:
        for feature in track["features"]:
            add_feature(
                track["label"],
                feature["id"],
                feature.get("name") or feature["id"],
                feature["start"],
                feature["end"],
                "",
                feature.get("class", track["label"]),
            )

    gene_data = []
    for gene in genes:
        counts = gene_count_rows.get(gene.gene, {})
        gene_data.append({
            "gene": gene.gene,
            "id": gene.gene_id,
            "start": gene.start,
            "end": gene.end,
            "strand": gene.strand,
            "tx": gene.transcript_id,
            "exons": [{"start": start, "end": end} for start, end in gene.exons],
            "cds": [{"start": start, "end": end} for start, end in gene.cds],
            "introns": [{"start": start, "end": end} for start, end in gene.introns],
            "promoter": (
                {"start": gene.promoter[0], "end": gene.promoter[1]}
                if gene.promoter
                else None
            ),
            "targetTotal": counts.get("target_private_total", 0),
            "targetSNP": counts.get("target_private_snp", 0),
            "targetIndel": counts.get("target_private_indel_or_complex", 0),
            "targetSV": counts.get("target_private_sv_like", 0),
        })
    return gene_data, rows


def _count_variants(
    positions: list[int],
    variants: list[dict[str, Any]],
    start: int,
    end: int,
) -> dict[str, int]:
    counts = Counter()
    for variant in variants[bisect_left(positions, start):bisect_right(positions, end)]:
        if variant["target"]:
            counts["target_total"] += 1
            counts[f"target_{variant['type']}"] += 1
        counts["all_total"] += 1
        counts[f"all_{variant['type']}"] += 1
    return dict(counts)


def _select_candidates(
    rows: list[dict[str, Any]],
    focus: FocusRegion,
    limit: int,
) -> list[dict[str, Any]]:
    candidates: list[dict[str, Any]] = []
    for row in rows:
        if row["target_private_total"] <= 0:
            continue
        bonus = 0
        priority = str(row.get("screening_priority", "")).lower()
        if "high" in priority:
            bonus += 8
        elif "moderate" in priority:
            bonus += 3
        if row.get("functional_category"):
            bonus += 2
        if row["feature_type"] in {"promoter", "exon", "CDS"}:
            bonus += 2
        if row["feature_type"] not in {"gene", "promoter", "exon", "intron", "CDS"}:
            bonus += 1
        out = dict(row)
        out["rank_score"] = int(row["score"]) + bonus
        out["region"] = focus.label
        candidates.append(out)
    candidates.sort(
        key=lambda item: (
            -item["rank_score"],
            item["start"],
            item["feature_type"],
            item["feature_id"],
        )
    )
    return candidates[:limit]


def _add_annotation_to_gene_data(
    gene_data: list[dict[str, Any]],
    annotations: dict[str, dict[str, str]],
) -> None:
    for gene in gene_data:
        annot = annotations.get(str(gene["gene"]), annotations.get(str(gene["id"]), {}))
        gene["function"] = (
            annot.get("representative_predicted_function")
            or annot.get("product")
            or annot.get("description")
            or ""
        )
        gene["categories"] = annot.get("functional_category_keywords", "")
        gene["priority"] = annot.get("screening_priority", "")


def _write_feature_tsv(path: Path, candidates: list[dict[str, Any]]) -> None:
    fields = [
        "rank",
        "region",
        "feature_type",
        "feature_id",
        "label",
        "gene",
        "contig",
        "start",
        "end",
        "target_private_total",
        "target_private_snp",
        "target_private_indel_or_complex",
        "target_private_sv_like",
        "all_variant_total",
        "rank_score",
        "screening_priority",
        "functional_category",
        "representative_function",
        "screening_note",
        "extra",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for idx, row in enumerate(candidates, start=1):
            out = {field: row.get(field, "") for field in fields}
            out["rank"] = idx
            writer.writerow(out)


def _build_browser_data(
    *,
    focus: FocusRegion,
    variants: list[dict[str, Any]],
    genes: list[dict[str, Any]],
    tracks: list[dict[str, Any]],
    candidates: list[dict[str, Any]],
    samples: tuple[str, str, str],
    sample_abbrev: dict[str, str],
    sites_tsv: Path,
    gff3: Path,
    functional_tsv: Path | None,
    track_gff: list[tuple[str, Path]],
    promoter_bp: int,
    sv_size_threshold: int,
    title: str | None,
    subtitle: str | None,
    extraction_summary: ExtractionSummary | None,
    extraction_filters: dict[str, object],
    keyword_groups: list[tuple[str, list[str]]],
) -> dict[str, Any]:
    target_counts = Counter(variant["type"] for variant in variants if variant["target"])
    type_counts = Counter(variant["type"] for variant in variants)
    groups = _candidate_groups(candidates, keyword_groups)
    return {
        "summary": {
            "title": title or f"Privy Interactive Focus: {focus.label}",
            "subtitle": subtitle or "Self-contained focus-region browser.",
            "contig": focus.contig,
            "start": focus.start,
            "end": focus.end,
            "length_bp": focus.length,
            "recommended_focus_bp": RECOMMENDED_FOCUS_BP,
            "promoter_bp": promoter_bp,
            "sv_size_threshold": sv_size_threshold,
            "samples": list(samples),
            "sample_display": {sample: sample_abbrev.get(sample, sample) for sample in samples},
            "gene_count": len(genes),
            "variant_count": len(variants),
            "variant_counts": dict(type_counts),
            "target_variant_counts": dict(target_counts),
            "track_counts": {track["label"]: len(track["features"]) for track in tracks},
            "inputs": {
                "sites_tsv": str(sites_tsv),
                "gff3": str(gff3),
                "functional_tsv": str(functional_tsv) if functional_tsv else "",
                "track_gff": [f"{label}={path}" for label, path in track_gff],
            },
            "extraction_filters": extraction_filters,
            "keyword_groups": [
                {"name": name, "terms": terms}
                for name, terms in keyword_groups
            ],
            "extraction_summary": (
                {
                    "records_seen": extraction_summary.records_seen,
                    "records_written": extraction_summary.records_written,
                    "skipped_filter": extraction_summary.skipped_filter,
                    "skipped_variant_type": extraction_summary.skipped_variant_type,
                    "skipped_missing": extraction_summary.skipped_missing,
                    "skipped_biallelic": extraction_summary.skipped_biallelic,
                }
                if extraction_summary is not None
                else {}
            ),
        },
        "genes": genes,
        "tracks": tracks,
        "variants": variants,
        "bins": _make_bins(focus, variants),
        "candidate_groups": groups,
    }


def _candidate_groups(
    candidates: list[dict[str, Any]],
    keyword_groups: list[tuple[str, list[str]]],
) -> list[dict[str, Any]]:
    gene_model_features = {"gene", "exon", "intron", "CDS"}
    built_in_features = {*gene_model_features, "promoter"}
    base = [
        ("all", "All Variant-Supported Features", lambda row: True),
        (
            "gene_models",
            "Gene / Exon / Intron / CDS",
            lambda row: row["feature_type"] in gene_model_features,
        ),
        ("promoters", "Promoters", lambda row: row["feature_type"] == "promoter"),
        (
            "extra_tracks",
            "Additional GFF3 Tracks",
            lambda row: row["feature_type"] not in built_in_features,
        ),
    ]
    groups = []
    used_keys = {key for key, _title, _predicate in base}
    for key, title, predicate in base:
        rows = [row for row in candidates if predicate(row)]
        groups.append({"key": key, "title": title, "rows": rows})
    for name, terms in keyword_groups:
        normalized_terms = [term.lower() for term in terms]
        rows = [
            row
            for row in candidates
            if _candidate_matches_terms(row, normalized_terms)
        ]
        groups.append(
            {
                "key": _unique_key(_slug(name), used_keys),
                "title": name,
                "rows": rows,
            }
        )
    return groups


def _candidate_matches_terms(row: dict[str, Any], terms: list[str]) -> bool:
    text = " ".join(
        str(row.get(field, ""))
        for field in (
            "feature_type",
            "label",
            "gene",
            "functional_category",
            "screening_priority",
            "representative_function",
            "screening_note",
            "extra",
        )
    ).lower()
    return any(term in text for term in terms)


def _slug(text: str) -> str:
    slug = re.sub(r"[^a-z0-9]+", "_", text.lower()).strip("_")
    return slug or "keyword_group"


def _unique_key(base: str, used: set[str]) -> str:
    key = base
    suffix = 2
    while key in used:
        key = f"{base}_{suffix}"
        suffix += 1
    used.add(key)
    return key


def _make_bins(
    region: FocusRegion,
    variants: list[dict[str, Any]],
    bin_bp: int = 50_000,
) -> list[dict[str, int]]:
    bins: list[dict[str, int]] = []
    for start in range(region.start, region.end + 1, bin_bp):
        bins.append({
            "start": start,
            "end": min(region.end, start + bin_bp - 1),
            "target": 0,
            "background": 0,
            "all_same": 0,
        })
    if not bins:
        return bins
    for variant in variants:
        idx = min(len(bins) - 1, max(0, (variant["p"] - region.start) // bin_bp))
        if variant["target"]:
            bins[idx]["target"] += 1
        elif variant["pat"] == "background_like":
            bins[idx]["background"] += 1
        elif variant["pat"] == "all_same":
            bins[idx]["all_same"] += 1
    return bins


def _write_run_metadata(
    path: Path,
    outputs: list[FocusOutput],
    sites_tsv: Path | None,
    vcf: Path | None,
    gff3: Path,
    track_gff: list[tuple[str, Path]],
) -> None:
    payload = {
        "analysis": "interactive_focus",
        "inputs": {
            "sites_tsv": str(sites_tsv) if sites_tsv is not None else "",
            "vcf": str(vcf) if vcf is not None else "",
            "gff3": str(gff3),
            "track_gff": [f"{label}={track_path}" for label, track_path in track_gff],
        },
        "outputs": [
            {
                "region": output.region.label,
                "html": str(output.html),
                "features_tsv": str(output.features_tsv),
                "metadata_json": str(output.metadata_json),
                "sites_tsv": str(output.sites_tsv) if output.sites_tsv else "",
            }
            for output in outputs
        ],
    }
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
