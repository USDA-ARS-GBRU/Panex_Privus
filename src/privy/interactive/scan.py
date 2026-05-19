"""Build interactive dashboards from ``privy scan`` output directories."""

from __future__ import annotations

import csv
import json
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from privy.interactive.branding import default_dashboard_title
from privy.interactive.scan_render import render_scan_html


@dataclass(frozen=True)
class ScanSourcePaths:
    key: str
    label: str
    path: Path
    hits: Path
    regions: Path | None = None
    evidence: Path | None = None
    qc: Path | None = None
    run_json: Path | None = None


SCORE_BINS: tuple[tuple[float, str], ...] = (
    (0.25, "<0.25"),
    (0.50, "0.25-0.50"),
    (0.75, "0.50-0.75"),
    (1.00, "0.75-1.00"),
    (1.25, "1.00-1.25"),
    (1.50, "1.25-1.50"),
)


def run_scan_dashboard(
    *,
    scan_dir: Path,
    outdir: Path,
    title: str | None = None,
    subtitle: str | None = None,
    max_hits: int = 5000,
    max_regions: int = 1000,
) -> list[Path]:
    """Write a self-contained HTML dashboard for existing ``privy scan`` outputs."""
    if not scan_dir.exists():
        raise FileNotFoundError(f"--scan path not found: {scan_dir}")
    if not scan_dir.is_dir():
        raise ValueError(f"--scan must be a directory: {scan_dir}")
    if max_hits < 1:
        raise ValueError("--max-hits must be at least 1.")
    if max_regions < 1:
        raise ValueError("--max-regions must be at least 1.")

    sources = _discover_scan_sources(scan_dir)
    if not sources:
        raise ValueError(
            "--scan must point to a directory containing hits.tsv or source "
            "subdirectories such as vcf/ and gfa/."
        )

    outdir.mkdir(parents=True, exist_ok=True)
    data = _build_dashboard_data(
        scan_dir=scan_dir,
        sources=sources,
        title=title,
        subtitle=subtitle,
        max_hits=max_hits,
        max_regions=max_regions,
    )

    html_path = outdir / "scan_dashboard.html"
    metadata_path = outdir / "scan_dashboard.json"
    html_path.write_text(render_scan_html(data), encoding="utf-8")
    metadata_path.write_text(json.dumps(data["summary"], indent=2) + "\n", encoding="utf-8")
    return [html_path, metadata_path]


def _discover_scan_sources(scan_dir: Path) -> list[ScanSourcePaths]:
    if (scan_dir / "hits.tsv").exists():
        return [_source_paths(scan_dir, scan_dir.name or "scan")]

    preferred = ["vcf", "gfa", "xmfa"]
    found: list[ScanSourcePaths] = []
    seen: set[Path] = set()
    for name in preferred:
        child = scan_dir / name
        if (child / "hits.tsv").exists():
            found.append(_source_paths(child, name.upper()))
            seen.add(child.resolve())
    for child in sorted(path for path in scan_dir.iterdir() if path.is_dir()):
        if child.resolve() in seen:
            continue
        if (child / "hits.tsv").exists():
            found.append(_source_paths(child, child.name))
    return found


def _source_paths(path: Path, label: str) -> ScanSourcePaths:
    return ScanSourcePaths(
        key=_slug(label),
        label=label,
        path=path,
        hits=path / "hits.tsv",
        regions=_existing(path / "regions.tsv"),
        evidence=_existing(path / "evidence.tsv"),
        qc=_existing(path / "qc.tsv"),
        run_json=_existing(path / "run.json"),
    )


def _existing(path: Path) -> Path | None:
    return path if path.exists() else None


def _build_dashboard_data(
    *,
    scan_dir: Path,
    sources: list[ScanSourcePaths],
    title: str | None,
    subtitle: str | None,
    max_hits: int,
    max_regions: int,
) -> dict[str, Any]:
    source_payloads = [
        _read_scan_source(source, max_hits=max_hits, max_regions=max_regions)
        for source in sources
    ]
    compare = _read_compare(scan_dir)
    hit_rows = [
        row
        for source in source_payloads
        for row in source["top_hits"]
    ]
    region_rows = [
        row
        for source in source_payloads
        for row in source["top_regions"]
    ]
    total_hits = sum(int(source["summary"]["hit_count"]) for source in source_payloads)
    total_regions = sum(int(source["summary"]["region_count"]) for source in source_payloads)
    strictness = Counter[str]()
    variant_types = Counter[str]()
    contigs = Counter[str]()
    scores = Counter[str]()
    for source in source_payloads:
        strictness.update(_counter_from_rows(source["strictness_counts"]))
        variant_types.update(_counter_from_rows(source["variant_type_counts"]))
        contigs.update(_counter_from_rows(source["contig_counts"]))
        scores.update(_counter_from_rows(source["score_bins"]))
    return {
        "summary": {
            "title": title or default_dashboard_title("Scan"),
            "subtitle": subtitle or "Self-contained dashboard from existing privy scan outputs.",
            "analysis": "interactive_scan",
            "scan_dir": str(scan_dir),
            "source_count": len(source_payloads),
            "hit_count": total_hits,
            "region_count": total_regions,
            "embedded_hit_rows": len(hit_rows),
            "embedded_region_rows": len(region_rows),
            "max_hits_per_source": max_hits,
            "max_regions_per_source": max_regions,
            "inputs": [
                {
                    "label": source.label,
                    "path": str(source.path),
                    "hits": str(source.hits),
                    "regions": str(source.regions) if source.regions else "",
                    "evidence": str(source.evidence) if source.evidence else "",
                    "qc": str(source.qc) if source.qc else "",
                    "run_json": str(source.run_json) if source.run_json else "",
                }
                for source in sources
            ],
        },
        "sources": source_payloads,
        "hits": hit_rows,
        "regions": region_rows,
        "compare": compare,
        "aggregate": {
            "strictness_counts": _counter_rows(strictness),
            "variant_type_counts": _counter_rows(variant_types),
            "contig_counts": _counter_rows(contigs, limit=40),
            "score_bins": _ordered_score_bins(scores),
        },
    }


def _read_scan_source(
    source: ScanSourcePaths,
    *,
    max_hits: int,
    max_regions: int,
) -> dict[str, Any]:
    strictness = Counter[str]()
    variant_types = Counter[str]()
    contigs = Counter[str]()
    scores = Counter[str]()
    top_hits: list[dict[str, Any]] = []
    hit_count = 0
    max_score: float | None = None
    min_score: float | None = None

    with source.hits.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            hit_count += 1
            strictness.update([_clean_label(row.get("strictness_class"), "unknown")])
            variant_types.update([_clean_label(row.get("variant_type"), "unknown")])
            contigs.update([_clean_label(row.get("contig"), "unknown")])
            score = _float(row.get("final_score"))
            scores.update([_score_bin(score)])
            max_score = score if max_score is None else max(max_score, score)
            min_score = score if min_score is None else min(min_score, score)
            if len(top_hits) < max_hits:
                top_hits.append(_hit_payload(row, source, hit_count, score))

    top_regions: list[dict[str, Any]] = []
    region_count = 0
    if source.regions is not None:
        with source.regions.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                region_count += 1
                if len(top_regions) < max_regions:
                    top_regions.append(_region_payload(row, source, region_count))

    qc_rows = _read_qc(source.qc)
    evidence_counts = _read_evidence_counts(source.evidence)
    return {
        "key": source.key,
        "label": source.label,
        "summary": {
            "hit_count": hit_count,
            "region_count": region_count,
            "embedded_hit_rows": len(top_hits),
            "embedded_region_rows": len(top_regions),
            "min_score": min_score if min_score is not None else 0.0,
            "max_score": max_score if max_score is not None else 0.0,
        },
        "strictness_counts": _counter_rows(strictness),
        "variant_type_counts": _counter_rows(variant_types),
        "contig_counts": _counter_rows(contigs, limit=40),
        "score_bins": _ordered_score_bins(scores),
        "evidence_counts": _counter_rows(evidence_counts),
        "qc": qc_rows,
        "top_hits": top_hits,
        "top_regions": top_regions,
    }


def _hit_payload(
    row: dict[str, str],
    source: ScanSourcePaths,
    rank: int,
    score: float,
) -> dict[str, Any]:
    return {
        "source": source.key,
        "source_label": source.label,
        "rank": rank,
        "locus_id": row.get("locus_id", ""),
        "contig": row.get("contig", ""),
        "start": _int(row.get("start")),
        "end": _int(row.get("end")),
        "variant_type": row.get("variant_type", ""),
        "allele_key": row.get("allele_key", ""),
        "target_support_n": _int(row.get("target_support_n")),
        "target_total_n": _int(row.get("target_total_n")),
        "offtarget_support_n": _int(row.get("offtarget_support_n")),
        "offtarget_total_n": _int(row.get("offtarget_total_n")),
        "target_missing_n": _int(row.get("target_missing_n")),
        "offtarget_missing_n": _int(row.get("offtarget_missing_n")),
        "strictness_class": row.get("strictness_class", ""),
        "discovery_score": _float(row.get("discovery_score")),
        "support_score": _float(row.get("support_score")),
        "penalty_score": _float(row.get("penalty_score")),
        "final_score": score,
    }


def _region_payload(
    row: dict[str, str],
    source: ScanSourcePaths,
    rank: int,
) -> dict[str, Any]:
    return {
        "source": source.key,
        "source_label": source.label,
        "rank": rank,
        "region_id": row.get("region_id", ""),
        "contig": row.get("contig", ""),
        "start": _int(row.get("start")),
        "end": _int(row.get("end")),
        "n_loci": _int(row.get("n_loci")),
        "variant_types": row.get("variant_types", ""),
        "dominant_strictness_class": row.get("dominant_strictness_class", ""),
        "target_consistency": _float(row.get("target_consistency")),
        "offtarget_exclusion": _float(row.get("offtarget_exclusion")),
        "final_score": _float(row.get("final_score")),
    }


def _read_qc(path: Path | None) -> list[dict[str, str]]:
    if path is None:
        return []
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [
            {
                "metric": row.get("metric", ""),
                "value": row.get("value", ""),
                "description": row.get("description", ""),
            }
            for row in reader
        ]


def _read_evidence_counts(path: Path | None) -> Counter[str]:
    counts = Counter[str]()
    if path is None:
        return counts
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            source_type = _clean_label(row.get("source_type"), "source")
            evidence_class = _clean_label(row.get("evidence_class"), "unknown")
            counts.update([f"{source_type}:{evidence_class}"])
    return counts


def _read_compare(scan_dir: Path) -> dict[str, Any]:
    compare_dir = scan_dir / "compare"
    summary_path = compare_dir / "compare_summary.tsv"
    compare_path = compare_dir / "compare.tsv"
    summary_rows: list[dict[str, str]] = []
    match_counts = Counter[str]()
    if summary_path.exists():
        with summary_path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            summary_rows = [dict(row) for row in reader]
    if compare_path.exists():
        with compare_path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                match_counts.update([_clean_label(row.get("match_class"), "unknown")])
    return {
        "summary_rows": summary_rows,
        "match_counts": _counter_rows(match_counts),
        "compare_tsv": str(compare_path) if compare_path.exists() else "",
        "compare_summary_tsv": str(summary_path) if summary_path.exists() else "",
    }


def _counter_from_rows(rows: object) -> Counter[str]:
    counts = Counter[str]()
    if not isinstance(rows, list):
        return counts
    for row in rows:
        if not isinstance(row, dict):
            continue
        label = row.get("label")
        count = row.get("count")
        if isinstance(label, str) and isinstance(count, int):
            counts[label] += count
    return counts


def _counter_rows(counter: Counter[str], limit: int | None = None) -> list[dict[str, Any]]:
    items = sorted(counter.items(), key=lambda item: (-item[1], item[0]))
    rows = [
        {"label": label, "count": count}
        for label, count in items
    ]
    if limit is None or len(rows) <= limit:
        return rows
    kept = rows[:limit]
    other = sum(count for _label, count in items[limit:])
    kept.append({"label": "other", "count": other})
    return kept


def _ordered_score_bins(counter: Counter[str]) -> list[dict[str, Any]]:
    labels = [label for _upper, label in SCORE_BINS] + [">=1.50"]
    return [{"label": label, "count": counter.get(label, 0)} for label in labels]


def _score_bin(value: float) -> str:
    for upper, label in SCORE_BINS:
        if value < upper:
            return label
    return ">=1.50"


def _clean_label(value: str | None, fallback: str) -> str:
    if value is None:
        return fallback
    value = value.strip()
    return value or fallback


def _float(value: str | None) -> float:
    try:
        return float(value or 0)
    except ValueError:
        return 0.0


def _int(value: str | None) -> int:
    try:
        return int(float(value or 0))
    except ValueError:
        return 0


def _slug(value: str) -> str:
    slug = re.sub(r"[^a-z0-9]+", "_", value.lower()).strip("_")
    return slug or "scan"
