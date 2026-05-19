"""Build interactive dashboards from ``privy pangenome`` output directories."""

from __future__ import annotations

import csv
import json
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from privy.interactive.branding import default_dashboard_title
from privy.interactive.pangenome_render import render_pangenome_html


@dataclass(frozen=True)
class PangenomeSourcePaths:
    key: str
    label: str
    path: Path
    feature_summary: Path
    composition: Path | None = None
    coverage: Path | None = None
    growth: Path | None = None
    run_json: Path | None = None


def run_pangenome_dashboard(
    *,
    pangenome_dir: Path,
    outdir: Path,
    title: str | None = None,
    subtitle: str | None = None,
    max_features: int = 10000,
    max_private_features: int = 5000,
) -> list[Path]:
    """Write a self-contained HTML dashboard for existing ``privy pangenome`` outputs."""
    if not pangenome_dir.exists():
        raise FileNotFoundError(f"--pangenome path not found: {pangenome_dir}")
    if not pangenome_dir.is_dir():
        raise ValueError(f"--pangenome must be a directory: {pangenome_dir}")
    if max_features < 1:
        raise ValueError("--max-features must be at least 1.")
    if max_private_features < 1:
        raise ValueError("--max-private-features must be at least 1.")

    sources = _discover_pangenome_sources(pangenome_dir)
    if not sources:
        raise ValueError(
            "--pangenome must point to a directory containing feature_summary.tsv "
            "or source subdirectories such as vcf/ and gfa/."
        )

    outdir.mkdir(parents=True, exist_ok=True)
    data = _build_dashboard_data(
        pangenome_dir=pangenome_dir,
        sources=sources,
        title=title,
        subtitle=subtitle,
        max_features=max_features,
        max_private_features=max_private_features,
    )

    html_path = outdir / "pangenome_dashboard.html"
    metadata_path = outdir / "pangenome_dashboard.json"
    html_path.write_text(render_pangenome_html(data), encoding="utf-8")
    metadata_path.write_text(json.dumps(data["summary"], indent=2) + "\n", encoding="utf-8")
    return [html_path, metadata_path]


def _discover_pangenome_sources(pangenome_dir: Path) -> list[PangenomeSourcePaths]:
    if (pangenome_dir / "feature_summary.tsv").exists():
        return [_source_paths(pangenome_dir, pangenome_dir.name or "pangenome")]

    preferred = ["vcf", "gfa"]
    found: list[PangenomeSourcePaths] = []
    seen: set[Path] = set()
    for name in preferred:
        child = pangenome_dir / name
        if (child / "feature_summary.tsv").exists():
            found.append(_source_paths(child, name.upper()))
            seen.add(child.resolve())
    for child in sorted(path for path in pangenome_dir.iterdir() if path.is_dir()):
        if child.resolve() in seen:
            continue
        if (child / "feature_summary.tsv").exists():
            found.append(_source_paths(child, child.name))
    return found


def _source_paths(path: Path, label: str) -> PangenomeSourcePaths:
    return PangenomeSourcePaths(
        key=_slug(label),
        label=label,
        path=path,
        feature_summary=path / "feature_summary.tsv",
        composition=_existing(path / "composition.tsv"),
        coverage=_existing(path / "coverage_histogram.tsv"),
        growth=_existing(path / "growth_curves.tsv"),
        run_json=_existing(path / "pangenome.json"),
    )


def _existing(path: Path) -> Path | None:
    return path if path.exists() else None


def _build_dashboard_data(
    *,
    pangenome_dir: Path,
    sources: list[PangenomeSourcePaths],
    title: str | None,
    subtitle: str | None,
    max_features: int,
    max_private_features: int,
) -> dict[str, Any]:
    source_payloads = [
        _read_pangenome_source(
            source,
            max_features=max_features,
            max_private_features=max_private_features,
        )
        for source in sources
    ]
    feature_rows = [
        row
        for source in source_payloads
        for row in source["features"]
    ]
    private_rows = [
        row
        for source in source_payloads
        for row in source["target_private_features"]
    ]
    total_features = sum(int(source["summary"]["feature_count"]) for source in source_payloads)
    target_private = sum(
        int(source["summary"]["target_private_count"]) for source in source_payloads
    )
    offtarget_private = sum(
        int(source["summary"]["offtarget_private_count"]) for source in source_payloads
    )
    return {
        "summary": {
            "title": title or default_dashboard_title("Pangenome"),
            "subtitle": subtitle
            or "Self-contained dashboard from existing privy pangenome outputs.",
            "analysis": "interactive_pangenome",
            "pangenome_dir": str(pangenome_dir),
            "source_count": len(source_payloads),
            "feature_count": total_features,
            "target_private_count": target_private,
            "offtarget_private_count": offtarget_private,
            "embedded_features": len(feature_rows),
            "embedded_target_private_features": len(private_rows),
            "max_features_per_source": max_features,
            "max_private_features_per_source": max_private_features,
            "inputs": [
                {
                    "label": source.label,
                    "path": str(source.path),
                    "feature_summary": str(source.feature_summary),
                    "composition": str(source.composition) if source.composition else "",
                    "coverage": str(source.coverage) if source.coverage else "",
                    "growth": str(source.growth) if source.growth else "",
                    "pangenome_json": str(source.run_json) if source.run_json else "",
                }
                for source in sources
            ],
        },
        "sources": source_payloads,
        "features": feature_rows,
        "target_private_features": private_rows,
        "aggregate": _aggregate_sources(source_payloads),
    }


def _read_pangenome_source(
    source: PangenomeSourcePaths,
    *,
    max_features: int,
    max_private_features: int,
) -> dict[str, Any]:
    run_meta = _read_json(source.run_json)
    feature_rows: list[dict[str, Any]] = []
    private_rows: list[dict[str, Any]] = []
    feature_count = 0
    target_private_count = 0
    offtarget_private_count = 0
    feature_types = Counter[str]()
    full_categories = Counter[str]()
    target_categories = Counter[str]()
    offtarget_categories = Counter[str]()
    contigs = Counter[str]()
    total_bp = 0
    target_private_bp = 0

    with source.feature_summary.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for raw in reader:
            row = _feature_payload(raw, source)
            feature_count += 1
            total_bp += int(row["length"])
            feature_types.update([str(row["feature_type"] or "unknown")])
            full_categories.update([str(row["full_category"] or "unknown")])
            target_categories.update([str(row["target_category"] or "unknown")])
            offtarget_categories.update([str(row["offtarget_category"] or "unknown")])
            if row["contig"]:
                contigs.update([str(row["contig"])])
            if row["target_private"]:
                target_private_count += 1
                target_private_bp += int(row["length"])
                if len(private_rows) < max_private_features:
                    private_rows.append(row)
            if row["offtarget_private"]:
                offtarget_private_count += 1
            if len(feature_rows) < max_features:
                feature_rows.append(row)

    private_rows.sort(
        key=lambda item: (
            -int(item["length"]),
            str(item["contig"]),
            int(item["start"]),
            str(item["feature_id"]),
        )
    )
    composition = _read_composition(source.composition)
    coverage = _read_coverage(source.coverage)
    growth = _read_growth(source.growth)
    return {
        "key": source.key,
        "label": source.label,
        "summary": {
            "feature_count": feature_count,
            "target_private_count": target_private_count,
            "offtarget_private_count": offtarget_private_count,
            "total_bp": total_bp,
            "target_private_bp": target_private_bp,
            "source_type": run_meta.get("source_type", ""),
            "n_samples": run_meta.get("summary", {}).get("n_samples", 0),
            "n_target_samples": run_meta.get("summary", {}).get("n_target_samples", 0),
            "n_offtarget_samples": run_meta.get("summary", {}).get("n_offtarget_samples", 0),
            "parameters": run_meta.get("parameters", {}),
            "samples": run_meta.get("samples", {}),
            "embedded_features": len(feature_rows),
            "embedded_target_private_features": len(private_rows),
        },
        "feature_type_counts": _counter_rows(feature_types),
        "full_category_counts": _counter_rows(full_categories),
        "target_category_counts": _counter_rows(target_categories),
        "offtarget_category_counts": _counter_rows(offtarget_categories),
        "contig_counts": _counter_rows(contigs, limit=40),
        "composition": composition,
        "coverage": coverage,
        "growth": growth,
        "features": feature_rows,
        "target_private_features": private_rows,
    }


def _feature_payload(raw: dict[str, str], source: PangenomeSourcePaths) -> dict[str, Any]:
    return {
        "source": source.key,
        "source_label": source.label,
        "feature_id": raw.get("feature_id", ""),
        "source_type": raw.get("source_type", ""),
        "feature_type": raw.get("feature_type", ""),
        "contig": raw.get("contig", ""),
        "start": _int(raw.get("start")),
        "end": _int(raw.get("end")),
        "length": _int(raw.get("length")),
        "total_present_n": _int(raw.get("total_present_n")),
        "target_present_n": _int(raw.get("target_present_n")),
        "target_total_n": _int(raw.get("target_total_n")),
        "offtarget_present_n": _int(raw.get("offtarget_present_n")),
        "offtarget_total_n": _int(raw.get("offtarget_total_n")),
        "full_category": raw.get("full_category", ""),
        "target_category": raw.get("target_category", ""),
        "offtarget_category": raw.get("offtarget_category", ""),
        "target_private": _bool(raw.get("target_private")),
        "offtarget_private": _bool(raw.get("offtarget_private")),
    }


def _read_composition(path: Path | None) -> list[dict[str, Any]]:
    if path is None:
        return []
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [
            {
                "group": row.get("group", ""),
                "category": row.get("category", ""),
                "n_features": _int(row.get("n_features")),
                "n_bp": _int(row.get("n_bp")),
            }
            for row in reader
        ]


def _read_coverage(path: Path | None) -> list[dict[str, Any]]:
    if path is None:
        return []
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [
            {
                "group": row.get("group", ""),
                "coverage": _int(row.get("coverage")),
                "n_features": _int(row.get("n_features")),
                "n_bp": _int(row.get("n_bp")),
            }
            for row in reader
        ]


def _read_growth(path: Path | None) -> list[dict[str, Any]]:
    if path is None:
        return []
    grouped: dict[tuple[str, int], dict[str, float]] = defaultdict(lambda: defaultdict(float))
    counts: Counter[tuple[str, int]] = Counter()
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            key = (row.get("group", ""), _int(row.get("n")))
            counts.update([key])
            grouped[key]["features"] += _float(row.get("features"))
            grouped[key]["bp"] += _float(row.get("bp"))
            grouped[key]["new_features"] += _float(row.get("new_features"))
            grouped[key]["singleton_features"] += _float(row.get("singleton_features"))
    rows = []
    for (group, n), sums in sorted(grouped.items(), key=lambda item: (item[0][0], item[0][1])):
        denom = max(1, counts[(group, n)])
        rows.append({
            "group": group,
            "n": n,
            "mean_features": sums["features"] / denom,
            "mean_bp": sums["bp"] / denom,
            "mean_new_features": sums["new_features"] / denom,
            "mean_singleton_features": sums["singleton_features"] / denom,
        })
    return rows


def _aggregate_sources(sources: list[dict[str, Any]]) -> dict[str, Any]:
    feature_types = Counter[str]()
    full_categories = Counter[str]()
    target_categories = Counter[str]()
    offtarget_categories = Counter[str]()
    contigs = Counter[str]()
    composition: list[dict[str, Any]] = []
    coverage: list[dict[str, Any]] = []
    growth: list[dict[str, Any]] = []
    for source in sources:
        feature_types.update(_counter_from_rows(source["feature_type_counts"]))
        full_categories.update(_counter_from_rows(source["full_category_counts"]))
        target_categories.update(_counter_from_rows(source["target_category_counts"]))
        offtarget_categories.update(_counter_from_rows(source["offtarget_category_counts"]))
        contigs.update(_counter_from_rows(source["contig_counts"]))
        for row in source["composition"]:
            out = dict(row)
            out["source"] = source["key"]
            out["source_label"] = source["label"]
            composition.append(out)
        for row in source["coverage"]:
            out = dict(row)
            out["source"] = source["key"]
            out["source_label"] = source["label"]
            coverage.append(out)
        for row in source["growth"]:
            out = dict(row)
            out["source"] = source["key"]
            out["source_label"] = source["label"]
            growth.append(out)
    return {
        "feature_type_counts": _counter_rows(feature_types),
        "full_category_counts": _counter_rows(full_categories),
        "target_category_counts": _counter_rows(target_categories),
        "offtarget_category_counts": _counter_rows(offtarget_categories),
        "contig_counts": _counter_rows(contigs, limit=40),
        "composition": composition,
        "coverage": coverage,
        "growth": growth,
    }


def _read_json(path: Path | None) -> dict[str, Any]:
    if path is None:
        return {}
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return {}
    return payload if isinstance(payload, dict) else {}


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
    rows = [{"label": label, "count": count} for label, count in items]
    if limit is None or len(rows) <= limit:
        return rows
    kept = rows[:limit]
    kept.append({"label": "other", "count": sum(count for _label, count in items[limit:])})
    return kept


def _bool(value: str | None) -> bool:
    return str(value or "").strip().lower() in {"true", "1", "yes", "y"}


def _float(value: str | None) -> float:
    if value is None:
        return 0.0
    text = value.strip()
    if text in {"", "NA", "nan"}:
        return 0.0
    try:
        return float(text)
    except ValueError:
        return 0.0


def _int(value: str | None) -> int:
    if value is None:
        return 0
    text = value.strip()
    if text in {"", "NA", "nan"}:
        return 0
    try:
        return int(float(text))
    except ValueError:
        return 0


def _slug(value: str) -> str:
    slug = re.sub(r"[^a-z0-9]+", "_", value.lower()).strip("_")
    return slug or "pangenome"
