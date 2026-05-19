"""Build interactive dashboards from ``privy landscape`` output directories."""

from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from privy.interactive.landscape_render import render_landscape_html


@dataclass(frozen=True)
class LandscapePaths:
    path: Path
    windows: Path
    sample_windows: Path
    candidate_blocks: Path | None
    background_blocks: Path | None
    filter_summary: Path | None
    similarity: Path | None
    run_json: Path | None


def run_landscape_dashboard(
    *,
    landscape_dir: Path,
    outdir: Path,
    title: str | None = None,
    subtitle: str | None = None,
    max_windows: int = 20000,
    max_sample_windows: int = 80000,
    max_blocks: int = 5000,
) -> list[Path]:
    """Write a self-contained HTML dashboard for existing ``privy landscape`` outputs."""
    if not landscape_dir.exists():
        raise FileNotFoundError(f"--landscape path not found: {landscape_dir}")
    if not landscape_dir.is_dir():
        raise ValueError(f"--landscape must be a directory: {landscape_dir}")
    if max_windows < 1:
        raise ValueError("--max-windows must be at least 1.")
    if max_sample_windows < 1:
        raise ValueError("--max-sample-windows must be at least 1.")
    if max_blocks < 1:
        raise ValueError("--max-blocks must be at least 1.")

    paths = _landscape_paths(landscape_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    data = _build_dashboard_data(
        paths=paths,
        title=title,
        subtitle=subtitle,
        max_windows=max_windows,
        max_sample_windows=max_sample_windows,
        max_blocks=max_blocks,
    )

    html_path = outdir / "landscape_dashboard.html"
    metadata_path = outdir / "landscape_dashboard.json"
    html_path.write_text(render_landscape_html(data), encoding="utf-8")
    metadata_path.write_text(json.dumps(data["summary"], indent=2) + "\n", encoding="utf-8")
    return [html_path, metadata_path]


def _landscape_paths(path: Path) -> LandscapePaths:
    windows = path / "windows.tsv"
    sample_windows = path / "sample_windows.tsv"
    if not windows.exists():
        raise FileNotFoundError(f"Landscape directory is missing windows.tsv: {path}")
    if not sample_windows.exists():
        raise FileNotFoundError(f"Landscape directory is missing sample_windows.tsv: {path}")
    return LandscapePaths(
        path=path,
        windows=windows,
        sample_windows=sample_windows,
        candidate_blocks=_existing(path / "candidate_introgression_blocks.tsv"),
        background_blocks=_existing(path / "background_blocks.tsv"),
        filter_summary=_existing(path / "filter_summary.tsv"),
        similarity=_existing(path / "similarity.tsv"),
        run_json=_existing(path / "landscape.json"),
    )


def _existing(path: Path) -> Path | None:
    return path if path.exists() else None


def _build_dashboard_data(
    *,
    paths: LandscapePaths,
    title: str | None,
    subtitle: str | None,
    max_windows: int,
    max_sample_windows: int,
    max_blocks: int,
) -> dict[str, Any]:
    run_meta = _read_json(paths.run_json)
    windows, window_summary = _read_windows(paths.windows, max_windows)
    sample_windows, sample_summary = _read_sample_windows(paths.sample_windows, max_sample_windows)
    candidate_blocks, candidate_count = _read_candidate_blocks(paths.candidate_blocks, max_blocks)
    background_summary, background_count = _read_background_blocks(paths.background_blocks)
    filters = _read_filter_summary(paths.filter_summary)
    similarity_pairs = _read_similarity_pairs(paths.similarity)

    contigs = sorted(window_summary["contigs"])
    samples = sorted(sample_summary["samples"])
    return {
        "summary": {
            "title": title or "Privy Interactive Landscape Dashboard",
            "subtitle": subtitle
            or "Self-contained dashboard from existing privy landscape outputs.",
            "analysis": "interactive_landscape",
            "landscape_dir": str(paths.path),
            "window_count": window_summary["window_count"],
            "sample_window_count": sample_summary["sample_window_count"],
            "candidate_block_count": candidate_count,
            "background_block_count": background_count,
            "embedded_windows": len(windows),
            "embedded_sample_windows": len(sample_windows),
            "embedded_candidate_blocks": len(candidate_blocks),
            "max_windows": max_windows,
            "max_sample_windows": max_sample_windows,
            "max_blocks": max_blocks,
            "contigs": contigs,
            "samples": samples,
            "cohort_roles": sample_summary["roles"],
            "parameters": run_meta.get("parameters", {}),
            "inputs": {
                "windows": str(paths.windows),
                "sample_windows": str(paths.sample_windows),
                "candidate_blocks": str(paths.candidate_blocks) if paths.candidate_blocks else "",
                "background_blocks": (
                    str(paths.background_blocks) if paths.background_blocks else ""
                ),
                "filter_summary": str(paths.filter_summary) if paths.filter_summary else "",
                "similarity": str(paths.similarity) if paths.similarity else "",
                "landscape_json": str(paths.run_json) if paths.run_json else "",
            },
        },
        "windows": windows,
        "sample_windows": sample_windows,
        "candidate_blocks": candidate_blocks,
        "background_summary": background_summary,
        "filter_summary": filters,
        "similarity_pairs": similarity_pairs,
        "contig_summary": _contig_summary(windows),
        "sample_summary": sample_summary["sample_rows"],
    }


def _read_windows(path: Path, limit: int) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    contigs: set[str] = set()
    count = 0
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for raw in reader:
            count += 1
            contig = raw.get("contig", "")
            contigs.add(contig)
            if len(rows) < limit:
                rows.append({
                    "window_id": raw.get("window_id", ""),
                    "contig": contig,
                    "window_index": _int(raw.get("window_index")),
                    "start": _int(raw.get("start")),
                    "end": _int(raw.get("end")),
                    "midpoint": _int(raw.get("midpoint")),
                    "n_variants": _int(raw.get("n_variants")),
                    "density_variants_per_kb": _float(raw.get("density_variants_per_kb")),
                    "target_mean_missing_rate": _float(raw.get("target_mean_missing_rate")),
                    "offtarget_mean_missing_rate": _float(raw.get("offtarget_mean_missing_rate")),
                    "target_private_alt_n": _int(raw.get("target_private_alt_n")),
                    "offtarget_private_alt_n": _int(raw.get("offtarget_private_alt_n")),
                    "target_private_alt_rate": _float(raw.get("target_private_alt_rate")),
                    "offtarget_private_alt_rate": _float(raw.get("offtarget_private_alt_rate")),
                    "top_nearest_background": raw.get("top_nearest_background", ""),
                })
    return rows, {"window_count": count, "contigs": contigs}


def _read_sample_windows(path: Path, limit: int) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    samples: set[str] = set()
    roles: dict[str, str] = {}
    sample_counts: dict[str, int] = defaultdict(int)
    sample_sums: dict[str, dict[str, float]] = defaultdict(lambda: defaultdict(float))
    count = 0
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for raw in reader:
            count += 1
            sample = raw.get("sample", "")
            role = raw.get("cohort_role", "")
            samples.add(sample)
            roles[sample] = role
            sample_counts[sample] += 1
            metrics = {
                "missing_rate": _float(raw.get("missing_rate")),
                "nonref_rate": _float(raw.get("nonref_rate")),
                "private_alt_rate": _float(raw.get("private_alt_rate")),
                "nearest_similarity": _float(raw.get("nearest_similarity")),
            }
            for metric, value in metrics.items():
                sample_sums[sample][metric] += value
            if len(rows) < limit:
                rows.append({
                    "window_id": raw.get("window_id", ""),
                    "contig": raw.get("contig", ""),
                    "window_index": _int(raw.get("window_index")),
                    "start": _int(raw.get("start")),
                    "end": _int(raw.get("end")),
                    "sample": sample,
                    "cohort_role": role,
                    "missing_rate": metrics["missing_rate"],
                    "nonref_rate": metrics["nonref_rate"],
                    "private_alt_rate": metrics["private_alt_rate"],
                    "nearest_background": raw.get("nearest_background", ""),
                    "nearest_background_role": raw.get("nearest_background_role", ""),
                    "nearest_similarity": metrics["nearest_similarity"],
                })
    sample_rows = []
    for sample in sorted(samples):
        n = max(1, sample_counts[sample])
        sample_rows.append({
            "sample": sample,
            "cohort_role": roles.get(sample, ""),
            "mean_missing_rate": sample_sums[sample]["missing_rate"] / n,
            "mean_nonref_rate": sample_sums[sample]["nonref_rate"] / n,
            "mean_private_alt_rate": sample_sums[sample]["private_alt_rate"] / n,
            "mean_nearest_similarity": sample_sums[sample]["nearest_similarity"] / n,
            "n_windows": sample_counts[sample],
        })
    return rows, {
        "sample_window_count": count,
        "samples": samples,
        "roles": roles,
        "sample_rows": sample_rows,
    }


def _read_candidate_blocks(
    path: Path | None,
    limit: int,
) -> tuple[list[dict[str, Any]], int]:
    if path is None:
        return [], 0
    rows: list[dict[str, Any]] = []
    count = 0
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for raw in reader:
            count += 1
            if len(rows) < limit:
                rows.append({
                    "block_id": raw.get("block_id", ""),
                    "sample": raw.get("sample", ""),
                    "contig": raw.get("contig", ""),
                    "start": _int(raw.get("start")),
                    "end": _int(raw.get("end")),
                    "n_windows": _int(raw.get("n_windows")),
                    "candidate_donor": raw.get("candidate_donor", ""),
                    "candidate_donor_role": raw.get("candidate_donor_role", ""),
                    "mean_donor_similarity": _float(raw.get("mean_donor_similarity")),
                    "mean_nearest_target_similarity": _float(
                        raw.get("mean_nearest_target_similarity")
                    ),
                    "mean_similarity_delta": _float(raw.get("mean_similarity_delta")),
                    "max_missing_rate": _float(raw.get("max_missing_rate")),
                    "mean_private_alt_rate": _float(raw.get("mean_private_alt_rate")),
                    "evidence_class": raw.get("evidence_class", ""),
                    "interpretation": raw.get("interpretation", ""),
                })
    return rows, count


def _read_background_blocks(path: Path | None) -> tuple[list[dict[str, Any]], int]:
    if path is None:
        return [], 0
    counts = Counter[str]()
    count = 0
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for raw in reader:
            count += 1
            sample = raw.get("sample", "")
            background = raw.get("nearest_background", "unassigned")
            role = raw.get("nearest_background_role", "")
            counts.update([f"{sample} -> {background} ({role})"])
    return _counter_rows(counts, limit=50), count


def _read_filter_summary(path: Path | None) -> list[dict[str, str]]:
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


def _read_similarity_pairs(path: Path | None) -> list[dict[str, Any]]:
    if path is None:
        return []
    sums: dict[str, float] = defaultdict(float)
    counts: dict[str, int] = defaultdict(int)
    compared: dict[str, int] = defaultdict(int)
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for raw in reader:
            a = raw.get("sample_a", "")
            b = raw.get("sample_b", "")
            if not a or not b:
                continue
            key = " / ".join(sorted([a, b]))
            sums[key] += _float(raw.get("similarity"))
            counts[key] += 1
            compared[key] += _int(raw.get("compared_variants"))
    rows = [
        {
            "pair": key,
            "mean_similarity": sums[key] / max(1, counts[key]),
            "n_windows": counts[key],
            "compared_variants": compared[key],
        }
        for key in sorted(sums)
    ]
    rows.sort(key=lambda row: (-float(str(row["mean_similarity"])), str(row["pair"])))
    return rows


def _read_json(path: Path | None) -> dict[str, Any]:
    if path is None:
        return {}
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return {}
    return payload if isinstance(payload, dict) else {}


def _contig_summary(windows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[str, dict[str, float]] = defaultdict(lambda: defaultdict(float))
    counts: Counter[str] = Counter()
    for row in windows:
        contig = str(row["contig"])
        counts.update([contig])
        grouped[contig]["n_variants"] += float(row["n_variants"])
        grouped[contig]["target_private_alt_n"] += float(row["target_private_alt_n"])
        grouped[contig]["offtarget_private_alt_n"] += float(row["offtarget_private_alt_n"])
        grouped[contig]["density_variants_per_kb"] += float(row["density_variants_per_kb"])
    output = []
    for contig in sorted(grouped):
        n = max(1, counts[contig])
        output.append({
            "contig": contig,
            "n_windows": counts[contig],
            "n_variants": int(grouped[contig]["n_variants"]),
            "target_private_alt_n": int(grouped[contig]["target_private_alt_n"]),
            "offtarget_private_alt_n": int(grouped[contig]["offtarget_private_alt_n"]),
            "mean_density_variants_per_kb": grouped[contig]["density_variants_per_kb"] / n,
        })
    return output


def _counter_rows(counter: Counter[str], limit: int | None = None) -> list[dict[str, Any]]:
    items = sorted(counter.items(), key=lambda item: (-item[1], item[0]))
    rows = [{"label": label, "count": count} for label, count in items]
    if limit is None or len(rows) <= limit:
        return rows
    kept = rows[:limit]
    kept.append({"label": "other", "count": sum(count for _label, count in items[limit:])})
    return kept


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
