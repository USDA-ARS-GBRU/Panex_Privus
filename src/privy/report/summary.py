"""TSV summary and ranked-hit table generators for privy report.

Reads hits.tsv (required) plus optional regions, evidence, qc, and run.json
outputs from a previous ``privy scan`` run, then writes a set of summary TSVs
and a Markdown (with optional HTML) report to the output directory.

Output files written here:
    summary.tsv              — run-level key/value summary
    ranked_hits.tsv          — top-N hits with explicit rank column
    strictness_summary.tsv   — count and percentage per strictness class
    support_summary.tsv      — evidence breakdown by source and class (if evidence.tsv given)
    contradiction_summary.tsv — contradiction metrics from qc.tsv / compare.tsv
    report.md                — human-readable Markdown report
    report.html              — HTML version (when fmt includes 'html')
"""

from __future__ import annotations

import logging
from collections import Counter
from pathlib import Path
from typing import Any

from privy.core.config import PrivyConfig
from privy.io.jsonio import read_run_json
from privy.io.tsv import (
    HITS_COLUMNS,
    QC_COLUMNS,
    RANKED_HITS_COLUMNS,
    STRICTNESS_SUMMARY_COLUMNS,
    SUPPORT_SUMMARY_COLUMNS,
    TsvWriter,
    read_tsv,
)
from privy.report.html import render_html_report
from privy.report.markdown import render_markdown_report
from privy.utils.misc import now_iso

log = logging.getLogger("privy.report.summary")

# Canonical order for the strictness class table — contradicted is never in
# hits.tsv but listed for completeness.
_STRICTNESS_ORDER = [
    "strict_complete",
    "strict_target_missing",
    "strict_offtarget_missing",
    "strict_both_missing",
    "relaxed_threshold",
]


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def run_report(
    hits: Path,
    regions: Path | None,
    evidence: Path | None,
    compare: Path | None,
    qc: Path | None,
    run_json: Path | None,
    cfg: PrivyConfig,
    fmt: str,
    top_n: int,
    include_qc: bool,
    include_strictness: bool,
    include_compare: bool,
    include_regions: bool,
    title: str,
    outdir: Path,
) -> None:
    """Generate ranked summaries and a human-readable report.

    Args:
        hits:               Path to ``hits.tsv`` (required).
        regions:            Path to ``regions.tsv`` (optional).
        evidence:           Path to ``evidence.tsv`` (optional).
        compare:            Path to ``compare.tsv`` from ``privy compare`` (optional).
        qc:                 Path to ``qc.tsv`` (optional).
        run_json:           Path to ``run.json`` (optional).
        cfg:                Resolved configuration.
        fmt:                Output format — ``"markdown"``, ``"html"``, or ``"both"``.
        top_n:              Number of top loci to include in ranked_hits.tsv.
        include_qc:         Include QC section in report.
        include_strictness: Include strictness distribution section.
        include_compare:    Include support/compare section.
        include_regions:    Include candidate region section.
        title:              Report title.
        outdir:             Output directory (must already exist).
    """
    # ── Load required input ────────────────────────────────────────────────────
    hit_rows: list[dict[str, str]] = read_tsv(hits)
    if not hit_rows:
        log.warning("hits.tsv is empty — report will contain no hit information.")

    # ── Load optional inputs ──────────────────────────────────────────────────
    region_rows: list[dict[str, str]] | None = (
        read_tsv(regions) if (regions is not None and regions.exists()) else None
    )
    evidence_rows: list[dict[str, str]] | None = (
        read_tsv(evidence) if (evidence is not None and evidence.exists()) else None
    )
    qc_rows: list[dict[str, str]] | None = (
        read_tsv(qc) if (qc is not None and qc.exists()) else None
    )
    compare_rows: list[dict[str, str]] | None = (
        read_tsv(compare) if (compare is not None and compare.exists()) else None
    )
    run_meta: dict[str, Any] | None = (
        read_run_json(run_json) if (run_json is not None and run_json.exists()) else None
    )

    # ── Compute summaries ─────────────────────────────────────────────────────
    ranked_hits = _rank_hits(hit_rows, top_n)
    strictness_rows = _compute_strictness_summary(hit_rows)
    support_rows: list[dict[str, Any]] | None = (
        _compute_support_summary(evidence_rows) if evidence_rows is not None else None
    )
    contradiction_rows = _compute_contradiction_summary(qc_rows, compare_rows)
    run_summary_rows = _compute_run_summary(
        hit_rows, region_rows, qc_rows, run_meta, cfg
    )

    # ── Write TSV summaries ──────────────────────────────────────────────────
    _write_summary_tsv(run_summary_rows, outdir)
    _write_ranked_hits_tsv(ranked_hits, outdir)
    if include_strictness:
        _write_strictness_summary_tsv(strictness_rows, outdir)
    if support_rows is not None and include_compare:
        _write_support_summary_tsv(support_rows, outdir)
    if contradiction_rows:
        _write_contradiction_summary_tsv(contradiction_rows, outdir)

    # ── Assemble sections for the narrative report ────────────────────────────
    sections: dict[str, Any] = {
        "title": title,
        "generated_at": now_iso(),
        "run_summary_rows": run_summary_rows,
        "ranked_hits": ranked_hits,
        "n_total_hits": len(hit_rows),
        "top_n": top_n,
        "strictness_summary": strictness_rows if include_strictness else [],
        "qc_rows": qc_rows if include_qc else None,
        "region_rows": region_rows if include_regions else None,
        "support_summary": (
            support_rows if (include_compare and support_rows is not None) else None
        ),
        "contradiction_rows": contradiction_rows,
    }

    # ── Render report ─────────────────────────────────────────────────────────
    md_path = render_markdown_report(sections, title, outdir)
    if fmt in ("html", "both"):
        render_html_report(md_path, outdir)

    log.info("Report written to %s (format=%s)", outdir, fmt)


# ---------------------------------------------------------------------------
# Sorting and ranking
# ---------------------------------------------------------------------------

def _rank_hits(
    hit_rows: list[dict[str, str]],
    top_n: int,
) -> list[dict[str, str]]:
    """Sort *hit_rows* by ``final_score`` descending and return the top *top_n*."""
    def _score(row: dict[str, str]) -> float:
        try:
            return float(row.get("final_score") or "0")
        except ValueError:
            return 0.0

    return sorted(hit_rows, key=_score, reverse=True)[:top_n]


# ---------------------------------------------------------------------------
# Summary computations
# ---------------------------------------------------------------------------

def _compute_strictness_summary(
    hit_rows: list[dict[str, str]],
) -> list[dict[str, Any]]:
    """Count hits per strictness class and compute percentage of total.

    Returns rows in canonical order (strict_complete first).  Classes that are
    present in the data but not in the canonical order are appended at the end.
    """
    total = len(hit_rows)
    counts: Counter[str] = Counter(
        r.get("strictness_class", "unknown") for r in hit_rows
    )

    result: list[dict[str, Any]] = []
    seen: set[str] = set()
    for cls in _STRICTNESS_ORDER:
        n = counts.get(cls, 0)
        pct = round(100.0 * n / total, 1) if total > 0 else 0.0
        result.append({"strictness_class": cls, "n_loci": n, "pct_hits": pct})
        seen.add(cls)

    # Append any unexpected classes (shouldn't arise from normal scan output)
    for cls, n in sorted(counts.items()):
        if cls not in seen:
            pct = round(100.0 * n / total, 1) if total > 0 else 0.0
            result.append({"strictness_class": cls, "n_loci": n, "pct_hits": pct})

    return result


def _compute_support_summary(
    evidence_rows: list[dict[str, str]],
) -> list[dict[str, Any]]:
    """Group evidence.tsv rows by source_type × evidence_class and compute counts."""
    counts: Counter[tuple[str, str]] = Counter()
    source_totals: Counter[str] = Counter()

    for row in evidence_rows:
        src = row.get("source_type", "unknown")
        cls = row.get("evidence_class", "unknown")
        counts[(src, cls)] += 1
        source_totals[src] += 1

    result: list[dict[str, Any]] = []
    for (src, cls), n in sorted(counts.items()):
        total = source_totals[src]
        pct = round(100.0 * n / total, 1) if total > 0 else 0.0
        result.append({
            "source_type": src,
            "evidence_class": cls,
            "n_records": n,
            "pct_of_source": pct,
        })
    return result


def _compute_contradiction_summary(
    qc_rows: list[dict[str, str]] | None,
    compare_rows: list[dict[str, str]] | None,
) -> list[dict[str, str]]:
    """Collect contradiction metrics from qc.tsv and optionally compare.tsv."""
    rows: list[dict[str, str]] = []

    if qc_rows is not None:
        qc_by_metric = {r["metric"]: r for r in qc_rows if "metric" in r}
        if "alleles_contradicted" in qc_by_metric:
            rows.append({
                "metric": "alleles_contradicted",
                "value": qc_by_metric["alleles_contradicted"].get("value", "0"),
                "description": "Alleles where off-targets also carry the allele",
            })

    if compare_rows is not None:
        n_contradicted = sum(
            1 for r in compare_rows if r.get("match_class") == "contradicted"
        )
        rows.append({
            "metric": "compare_contradicted_loci",
            "value": str(n_contradicted),
            "description": "Loci classified as contradicted in privy compare output",
        })

    return rows


def _compute_run_summary(
    hit_rows: list[dict[str, str]],
    region_rows: list[dict[str, str]] | None,
    qc_rows: list[dict[str, str]] | None,
    run_meta: dict[str, Any] | None,
    cfg: PrivyConfig,
) -> list[dict[str, str]]:
    """Build run-level summary rows for summary.tsv."""
    rows: list[dict[str, str]] = []

    project = cfg.project_name or ""
    if not project and run_meta:
        project = str(run_meta.get("project_name", ""))
    rows.append({
        "metric": "project_name",
        "value": project,
        "description": "Project name from configuration",
    })

    rows.append({
        "metric": "n_hits",
        "value": str(len(hit_rows)),
        "description": "Total loci in hits.tsv",
    })

    if region_rows is not None:
        rows.append({
            "metric": "n_regions",
            "value": str(len(region_rows)),
            "description": "Total candidate regions in regions.tsv",
        })

    if hit_rows:
        top = max(
            hit_rows,
            key=lambda r: float(r.get("final_score") or "0"),
            default=None,
        )
        if top is not None:
            rows.append({
                "metric": "top_locus_id",
                "value": top.get("locus_id", ""),
                "description": "Locus ID with highest final_score",
            })
            rows.append({
                "metric": "top_final_score",
                "value": top.get("final_score", ""),
                "description": "Highest final_score among all hits",
            })

    if qc_rows is not None:
        qc_vals = {r["metric"]: r["value"] for r in qc_rows if "metric" in r}
        for key in ("records_evaluated", "alleles_passed", "alleles_contradicted"):
            if key in qc_vals:
                rows.append({
                    "metric": key,
                    "value": qc_vals[key],
                    "description": f"From qc.tsv: {key.replace('_', ' ')}",
                })

    if run_meta is not None:
        for key in ("start_time", "end_time", "privy_version"):
            val = run_meta.get(key)
            if val is not None:
                rows.append({
                    "metric": key,
                    "value": str(val),
                    "description": f"From run.json: {key.replace('_', ' ')}",
                })

    return rows


# ---------------------------------------------------------------------------
# TSV writers
# ---------------------------------------------------------------------------

def _write_summary_tsv(rows: list[dict[str, str]], outdir: Path) -> None:
    with TsvWriter(outdir / "summary.tsv", QC_COLUMNS) as w:
        w.write_rows(rows)


def _write_ranked_hits_tsv(
    ranked_hits: list[dict[str, str]], outdir: Path
) -> None:
    with TsvWriter(outdir / "ranked_hits.tsv", RANKED_HITS_COLUMNS) as w:
        for rank, row in enumerate(ranked_hits, start=1):
            out: dict[str, Any] = {"rank": rank}
            for col in HITS_COLUMNS:
                out[col] = row.get(col, "")
            w.write_row(out)


def _write_strictness_summary_tsv(
    rows: list[dict[str, Any]], outdir: Path
) -> None:
    with TsvWriter(outdir / "strictness_summary.tsv", STRICTNESS_SUMMARY_COLUMNS) as w:
        for row in rows:
            w.write_row({
                "strictness_class": row["strictness_class"],
                "n_loci": row["n_loci"],
                "pct_hits": row["pct_hits"],
            })


def _write_support_summary_tsv(
    rows: list[dict[str, Any]], outdir: Path
) -> None:
    with TsvWriter(outdir / "support_summary.tsv", SUPPORT_SUMMARY_COLUMNS) as w:
        for row in rows:
            w.write_row({
                "source_type": row["source_type"],
                "evidence_class": row["evidence_class"],
                "n_records": row["n_records"],
                "pct_of_source": row["pct_of_source"],
            })


def _write_contradiction_summary_tsv(
    rows: list[dict[str, str]], outdir: Path
) -> None:
    with TsvWriter(outdir / "contradiction_summary.tsv", QC_COLUMNS) as w:
        w.write_rows(rows)
