"""Cross-source comparison engine for ``privy compare``.

Compares two ``hits.tsv`` files produced by separate ``privy scan`` runs
(typically one VCF scan and one GFA scan) and classifies each locus pair
as supported, partially_supported, contradicted, source_specific,
uninformative, or missing_data.

Algorithm:
    1. Load source-A and source-B hits rows.
    2. Build a contig-bucketed index over source B, optionally canonicalizing
       minigraph-cactus names like ``SAMPLE#HAP#CONTIG`` to ``CONTIG``.
    3. For each source-A locus, find overlapping source-B loci using
       ``overlap_mode`` (``contained``, ``reciprocal``, or ``any``), or within
       ``breakpoint_tolerance_bp`` as a partial-support fallback.
    4. Classify the pair: SUPPORTED / PARTIALLY_SUPPORTED / CONTRADICTED.
    5. Emit unmatched source-A and source-B loci as SOURCE_SPECIFIC entries.
    6. Write compare.tsv, compare_summary.tsv, compare.json.
"""

from __future__ import annotations

import json
import logging
from bisect import bisect_right
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Literal

from privy.core.config import CompareConfig, PrivyConfig
from privy.core.evidence import MatchClass
from privy.io.tsv import (
    COMPARE_COLUMNS,
    COMPARE_SUMMARY_COLUMNS,
    TsvWriter,
    read_tsv,
)

log = logging.getLogger("privy.backends.compare")

_MATCH_SCORES: dict[MatchClass, float] = {
    MatchClass.SUPPORTED: 1.0,
    MatchClass.PARTIALLY_SUPPORTED: 0.5,
    MatchClass.CONTRADICTED: 0.0,
    MatchClass.SOURCE_SPECIFIC: 0.3,
    MatchClass.UNINFORMATIVE: 0.1,
    MatchClass.MISSING_DATA: 0.0,
}

_STRICT_CLASSES = frozenset(
    {
        "strict_complete",
        "strict_target_missing",
        "strict_offtarget_missing",
        "strict_both_missing",
    }
)


@dataclass
class HitsRow:
    """One row from a hits.tsv file with numeric types parsed."""

    locus_id: str
    contig: str
    start: int
    end: int
    variant_type: str
    allele_key: str
    strictness_class: str
    final_score: float
    raw: dict[str, str] = field(default_factory=dict, repr=False)


@dataclass(frozen=True)
class OverlapMetrics:
    """Interval overlap details for a candidate locus pair."""

    intersection_bp: int
    reciprocal: float
    containment: float
    query_fraction: float
    candidate_fraction: float


@dataclass(frozen=True)
class MatchResult:
    """One selected candidate for a source-A locus."""

    row: HitsRow
    overlap: float
    metrics: OverlapMetrics
    method: Literal["overlap", "breakpoint"]
    breakpoint_distance_bp: int = 0


@dataclass
class MatchSearchStats:
    """Counters collected while matching for compare.json diagnostics."""

    rows_a_with_candidate_contig: int = 0
    rows_a_without_candidate_contig: int = 0
    candidate_pairs_checked: int = 0
    overlapping_candidate_pairs: int = 0
    qualifying_overlap_pairs: int = 0
    breakpoint_candidate_pairs: int = 0
    breakpoint_fallback_matches: int = 0
    overlap_matches: int = 0


@dataclass(frozen=True)
class ContigCandidateIndex:
    """Sorted source-B rows for one canonical contig."""

    rows: list[HitsRow]
    starts: list[int]


def load_hits_tsv(path: Path) -> list[HitsRow]:
    """Parse a hits.tsv file into :class:`HitsRow` objects.

    Args:
        path: Path to a hits.tsv produced by ``privy scan``.

    Returns:
        List of parsed rows.

    Raises:
        FileNotFoundError: If *path* does not exist.
        ValueError: If a row is missing required columns or has unparseable values.
    """
    if not path.exists():
        raise FileNotFoundError(f"hits.tsv not found: {path}")
    raw_rows = read_tsv(path)
    result: list[HitsRow] = []
    for raw in raw_rows:
        try:
            result.append(
                HitsRow(
                    locus_id=raw["locus_id"],
                    contig=raw["contig"],
                    start=int(raw["start"]),
                    end=int(raw["end"]),
                    variant_type=raw.get("variant_type", ""),
                    allele_key=raw.get("allele_key", ""),
                    strictness_class=raw.get("strictness_class", ""),
                    final_score=float(raw.get("final_score", 0.0)),
                    raw=raw,
                )
            )
        except (KeyError, ValueError) as exc:
            raise ValueError(f"Malformed hits.tsv row in {path}: {raw!r}") from exc
    return result


def infer_source_label(rows: list[HitsRow], explicit: str | None) -> str:
    """Return the display label for an evidence source.

    Uses *explicit* when provided.  Otherwise infers from locus_id prefix:
    ``GPX`` → ``"gfa"``, anything else → ``"vcf"``.
    """
    if explicit:
        return explicit
    if rows and rows[0].locus_id.startswith("GPX"):
        return "gfa"
    return "vcf"


def canonicalize_contig(contig: str, *, normalize: bool = True) -> str:
    """Return the comparable contig name for a scan output row.

    Minigraph-cactus GFA walks commonly use ``SAMPLE#HAP#CONTIG`` names while
    VCF records use the bare reference contig.  With normalization enabled,
    compare uses the final ``#``-delimited field for matching.
    """
    if normalize and "#" in contig:
        return contig.rsplit("#", 1)[-1]
    return contig


def _same_contig(a: HitsRow, b: HitsRow, cfg: CompareConfig) -> bool:
    return canonicalize_contig(a.contig, normalize=cfg.normalize_contigs) == canonicalize_contig(
        b.contig, normalize=cfg.normalize_contigs
    )


def _interval_len(row: HitsRow) -> int:
    return max(0, row.end - row.start)


def _overlap_metrics(
    a: HitsRow,
    b: HitsRow,
    cfg: CompareConfig,
) -> OverlapMetrics:
    if not _same_contig(a, b, cfg):
        return OverlapMetrics(
            intersection_bp=0,
            reciprocal=0.0,
            containment=0.0,
            query_fraction=0.0,
            candidate_fraction=0.0,
        )

    ov_start = max(a.start, b.start)
    ov_end = min(a.end, b.end)
    intersection = max(0, ov_end - ov_start)
    if intersection == 0:
        return OverlapMetrics(
            intersection_bp=0,
            reciprocal=0.0,
            containment=0.0,
            query_fraction=0.0,
            candidate_fraction=0.0,
        )

    len_a = _interval_len(a)
    len_b = _interval_len(b)
    union = max(a.end, b.end) - min(a.start, b.start)
    reciprocal = intersection / union if union > 0 else 0.0
    query_fraction = intersection / len_a if len_a > 0 else 0.0
    candidate_fraction = intersection / len_b if len_b > 0 else 0.0
    containment = min(1.0, max(query_fraction, candidate_fraction))
    return OverlapMetrics(
        intersection_bp=intersection,
        reciprocal=reciprocal,
        containment=containment,
        query_fraction=query_fraction,
        candidate_fraction=candidate_fraction,
    )


def _overlap_score(metrics: OverlapMetrics, cfg: CompareConfig) -> float:
    if cfg.overlap_mode == "reciprocal":
        return metrics.reciprocal
    return metrics.containment


def _passes_overlap(metrics: OverlapMetrics, cfg: CompareConfig) -> bool:
    if metrics.intersection_bp <= 0:
        return False
    if cfg.overlap_mode == "any":
        return True
    return _overlap_score(metrics, cfg) >= cfg.min_reciprocal_overlap


def reciprocal_overlap_rows(
    a: HitsRow,
    b: HitsRow,
    cfg: CompareConfig | None = None,
) -> float:
    """Reciprocal overlap (intersection / union) between two HitsRow objects.

    Returns 0.0 if loci are on different contigs or do not overlap.
    """
    compare_cfg = cfg or CompareConfig(normalize_contigs=False)
    return _overlap_metrics(a, b, compare_cfg).reciprocal


def _breakpoint_distance(a: HitsRow, b: HitsRow, cfg: CompareConfig) -> int:
    """Minimum bp gap between two loci on the same contig."""
    if not _same_contig(a, b, cfg):
        return 10**9
    return max(0, max(a.start, b.start) - min(a.end, b.end))


def is_state_compatible(
    strictness_a: str,
    strictness_b: str,
    require: bool,
) -> bool:
    """Return True when two strictness classes are evidence-compatible.

    Without *require*: compatible if neither is ``"contradicted"``.
    With *require*: additionally, both must be in the same broad category
    (``strict_*`` vs ``relaxed_threshold``).
    """
    if strictness_a == "contradicted" or strictness_b == "contradicted":
        return False
    if not require:
        return True
    a_strict = strictness_a in _STRICT_CLASSES
    b_strict = strictness_b in _STRICT_CLASSES
    return a_strict == b_strict


def classify_match(
    overlap: float,
    state_compat: bool,
    strictness_a: str,
    strictness_b: str,
    cfg: CompareConfig,
    *,
    method: Literal["overlap", "breakpoint"] = "overlap",
) -> MatchClass:
    """Determine the MatchClass for a locus pair.

    Precedence:
        1. Either strictness is ``contradicted`` → CONTRADICTED
        2. No coordinate overlap but within breakpoint tolerance → PARTIALLY_SUPPORTED
        3. No coordinate overlap → SOURCE_SPECIFIC
        4. ``require_state_compatibility`` and states mismatch → PARTIALLY_SUPPORTED
        5. States compatible → SUPPORTED
        6. Otherwise → PARTIALLY_SUPPORTED
    """
    if strictness_a == "contradicted" or strictness_b == "contradicted":
        return MatchClass.CONTRADICTED
    if overlap <= 0.0:
        if method == "breakpoint":
            return MatchClass.PARTIALLY_SUPPORTED
        return MatchClass.SOURCE_SPECIFIC
    if cfg.require_state_compatibility and not state_compat:
        return MatchClass.PARTIALLY_SUPPORTED
    if state_compat:
        return MatchClass.SUPPORTED
    return MatchClass.PARTIALLY_SUPPORTED


def compute_comparison_score(match_class: MatchClass, overlap: float) -> float:
    """Numeric comparison score in [0, 1].

    For SUPPORTED and PARTIALLY_SUPPORTED the base score is scaled by
    the configured overlap score so higher-confidence matches score higher.
    """
    base = _MATCH_SCORES[match_class]
    if match_class in (MatchClass.SUPPORTED, MatchClass.PARTIALLY_SUPPORTED):
        return round(base * max(overlap, 0.01), 4)
    return base


def _build_contig_index(
    rows: list[HitsRow],
    cfg: CompareConfig,
) -> dict[str, ContigCandidateIndex]:
    grouped: dict[str, list[HitsRow]] = defaultdict(list)
    for r in rows:
        grouped[canonicalize_contig(r.contig, normalize=cfg.normalize_contigs)].append(r)

    idx: dict[str, ContigCandidateIndex] = {}
    for contig, contig_rows in grouped.items():
        sorted_rows = sorted(contig_rows, key=lambda row: (row.start, row.end))
        idx[contig] = ContigCandidateIndex(
            rows=sorted_rows,
            starts=[row.start for row in sorted_rows],
        )
    return idx


def _candidate_window(
    query: HitsRow,
    contig_index: ContigCandidateIndex,
    cfg: CompareConfig,
) -> list[HitsRow]:
    """Return source-B candidates close enough to overlap or be a breakpoint hit."""
    tolerance = cfg.breakpoint_tolerance_bp
    lower_bound = query.start - tolerance
    upper_bound = query.end + tolerance
    stop = bisect_right(contig_index.starts, upper_bound)
    return [row for row in contig_index.rows[:stop] if row.end >= lower_bound]


def _match_sort_key(result: MatchResult) -> tuple[float, float, int, float]:
    return (
        result.overlap,
        result.metrics.reciprocal,
        result.metrics.intersection_bp,
        result.row.final_score,
    )


def find_match_results(
    query: HitsRow,
    candidates: list[HitsRow],
    cfg: CompareConfig,
    stats: MatchSearchStats | None = None,
) -> list[MatchResult]:
    """Return candidate matches for a query locus.

    Primary search emits all candidates passing ``cfg.overlap_mode``.  Fallback uses
    ``cfg.breakpoint_tolerance_bp`` and returns a breakpoint match only when
    no coordinate-overlap candidate qualifies.

    Returns an empty list when no candidate qualifies.
    """
    results: list[MatchResult] = []

    for cand in candidates:
        if stats is not None:
            stats.candidate_pairs_checked += 1
        metrics = _overlap_metrics(query, cand, cfg)
        if metrics.intersection_bp > 0 and stats is not None:
            stats.overlapping_candidate_pairs += 1
        if not _passes_overlap(metrics, cfg):
            continue
        if stats is not None:
            stats.qualifying_overlap_pairs += 1
        result = MatchResult(
            row=cand,
            overlap=_overlap_score(metrics, cfg),
            metrics=metrics,
            method="overlap",
        )
        results.append(result)

    if results:
        results.sort(key=_match_sort_key, reverse=True)
        if stats is not None:
            stats.overlap_matches += len(results)
        return results

    if cfg.breakpoint_tolerance_bp > 0:
        best: MatchResult | None = None
        best_dist = cfg.breakpoint_tolerance_bp + 1
        for cand in candidates:
            dist = _breakpoint_distance(query, cand, cfg)
            if dist <= cfg.breakpoint_tolerance_bp and dist < best_dist:
                if stats is not None:
                    stats.breakpoint_candidate_pairs += 1
                metrics = _overlap_metrics(query, cand, cfg)
                best = MatchResult(
                    row=cand,
                    overlap=0.0,
                    metrics=metrics,
                    method="breakpoint",
                    breakpoint_distance_bp=dist,
                )
                best_dist = dist

        if best is not None:
            if stats is not None:
                stats.breakpoint_fallback_matches += 1
            return [best]

    return []


def find_best_match_result(
    query: HitsRow,
    candidates: list[HitsRow],
    cfg: CompareConfig,
    stats: MatchSearchStats | None = None,
) -> MatchResult | None:
    """Return the best candidate match for a query locus."""
    results = find_match_results(query, candidates, cfg, stats)
    if not results:
        return None
    return results[0]


def find_best_match(
    query: HitsRow,
    candidates: list[HitsRow],
    cfg: CompareConfig,
) -> tuple[HitsRow | None, float]:
    """Return ``(best_candidate, overlap)`` for a query locus."""
    result = find_best_match_result(query, candidates, cfg)
    if result is None:
        return None, 0.0
    return result.row, result.overlap


def _pair_to_row(
    compare_id: str,
    row_a: HitsRow,
    row_b: HitsRow,
    match_class: MatchClass,
    overlap: float,
    state_compat: bool,
    score: float,
    source_label_a: str,
    source_label_b: str,
    *,
    method: Literal["overlap", "breakpoint"] = "overlap",
    breakpoint_distance_bp: int = 0,
) -> dict[str, str]:
    method_detail = f"method={method}"
    if method == "breakpoint":
        method_detail += f"; distance_bp={breakpoint_distance_bp}"

    if match_class == MatchClass.SUPPORTED:
        support = (
            f"overlap={overlap:.4f}; {method_detail}; "
            f"a={row_a.strictness_class}; b={row_b.strictness_class}"
        )
        contradiction = "NA"
    elif match_class == MatchClass.PARTIALLY_SUPPORTED:
        support = f"overlap={overlap:.4f}; {method_detail}"
        contradiction = f"state_compat={state_compat}"
    elif match_class == MatchClass.CONTRADICTED:
        support = "NA"
        contradiction = f"a={row_a.strictness_class}; b={row_b.strictness_class}"
    else:
        support = "NA"
        contradiction = "NA"

    return {
        "compare_id": compare_id,
        "locus_id_a": row_a.locus_id,
        "locus_id_b": row_b.locus_id,
        "source_a": source_label_a,
        "source_b": source_label_b,
        "contig": row_a.contig,
        "start_a": str(row_a.start),
        "end_a": str(row_a.end),
        "start_b": str(row_b.start),
        "end_b": str(row_b.end),
        "match_class": match_class.value,
        "coordinate_overlap": f"{overlap:.4f}",
        "state_compatibility": str(state_compat),
        "strictness_a": row_a.strictness_class,
        "strictness_b": row_b.strictness_class,
        "support_summary": support,
        "contradiction_summary": contradiction,
        "comparison_score": str(score),
    }


def _source_specific_row(
    compare_id: str,
    row: HitsRow,
    source_label: str,
    *,
    is_a: bool,
) -> dict[str, str]:
    score = str(_MATCH_SCORES[MatchClass.SOURCE_SPECIFIC])
    if is_a:
        return {
            "compare_id": compare_id,
            "locus_id_a": row.locus_id,
            "locus_id_b": "NA",
            "source_a": source_label,
            "source_b": "NA",
            "contig": row.contig,
            "start_a": str(row.start),
            "end_a": str(row.end),
            "start_b": "NA",
            "end_b": "NA",
            "match_class": MatchClass.SOURCE_SPECIFIC.value,
            "coordinate_overlap": "0.0000",
            "state_compatibility": "False",
            "strictness_a": row.strictness_class,
            "strictness_b": "NA",
            "support_summary": "NA",
            "contradiction_summary": "NA",
            "comparison_score": score,
        }
    return {
        "compare_id": compare_id,
        "locus_id_a": "NA",
        "locus_id_b": row.locus_id,
        "source_a": "NA",
        "source_b": source_label,
        "contig": row.contig,
        "start_a": "NA",
        "end_a": "NA",
        "start_b": str(row.start),
        "end_b": str(row.end),
        "match_class": MatchClass.SOURCE_SPECIFIC.value,
        "coordinate_overlap": "0.0000",
        "state_compatibility": "False",
        "strictness_a": "NA",
        "strictness_b": row.strictness_class,
        "support_summary": "NA",
        "contradiction_summary": "NA",
        "comparison_score": score,
    }


def _write_compare_summary_tsv(
    compare_rows: list[dict[str, str]],
    outdir: Path,
) -> None:
    counts: Counter[str] = Counter(r["match_class"] for r in compare_rows)
    total = len(compare_rows)
    summary_rows = []
    for mc in MatchClass:
        n = counts.get(mc.value, 0)
        pct = round(100.0 * n / total, 1) if total > 0 else 0.0
        scores = [
            float(r["comparison_score"]) for r in compare_rows if r["match_class"] == mc.value
        ]
        overlaps = [
            float(r["coordinate_overlap"]) for r in compare_rows if r["match_class"] == mc.value
        ]
        summary_rows.append(
            {
                "match_class": mc.value,
                "n_loci": str(n),
                "pct_total": f"{pct:.1f}",
                "mean_overlap": f"{sum(overlaps) / len(overlaps):.4f}" if overlaps else "NA",
                "mean_score": f"{sum(scores) / len(scores):.4f}" if scores else "NA",
            }
        )
    with TsvWriter(outdir / "compare_summary.tsv", COMPARE_SUMMARY_COLUMNS) as w:
        w.write_rows(summary_rows)


def _sample_values(values: set[str], *, limit: int = 10) -> list[str]:
    return sorted(values)[:limit]


def _build_compare_diagnostics(
    rows_a: list[HitsRow],
    rows_b: list[HitsRow],
    compare_rows: list[dict[str, str]],
    cfg: CompareConfig,
    search_stats: MatchSearchStats,
) -> dict[str, object]:
    raw_contigs_a = {r.contig for r in rows_a}
    raw_contigs_b = {r.contig for r in rows_b}
    canonical_contigs_a = {
        canonicalize_contig(r.contig, normalize=cfg.normalize_contigs) for r in rows_a
    }
    canonical_contigs_b = {
        canonicalize_contig(r.contig, normalize=cfg.normalize_contigs) for r in rows_b
    }
    counts = Counter(r["match_class"] for r in compare_rows)
    matched_classes = {
        MatchClass.SUPPORTED.value,
        MatchClass.PARTIALLY_SUPPORTED.value,
        MatchClass.CONTRADICTED.value,
    }
    matched_pairs = sum(counts.get(mc, 0) for mc in matched_classes)

    return {
        "normalize_contigs": cfg.normalize_contigs,
        "overlap_mode": cfg.overlap_mode,
        "raw_contigs_a": len(raw_contigs_a),
        "raw_contigs_b": len(raw_contigs_b),
        "raw_shared_contigs": len(raw_contigs_a & raw_contigs_b),
        "canonical_contigs_a": len(canonical_contigs_a),
        "canonical_contigs_b": len(canonical_contigs_b),
        "canonical_shared_contigs": len(canonical_contigs_a & canonical_contigs_b),
        "raw_only_contigs_a_sample": _sample_values(raw_contigs_a - raw_contigs_b),
        "raw_only_contigs_b_sample": _sample_values(raw_contigs_b - raw_contigs_a),
        "canonical_only_contigs_a_sample": _sample_values(
            canonical_contigs_a - canonical_contigs_b
        ),
        "canonical_only_contigs_b_sample": _sample_values(
            canonical_contigs_b - canonical_contigs_a
        ),
        "rows_a_with_candidate_contig": search_stats.rows_a_with_candidate_contig,
        "rows_a_without_candidate_contig": search_stats.rows_a_without_candidate_contig,
        "candidate_pairs_checked": search_stats.candidate_pairs_checked,
        "overlapping_candidate_pairs": search_stats.overlapping_candidate_pairs,
        "qualifying_overlap_pairs": search_stats.qualifying_overlap_pairs,
        "overlap_matches": search_stats.overlap_matches,
        "breakpoint_candidate_pairs": search_stats.breakpoint_candidate_pairs,
        "breakpoint_fallback_matches": search_stats.breakpoint_fallback_matches,
        "matched_pairs": matched_pairs,
        "match_class_counts": dict(counts),
    }


def run_compare(
    hits_a: Path,
    hits_b: Path,
    outdir: Path,
    cfg: PrivyConfig,
    source_label_a: str | None = None,
    source_label_b: str | None = None,
    write_compare_tsv: bool = True,
    write_summary_tsv: bool = True,
    write_json: bool = True,
) -> list[dict[str, str]]:
    """Compare two hits.tsv files and write compare outputs.

    Args:
        hits_a: Path to the first (primary) hits.tsv file.
        hits_b: Path to the second (comparison) hits.tsv file.
        outdir: Output directory (created if it does not exist).
        cfg: Resolved :class:`~privy.core.config.PrivyConfig`; the
            ``compare`` sub-model controls matching thresholds.
        source_label_a: Display label for source A.  Inferred from the
            locus_id prefix if ``None``.
        source_label_b: Display label for source B.
        write_compare_tsv: Write ``compare.tsv``.
        write_summary_tsv: Write ``compare_summary.tsv``.
        write_json: Write ``compare.json``.

    Returns:
        List of compare row dicts (same content as compare.tsv).
    """
    rows_a = load_hits_tsv(hits_a)
    rows_b = load_hits_tsv(hits_b)
    log.info("Loaded %d rows from source A (%s)", len(rows_a), hits_a.name)
    log.info("Loaded %d rows from source B (%s)", len(rows_b), hits_b.name)

    label_a = infer_source_label(rows_a, source_label_a)
    label_b = infer_source_label(rows_b, source_label_b)
    compare_cfg = cfg.compare

    index_b = _build_contig_index(rows_b, compare_cfg)
    search_stats = MatchSearchStats()
    matched_b_ids: set[str] = set()
    compare_rows: list[dict[str, str]] = []
    counter = 0

    for row_a in rows_a:
        contig_key = canonicalize_contig(
            row_a.contig,
            normalize=compare_cfg.normalize_contigs,
        )
        contig_index = index_b.get(contig_key)
        if contig_index is not None:
            search_stats.rows_a_with_candidate_contig += 1
            candidates = _candidate_window(row_a, contig_index, compare_cfg)
        else:
            search_stats.rows_a_without_candidate_contig += 1
            candidates = []
        matches = find_match_results(row_a, candidates, compare_cfg, search_stats)

        if not matches:
            counter += 1
            compare_id = f"CMP{counter:06d}"
            compare_rows.append(_source_specific_row(compare_id, row_a, label_a, is_a=True))
            continue

        for match in matches:
            counter += 1
            compare_id = f"CMP{counter:06d}"
            state_compat = is_state_compatible(
                row_a.strictness_class,
                match.row.strictness_class,
                compare_cfg.require_state_compatibility,
            )
            mc = classify_match(
                match.overlap,
                state_compat,
                row_a.strictness_class,
                match.row.strictness_class,
                compare_cfg,
                method=match.method,
            )
            if mc != MatchClass.SOURCE_SPECIFIC:
                matched_b_ids.add(match.row.locus_id)
            score = compute_comparison_score(mc, match.overlap)
            compare_rows.append(
                _pair_to_row(
                    compare_id,
                    row_a,
                    match.row,
                    mc,
                    match.overlap,
                    state_compat,
                    score,
                    label_a,
                    label_b,
                    method=match.method,
                    breakpoint_distance_bp=match.breakpoint_distance_bp,
                )
            )

    for row_b in rows_b:
        if row_b.locus_id in matched_b_ids:
            continue
        counter += 1
        compare_id = f"CMP{counter:06d}"
        compare_rows.append(_source_specific_row(compare_id, row_b, label_b, is_a=False))

    outdir.mkdir(parents=True, exist_ok=True)

    if write_compare_tsv:
        with TsvWriter(outdir / "compare.tsv", COMPARE_COLUMNS) as w:
            w.write_rows(compare_rows)
        log.info("Wrote compare.tsv (%d rows)", len(compare_rows))

    if write_summary_tsv:
        _write_compare_summary_tsv(compare_rows, outdir)
        log.info("Wrote compare_summary.tsv")

    if write_json:
        diagnostics = _build_compare_diagnostics(
            rows_a,
            rows_b,
            compare_rows,
            compare_cfg,
            search_stats,
        )
        meta: dict[str, object] = {
            "tool": "privy compare",
            "hits_a": str(hits_a),
            "hits_b": str(hits_b),
            "source_a": label_a,
            "source_b": label_b,
            "n_rows_a": len(rows_a),
            "n_rows_b": len(rows_b),
            "n_compare_rows": len(compare_rows),
            "config": cfg.compare.model_dump(),
            "diagnostics": diagnostics,
            "timestamp": datetime.now(timezone.utc).isoformat(),
        }
        (outdir / "compare.json").write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
        log.info("Wrote compare.json")

    return compare_rows
