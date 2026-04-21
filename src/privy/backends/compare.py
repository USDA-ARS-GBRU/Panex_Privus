"""Cross-source comparison engine for ``privy compare``.

Compares two ``hits.tsv`` files produced by separate ``privy scan`` runs
(typically one VCF scan and one GFA scan) and classifies each locus pair
as supported, partially_supported, contradicted, source_specific,
uninformative, or missing_data.

Algorithm:
    1. Load source-A and source-B hits rows.
    2. Build a contig-bucketed index over source B for fast lookup.
    3. For each source-A locus, find the best overlapping source-B locus
       (reciprocal overlap ≥ ``min_reciprocal_overlap``, or within
       ``breakpoint_tolerance_bp`` as a fallback).
    4. Classify the pair: SUPPORTED / PARTIALLY_SUPPORTED / CONTRADICTED.
    5. Emit unmatched source-B loci as SOURCE_SPECIFIC entries.
    6. Write compare.tsv, compare_summary.tsv, compare.json.
"""

from __future__ import annotations

import json
import logging
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

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

_STRICT_CLASSES = frozenset({
    "strict_complete",
    "strict_target_missing",
    "strict_offtarget_missing",
    "strict_both_missing",
})


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
            result.append(HitsRow(
                locus_id=raw["locus_id"],
                contig=raw["contig"],
                start=int(raw["start"]),
                end=int(raw["end"]),
                variant_type=raw.get("variant_type", ""),
                allele_key=raw.get("allele_key", ""),
                strictness_class=raw.get("strictness_class", ""),
                final_score=float(raw.get("final_score", 0.0)),
                raw=raw,
            ))
        except (KeyError, ValueError) as exc:
            raise ValueError(
                f"Malformed hits.tsv row in {path}: {raw!r}"
            ) from exc
    return result


def infer_source_label(rows: list[HitsRow], explicit: Optional[str]) -> str:
    """Return the display label for an evidence source.

    Uses *explicit* when provided.  Otherwise infers from locus_id prefix:
    ``GPX`` → ``"gfa"``, anything else → ``"vcf"``.
    """
    if explicit:
        return explicit
    if rows and rows[0].locus_id.startswith("GPX"):
        return "gfa"
    return "vcf"


def reciprocal_overlap_rows(a: HitsRow, b: HitsRow) -> float:
    """Reciprocal overlap (intersection / union) between two HitsRow objects.

    Returns 0.0 if loci are on different contigs or do not overlap.
    """
    if a.contig != b.contig:
        return 0.0
    ov_start = max(a.start, b.start)
    ov_end = min(a.end, b.end)
    if ov_end <= ov_start:
        return 0.0
    intersection = ov_end - ov_start
    union = max(a.end, b.end) - min(a.start, b.start)
    return intersection / union if union > 0 else 0.0


def _breakpoint_distance(a: HitsRow, b: HitsRow) -> int:
    """Minimum bp gap between two loci on the same contig."""
    if a.contig != b.contig:
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
) -> MatchClass:
    """Determine the MatchClass for a locus pair.

    Precedence:
        1. Either strictness is ``contradicted`` → CONTRADICTED
        2. No coordinate overlap → SOURCE_SPECIFIC
        3. ``require_state_compatibility`` and states mismatch → PARTIALLY_SUPPORTED
        4. States compatible → SUPPORTED
        5. Otherwise → PARTIALLY_SUPPORTED
    """
    if strictness_a == "contradicted" or strictness_b == "contradicted":
        return MatchClass.CONTRADICTED
    if overlap <= 0.0:
        return MatchClass.SOURCE_SPECIFIC
    if cfg.require_state_compatibility and not state_compat:
        return MatchClass.PARTIALLY_SUPPORTED
    if state_compat:
        return MatchClass.SUPPORTED
    return MatchClass.PARTIALLY_SUPPORTED


def compute_comparison_score(match_class: MatchClass, overlap: float) -> float:
    """Numeric comparison score in [0, 1].

    For SUPPORTED and PARTIALLY_SUPPORTED the base score is scaled by
    the reciprocal overlap so higher-confidence matches score higher.
    """
    base = _MATCH_SCORES[match_class]
    if match_class in (MatchClass.SUPPORTED, MatchClass.PARTIALLY_SUPPORTED):
        return round(base * max(overlap, 0.01), 4)
    return base


def _build_contig_index(rows: list[HitsRow]) -> dict[str, list[HitsRow]]:
    idx: dict[str, list[HitsRow]] = defaultdict(list)
    for r in rows:
        idx[r.contig].append(r)
    return dict(idx)


def find_best_match(
    query: HitsRow,
    candidates: list[HitsRow],
    cfg: CompareConfig,
) -> tuple[Optional[HitsRow], float]:
    """Return ``(best_candidate, overlap)`` for a query locus.

    Primary search: reciprocal overlap ≥ ``cfg.min_reciprocal_overlap``.
    Fallback: gap ≤ ``cfg.breakpoint_tolerance_bp`` (returns the nearest
    candidate even if overlap is zero).

    Returns ``(None, 0.0)`` when no candidate qualifies.
    """
    best: Optional[HitsRow] = None
    best_overlap = 0.0

    for cand in candidates:
        ov = reciprocal_overlap_rows(query, cand)
        if ov >= cfg.min_reciprocal_overlap and ov > best_overlap:
            best = cand
            best_overlap = ov

    if best is None and cfg.breakpoint_tolerance_bp > 0:
        best_dist = cfg.breakpoint_tolerance_bp + 1
        for cand in candidates:
            dist = _breakpoint_distance(query, cand)
            if dist <= cfg.breakpoint_tolerance_bp and dist < best_dist:
                best = cand
                best_dist = dist
                best_overlap = reciprocal_overlap_rows(query, cand)

    return best, best_overlap


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
) -> dict[str, str]:
    if match_class == MatchClass.SUPPORTED:
        support = f"overlap={overlap:.4f}; a={row_a.strictness_class}; b={row_b.strictness_class}"
        contradiction = "NA"
    elif match_class == MatchClass.PARTIALLY_SUPPORTED:
        support = f"overlap={overlap:.4f}"
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
        scores = [float(r["comparison_score"]) for r in compare_rows if r["match_class"] == mc.value]
        overlaps = [float(r["coordinate_overlap"]) for r in compare_rows if r["match_class"] == mc.value]
        summary_rows.append({
            "match_class": mc.value,
            "n_loci": str(n),
            "pct_total": f"{pct:.1f}",
            "mean_overlap": f"{sum(overlaps) / len(overlaps):.4f}" if overlaps else "NA",
            "mean_score": f"{sum(scores) / len(scores):.4f}" if scores else "NA",
        })
    with TsvWriter(outdir / "compare_summary.tsv", COMPARE_SUMMARY_COLUMNS) as w:
        w.write_rows(summary_rows)


def run_compare(
    hits_a: Path,
    hits_b: Path,
    outdir: Path,
    cfg: PrivyConfig,
    source_label_a: Optional[str] = None,
    source_label_b: Optional[str] = None,
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

    index_b = _build_contig_index(rows_b)
    matched_b_ids: set[str] = set()
    compare_rows: list[dict[str, str]] = []
    counter = 0

    for row_a in rows_a:
        counter += 1
        compare_id = f"CMP{counter:06d}"
        candidates = index_b.get(row_a.contig, [])
        best_b, overlap = find_best_match(row_a, candidates, compare_cfg)

        if best_b is None:
            compare_rows.append(
                _source_specific_row(compare_id, row_a, label_a, is_a=True)
            )
            continue

        matched_b_ids.add(best_b.locus_id)
        state_compat = is_state_compatible(
            row_a.strictness_class,
            best_b.strictness_class,
            compare_cfg.require_state_compatibility,
        )
        mc = classify_match(
            overlap, state_compat,
            row_a.strictness_class, best_b.strictness_class,
            compare_cfg,
        )
        score = compute_comparison_score(mc, overlap)
        compare_rows.append(
            _pair_to_row(
                compare_id, row_a, best_b, mc, overlap,
                state_compat, score, label_a, label_b,
            )
        )

    for row_b in rows_b:
        if row_b.locus_id in matched_b_ids:
            continue
        counter += 1
        compare_id = f"CMP{counter:06d}"
        compare_rows.append(
            _source_specific_row(compare_id, row_b, label_b, is_a=False)
        )

    outdir.mkdir(parents=True, exist_ok=True)

    if write_compare_tsv:
        with TsvWriter(outdir / "compare.tsv", COMPARE_COLUMNS) as w:
            w.write_rows(compare_rows)
        log.info("Wrote compare.tsv (%d rows)", len(compare_rows))

    if write_summary_tsv:
        _write_compare_summary_tsv(compare_rows, outdir)
        log.info("Wrote compare_summary.tsv")

    if write_json:
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
            "timestamp": datetime.now(timezone.utc).isoformat(),
        }
        (outdir / "compare.json").write_text(
            json.dumps(meta, indent=2) + "\n", encoding="utf-8"
        )
        log.info("Wrote compare.json")

    return compare_rows
