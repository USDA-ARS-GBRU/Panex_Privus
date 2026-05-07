"""Standalone GFA private-segment discovery backend.

Discovers target-private segments from a pangenome graph, using the same
StrictnessClass framework and output schema as the VCF backend.  When the
user provides only ``--gfa`` (no ``--vcf``), this backend is the primary
discovery engine.

Discovery model
---------------
A graph segment is "target-private" when:

- Target-sample paths/walks traverse it, and
- Off-target paths/walks do not traverse it (or traverse a different segment
  at the same locus — an alternative bubble arm).

Missingness is defined the same way as in VCF analysis: a sample is *missing*
for a segment's locus if it has no path or walk coverage at that genomic
position at all (not merely traversing an alternative arm).

Coordinate requirement
-----------------------
Segments must carry SN/SO/LN optional tags (minigraph/PGGB output convention)
to be included in the scan.  Segments without coordinates cannot be placed on
the output grid and are skipped with a warning.  W-line coordinate ranges are
used to detect sample presence/absence at each locus.

Output schema
-------------
Produces the same files as :mod:`~privy.backends.vcf_scan`:
hits.tsv, regions.tsv, evidence.tsv, sample_support.tsv, qc.tsv, run.json.
All locus IDs use the ``GPX`` prefix (Graph Private region X).
"""

from __future__ import annotations

import logging
import time
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

from privy.core.cohort import CohortDefinition
from privy.core.config import PrivyConfig
from privy.core.intervals import merge_loci_to_regions
from privy.core.locus import Locus, LocusType, PrimarySource
from privy.core.patterns import AllelePattern, build_allele_pattern
from privy.core.scoring import (
    ScoredHit,
    compute_discovery_score,
    compute_final_score,
    compute_penalty_score,
    rank_scored_hits,
)
from privy.io.gfa import (
    GfaScanIndex,
    build_gfa_scan_index,
    load_gfa_scan_index,
)
from privy.io.jsonio import write_run_json as _write_json
from privy.io.tsv import (
    EVIDENCE_COLUMNS,
    GFA_SEGMENT_COLUMNS,
    HITS_COLUMNS,
    QC_COLUMNS,
    REGIONS_COLUMNS,
    SAMPLE_SUPPORT_COLUMNS,
    TsvWriter,
)
from privy.utils.metrics import ScanStats
from privy.utils.misc import now_iso

log = logging.getLogger("privy.backends.gfa_scan")


# ---------------------------------------------------------------------------
# Internal data structures
# ---------------------------------------------------------------------------


@dataclass
class GfaHitRecord:
    """Internal record for one passing segment, produced during the scan loop.

    Mirrors :class:`~privy.backends.vcf_scan.HitRecord` but carries
    GFA-specific fields (segment name, traversal map) instead of VCF fields.
    """

    locus_id: str
    locus: Locus
    allele_key: str          # "contig:start+1:REF:seg_name"  (1-based for display)
    segment_name: str
    pattern: AllelePattern
    # sample_name → "traverses" | "absent" | "missing"
    sample_traversal: dict[str, str]
    discovery_score: float = 0.0
    support_score: float = 0.0
    penalty_score: float = 0.0
    final_score: float = 0.0
    rank: int = 0


@dataclass
class _PresenceTracker:
    """Mutable cursor for one sample's sorted presence intervals on a contig."""

    bit: int
    starts: tuple[int, ...]
    ends: tuple[int, ...]
    pos: int = 0


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------


def run_gfa_scan(
    gfa: Path,
    cohort: CohortDefinition,
    cfg: PrivyConfig,
    outdir: Path,
    mode: str = "private_allele",
    region: str | None = None,
    contig: str | None = None,
    write_hits: bool = True,
    write_regions: bool = True,
    write_evidence: bool = True,
    write_sample_support: bool = True,
    write_qc: bool = True,
    write_run_json: bool = True,
    threads: int = 1,
    gfa_index: Path | None = None,
) -> None:
    """Run the standalone GFA private-segment scan.

    Args:
        gfa: Path to a GFA1/1.1 file.
        cohort: Validated cohort definition.
        cfg: Resolved configuration (``cfg.gfa`` controls GFA-specific options,
             ``cfg.scan`` controls discovery thresholds).
        outdir: Output directory (must already exist).
        mode: Discovery mode.  Only ``"private_allele"`` is implemented.
        region: Optional region filter ``"contig:start-end"`` (1-based samtools).
        contig: Optional single-contig restriction.
        write_hits: Write ``hits.tsv``.
        write_regions: Write ``regions.tsv``.
        write_evidence: Write ``evidence.tsv``.
        write_sample_support: Write ``sample_support.tsv``.
        write_qc: Write ``qc.tsv``.
        write_run_json: Write ``run.json``.
        threads: Worker threads. Values greater than 1 currently run serially.
        gfa_index: Optional prebuilt Privy GFA scan index.

    Raises:
        FileNotFoundError: If *gfa* does not exist.
        ValueError: If no cohort samples are found in the GFA.
        NotImplementedError: If *mode* is not ``"private_allele"``.
    """
    start_time = now_iso()
    if threads > 1:
        log.warning("GFA scan parallel execution is not implemented; running serially.")

    if mode != "private_allele":
        raise NotImplementedError(
            f"Discovery mode {mode!r} is not yet implemented for GFA. "
            "Only 'private_allele' is available."
        )

    if not gfa.exists():
        raise FileNotFoundError(f"GFA file not found: {gfa}")

    # ── Step 1: Parse region / contig filter ─────────────────────────────────
    filter_contig: str | None = None
    filter_start: int | None = None
    filter_end: int | None = None

    if region is not None:
        filter_contig, filter_start, filter_end = _parse_region(region)
    elif contig is not None:
        filter_contig = contig

    all_targets = list(cohort.targets)
    all_offtargets = list(cohort.off_targets)
    if gfa_index is not None:
        log.info("Loading GFA scan index: %s", gfa_index)
        scan_index = load_gfa_scan_index(gfa_index, gfa_path=gfa)
    else:
        log.info("Indexing GFA for cohort scan: %s", gfa)
        scan_index = build_gfa_scan_index(
            gfa_path=gfa,
            sample_names=all_targets + all_offtargets,
            filter_contig=filter_contig,
            filter_start=filter_start,
            filter_end=filter_end,
        )

    # ── Step 3: Validate samples ──────────────────────────────────────────────
    gfa_samples = scan_index.samples_seen

    active_targets = [s for s in cohort.targets if s in gfa_samples]
    active_offtargets = [s for s in cohort.off_targets if s in gfa_samples]
    missing_targets = [s for s in cohort.targets if s not in gfa_samples]
    missing_offtargets = [s for s in cohort.off_targets if s not in gfa_samples]

    if missing_targets:
        log.warning(
            "Target samples not found in GFA (will count as missing): %s",
            ", ".join(missing_targets),
        )
    if missing_offtargets:
        log.warning(
            "Off-target samples not found in GFA (will count as missing): %s",
            ", ".join(missing_offtargets),
        )
    if not active_targets:
        raise ValueError(
            "No target samples from the cohort definition were found in the GFA. "
            f"Cohort targets: {list(cohort.targets)}. "
            f"GFA samples: {sorted(gfa_samples)!r}."
        )
    if not active_offtargets:
        raise ValueError(
            "No off-target samples from the cohort definition were found in the GFA. "
            f"Cohort off-targets: {list(cohort.off_targets)}."
        )

    log.info(
        "Cohort validated | targets=%d/%d in GFA | off-targets=%d/%d in GFA",
        len(active_targets), cohort.n_targets,
        len(active_offtargets), cohort.n_off_targets,
    )
    log.info(
        "Starting GFA scan | contigs=%d | coordinate_segments=%d | mode=%s | "
        "targets=%d | off-targets=%d",
        len(scan_index.contig_names()),
        scan_index.coordinate_segment_count(),
        mode,
        len(active_targets),
        len(active_offtargets),
    )

    # ── Step 4: Scan segments ─────────────────────────────────────────────────
    stats = ScanStats(
        n_target_samples=len(active_targets),
        n_offtarget_samples=len(active_offtargets),
    )

    hit_records = _scan_segments(
        scan_index=scan_index,
        all_targets=all_targets,
        all_offtargets=all_offtargets,
        cfg=cfg,
        stats=stats,
        filter_contig=filter_contig,
        filter_start=filter_start,
        filter_end=filter_end,
    )

    stats.loci_emitted = len(hit_records)
    log.info(
        "Scan complete | hits=%d | contradicted=%d | segments_evaluated=%d",
        len(hit_records), stats.alleles_contradicted, stats.alleles_evaluated,
    )

    # ── Step 5: Score hits ────────────────────────────────────────────────────
    _score_hit_records(hit_records, cfg)

    # ── Step 6: Merge to regions ──────────────────────────────────────────────
    hit_loci = [hr.locus for hr in hit_records]
    hit_by_id: dict[str, GfaHitRecord] = {hr.locus_id: hr for hr in hit_records}

    region_loci: list[Locus] = []
    if cfg.scan.merge_distance >= 0:
        region_loci = merge_loci_to_regions(
            hit_loci,
            merge_distance=cfg.scan.merge_distance,
            same_variant_class_only=cfg.scan.same_variant_class_only,
            region_id_prefix="GFAR",
        )

    stats.regions_emitted = len(region_loci)
    log.info(
        "Region merging complete | regions=%d (merge_distance=%d bp)",
        len(region_loci), cfg.scan.merge_distance,
    )

    # ── Step 7: Write outputs ─────────────────────────────────────────────────
    if write_hits:
        _write_hits_tsv(hit_records, outdir)
        _write_graph_segments_tsv(hit_records, outdir)
    if write_regions:
        _write_regions_tsv(region_loci, hit_by_id, outdir)
    if write_evidence:
        _write_evidence_tsv(hit_records, outdir)
    if write_sample_support:
        _write_sample_support_tsv(hit_records, cohort, outdir)
    if write_qc:
        _write_qc_tsv(stats, outdir)
    if write_run_json:
        _write_run_json_file(
            cfg=cfg,
            cohort=cohort,
            stats=stats,
            gfa_path=gfa,
            gfa_index_path=gfa_index,
            outdir=outdir,
            start_time=start_time,
            end_time=now_iso(),
        )

    log.info("Outputs written to %s", outdir)


# ---------------------------------------------------------------------------
# Scan loop
# ---------------------------------------------------------------------------


def _scan_segments(
    scan_index: GfaScanIndex,
    all_targets: list[str],
    all_offtargets: list[str],
    cfg: PrivyConfig,
    stats: ScanStats,
    filter_contig: str | None,
    filter_start: int | None,
    filter_end: int | None,
) -> list[GfaHitRecord]:
    """Evaluate every coordinate-tagged segment and return passing hits."""
    hits: list[GfaHitRecord] = []
    locus_n = 0
    skipped_too_short = 0
    min_len = cfg.gfa.min_segment_length
    target_mask = scan_index.sample_mask(all_targets)
    offtarget_mask = scan_index.sample_mask(all_offtargets)

    # Determine contigs to scan
    if filter_contig is not None:
        contigs = list(scan_index.matching_contigs(filter_contig))
    else:
        contigs = sorted(scan_index.contig_names())

    contigs_visited: set[str] = set()

    for ctg in contigs:
        if not scan_index.has_contig(ctg):
            continue
        contigs_visited.add(ctg)
        contig_hits_before = len(hits)
        contig_evaluated_before = stats.alleles_evaluated
        contig_rows_seen = 0
        contig_total_rows = scan_index.contig_segment_count(ctg)
        progress_next = 1_000_000
        last_time_progress = time.monotonic()
        presence_trackers = _build_presence_trackers(scan_index, ctg)

        for seg_start, seg_end, seg_name, support_mask in scan_index.iter_contig_segments(ctg):
            contig_rows_seen += 1
            now = time.monotonic()
            if now - last_time_progress >= 30.0:
                log.info(
                    "  %s: scanned %d/%d coordinate segments, evaluated %d "
                    "target-supported segments, %d hits so far",
                    ctg,
                    contig_rows_seen,
                    contig_total_rows,
                    stats.alleles_evaluated - contig_evaluated_before,
                    len(hits) - contig_hits_before,
                )
                last_time_progress = now

            # Apply region filter
            if filter_start is not None and seg_end <= filter_start:
                continue
            if filter_end is not None and seg_start >= filter_end:
                continue

            # Minimum segment length filter
            seg_len = seg_end - seg_start
            if seg_len < min_len:
                skipped_too_short += 1
                continue

            if not (support_mask & target_mask):
                continue

            stats.records_evaluated += 1
            stats.alleles_evaluated += 1

            present_mask = (
                _present_mask_for_sorted_locus(presence_trackers, seg_start, seg_end)
                | support_mask
            )
            ts_n = (support_mask & target_mask).bit_count()
            os_n = (support_mask & offtarget_mask).bit_count()
            tm_n = (target_mask & ~present_mask).bit_count()
            om_n = (offtarget_mask & ~present_mask).bit_count()
            tt_n = len(all_targets)
            ot_n = len(all_offtargets)

            allele_key = _format_segment_allele_key(ctg, seg_start, seg_name)

            pattern = build_allele_pattern(
                allele_key=allele_key,
                target_support_n=ts_n,
                target_total_n=tt_n,
                offtarget_support_n=os_n,
                offtarget_total_n=ot_n,
                target_missing_n=tm_n,
                offtarget_missing_n=om_n,
                min_target_support=cfg.scan.min_target_support,
                max_offtarget_support=cfg.scan.max_off_target_support,
                relaxed_target_missing=cfg.scan.relaxed_target_missing,
                relaxed_offtarget_missing=cfg.scan.relaxed_offtarget_missing,
            )

            stats.increment_strictness(pattern.strictness_class.value)
            if pattern.strictness_class.value == "contradicted":
                stats.alleles_contradicted += 1

            contig_evaluated = stats.alleles_evaluated - contig_evaluated_before
            if contig_evaluated >= progress_next:
                log.info(
                    "  %s: evaluated %d segments, %d hits so far",
                    ctg,
                    contig_evaluated,
                    len(hits) - contig_hits_before,
                )
                progress_next += 1_000_000

            if not pattern.pattern_pass:
                continue

            stats.alleles_passed += 1
            locus_n += 1
            locus_id = f"GPX{locus_n:08d}"

            locus = Locus(
                locus_id=locus_id,
                contig=ctg,
                start=seg_start,
                end=seg_end,
                locus_type=LocusType.GRAPH_REGION,
                primary_source=PrimarySource.GFA,
                source_ids=[locus_id],
                metadata={"segment": seg_name},
            )

            sample_traversal = scan_index.mask_to_statuses(
                support_mask=support_mask,
                present_mask=present_mask,
                samples=all_targets + all_offtargets,
            )

            hits.append(GfaHitRecord(
                locus_id=locus_id,
                locus=locus,
                allele_key=allele_key,
                segment_name=seg_name,
                pattern=pattern,
                sample_traversal=sample_traversal,
            ))

        log.info(
            "  %s: %d hits (%d segments evaluated)",
            ctg,
            len(hits) - contig_hits_before,
            stats.alleles_evaluated - contig_evaluated_before,
        )

    stats.n_contigs_scanned = len(contigs_visited)

    if skipped_too_short:
        log.debug(
            "Skipped %d segments shorter than min_segment_length=%d bp",
            skipped_too_short, min_len,
        )

    return hits


def _build_presence_trackers(
    scan_index: GfaScanIndex,
    contig: str,
) -> list[_PresenceTracker]:
    """Build per-sample interval cursors for one sorted contig scan."""
    trackers: list[_PresenceTracker] = []
    for sample_idx, intervals_by_contig in enumerate(scan_index.sample_intervals):
        intervals = intervals_by_contig.get(contig)
        if intervals is None:
            continue
        trackers.append(_PresenceTracker(
            bit=1 << sample_idx,
            starts=intervals.starts,
            ends=intervals.ends,
        ))
    return trackers


def _present_mask_for_sorted_locus(
    trackers: list[_PresenceTracker],
    start: int,
    end: int,
) -> int:
    """Return presence mask while advancing interval cursors for sorted loci."""
    mask = 0
    for tracker in trackers:
        while tracker.pos < len(tracker.ends) and tracker.ends[tracker.pos] <= start:
            tracker.pos += 1
        if tracker.pos < len(tracker.starts) and tracker.starts[tracker.pos] < end:
            mask |= tracker.bit
    return mask


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------


def _score_hit_records(hit_records: list[GfaHitRecord], cfg: PrivyConfig) -> None:
    """Score all hits in-place using the configured scoring weights."""
    scored: list[ScoredHit] = []

    for hr in hit_records:
        ds = compute_discovery_score(
            hr.pattern,
            variant_qual=None,
            discovery_weight=cfg.scoring.discovery_weight,
        )
        ps = compute_penalty_score(hr.pattern, penalty_weight=cfg.scoring.penalty_weight)
        ss = 0.0  # BAM/XMFA support not yet wired
        fs = compute_final_score(ds, ss, ps)

        hr.discovery_score = ds
        hr.support_score = ss
        hr.penalty_score = ps
        hr.final_score = fs

        scored.append(ScoredHit(
            locus_id=hr.locus_id,
            discovery_score=ds,
            support_score=ss,
            penalty_score=ps,
            final_score=fs,
            rank=0,
            strictness_class=hr.pattern.strictness_class.value,
            summary_label=hr.pattern.pattern_reason,
        ))

    ranked = rank_scored_hits(scored)
    rank_map = {s.locus_id: s.rank for s in ranked}
    for hr in hit_records:
        hr.rank = rank_map.get(hr.locus_id, 0)


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------


def _write_hits_tsv(hit_records: list[GfaHitRecord], outdir: Path) -> None:
    path = outdir / "hits.tsv"
    sorted_hits = sorted(hit_records, key=lambda h: h.rank)
    with TsvWriter(path, HITS_COLUMNS) as w:
        for hr in sorted_hits:
            w.write_row({
                "locus_id":            hr.locus_id,
                "contig":              hr.locus.contig,
                "start":               hr.locus.start,
                "end":                 hr.locus.end,
                "variant_type":        "graph_region",
                "allele_key":          hr.allele_key,
                "target_support_n":    hr.pattern.target_support_n,
                "target_total_n":      hr.pattern.target_total_n,
                "offtarget_support_n": hr.pattern.offtarget_support_n,
                "offtarget_total_n":   hr.pattern.offtarget_total_n,
                "target_missing_n":    hr.pattern.target_missing_n,
                "offtarget_missing_n": hr.pattern.offtarget_missing_n,
                "strictness_class":    hr.pattern.strictness_class.value,
                "discovery_score":     hr.discovery_score,
                "support_score":       hr.support_score,
                "penalty_score":       hr.penalty_score,
                "final_score":         hr.final_score,
            })
    log.info("Wrote %s (%d rows)", path, len(hit_records))


def _write_graph_segments_tsv(hit_records: list[GfaHitRecord], outdir: Path) -> None:
    """Write a GFA-specific companion table for private graph-node evidence."""
    path = outdir / "graph_segments.tsv"
    sorted_hits = sorted(hit_records, key=lambda h: h.rank)
    with TsvWriter(path, GFA_SEGMENT_COLUMNS) as w:
        for hr in sorted_hits:
            pattern = hr.pattern
            target_covered_n = pattern.target_total_n - pattern.target_missing_n
            offtarget_covered_n = pattern.offtarget_total_n - pattern.offtarget_missing_n
            offtarget_absent_n = (
                pattern.offtarget_total_n
                - pattern.offtarget_support_n
                - pattern.offtarget_missing_n
            )
            segment_length = hr.locus.length
            w.write_row({
                "locus_id": hr.locus_id,
                "contig": hr.locus.contig,
                "start": hr.locus.start,
                "end": hr.locus.end,
                "segment_name": hr.segment_name,
                "segment_length": segment_length,
                "segment_length_class": _segment_length_class(segment_length),
                "graph_signal_type": "target_traversed_graph_segment",
                "target_traverse_n": pattern.target_support_n,
                "target_total_n": pattern.target_total_n,
                "target_coordinate_covered_n": target_covered_n,
                "target_missing_n": pattern.target_missing_n,
                "offtarget_same_segment_traverse_n": pattern.offtarget_support_n,
                "offtarget_same_segment_absent_n": offtarget_absent_n,
                "offtarget_coordinate_covered_n": offtarget_covered_n,
                "offtarget_missing_n": pattern.offtarget_missing_n,
                "offtarget_total_n": pattern.offtarget_total_n,
                "strictness_class": pattern.strictness_class.value,
                "interpretation": _graph_segment_interpretation(pattern),
            })
    log.info("Wrote %s (%d rows)", path, len(hit_records))


def _write_regions_tsv(
    region_loci: list[Locus],
    hit_by_id: dict[str, GfaHitRecord],
    outdir: Path,
) -> None:
    path = outdir / "regions.tsv"
    with TsvWriter(path, REGIONS_COLUMNS) as w:
        for region in region_loci:
            constituent = [hit_by_id[sid] for sid in region.source_ids if sid in hit_by_id]
            if not constituent:
                continue

            strictness_counter: Counter[str] = Counter(
                c.pattern.strictness_class.value for c in constituent
            )
            dominant = strictness_counter.most_common(1)[0][0]

            n_target_complete = sum(1 for c in constituent if c.pattern.target_missing_n == 0)
            n_offtarget_excluded = sum(
                1 for c in constituent if c.pattern.offtarget_support_n == 0
            )
            mean_score = sum(c.final_score for c in constituent) / len(constituent)

            w.write_row({
                "region_id":               region.locus_id,
                "contig":                  region.contig,
                "start":                   region.start,
                "end":                     region.end,
                "n_loci":                  len(constituent),
                "variant_types":           "graph_region",
                "dominant_strictness_class": dominant,
                "target_consistency":      round(n_target_complete / len(constituent), 4),
                "offtarget_exclusion":     round(n_offtarget_excluded / len(constituent), 4),
                "final_score":             round(mean_score, 6),
            })
    log.info("Wrote %s (%d rows)", path, len(region_loci))


def _write_evidence_tsv(hit_records: list[GfaHitRecord], outdir: Path) -> None:
    path = outdir / "evidence.tsv"
    with TsvWriter(path, EVIDENCE_COLUMNS) as w:
        for hr in hit_records:
            w.write_row({
                "locus_id":      hr.locus_id,
                "source_type":   "gfa",
                "sample_id":     "",
                "evidence_class": "support" if hr.pattern.pattern_pass else "contradiction",
                "metric_name":   "segment_pattern",
                "metric_value":  hr.final_score,
                "details":       hr.pattern.pattern_reason,
            })
    log.info("Wrote %s (%d rows)", path, len(hit_records))


def _write_sample_support_tsv(
    hit_records: list[GfaHitRecord],
    cohort: CohortDefinition,
    outdir: Path,
) -> None:
    path = outdir / "sample_support.tsv"
    with TsvWriter(path, SAMPLE_SUPPORT_COLUMNS) as w:
        for hr in sorted(hit_records, key=lambda h: h.rank):
            for sample, status in hr.sample_traversal.items():
                if cohort.is_ignored(sample):
                    continue
                cohort_role = "target" if cohort.is_target(sample) else "off_target"
                allele_supported = str(status == "traverses").lower()
                evidence_class = _traversal_to_evidence_class(status, cohort_role)
                w.write_row({
                    "locus_id":        hr.locus_id,
                    "sample_id":       sample,
                    "cohort_role":     cohort_role,
                    "genotype":        status,    # traversal status stands in for genotype
                    "allele_supported": allele_supported,
                    "depth":           "NA",
                    "allele_fraction": "NA",
                    "evidence_class":  evidence_class,
                })
    n_rows = len(hit_records) * len(hit_records[0].sample_traversal) if hit_records else 0
    log.info("Wrote %s (~%d rows)", path, n_rows)


def _write_qc_tsv(stats: ScanStats, outdir: Path) -> None:
    path = outdir / "qc.tsv"
    with TsvWriter(path, QC_COLUMNS) as w:
        w.write_rows(stats.as_qc_rows(source="gfa"))
    log.info("Wrote %s", path)


def _write_run_json_file(
    cfg: PrivyConfig,
    cohort: CohortDefinition,
    stats: ScanStats,
    gfa_path: Path,
    gfa_index_path: Path | None,
    outdir: Path,
    start_time: str,
    end_time: str,
) -> None:
    from privy import __version__  # noqa: PLC0415

    run_cfg = cfg.model_copy(update={
        "gfa": cfg.gfa.model_copy(update={"enabled": True}),
    })
    data = {
        "privy_version": __version__,
        "start_time": start_time,
        "end_time": end_time,
        "project_name": cfg.project_name,
        "config": run_cfg.as_run_dict(),
        "cohort": {
            "targets": list(cohort.targets),
            "off_targets": list(cohort.off_targets),
            "ignored_samples": list(cohort.ignored_samples),
            "n_targets": cohort.n_targets,
            "n_off_targets": cohort.n_off_targets,
        },
        "inputs": {
            "gfa": str(gfa_path),
            "gfa_index": str(gfa_index_path) if gfa_index_path is not None else None,
        },
        "scan_stats": stats.as_summary_dict(),
        "output_dir": str(outdir),
    }
    _write_json(outdir / "run.json", data)
    log.info("Wrote %s/run.json", outdir)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _format_segment_allele_key(contig: str, seg_start: int, seg_name: str) -> str:
    """Return a compact allele key for the hits.tsv ``allele_key`` column.

    Format: ``contig:1based_start:SEG:seg_name``

    The ``SEG:`` prefix distinguishes GFA keys from VCF keys (which use
    ``ref:alt`` notation).  1-based start matches the VCF POS convention so
    coordinates are human-readable.
    """
    return f"{contig}:{seg_start + 1}:SEG:{seg_name}"


def _traversal_to_evidence_class(status: str, cohort_role: str) -> str:
    """Map traversal status and cohort role to an evidence class string."""
    if status == "missing":
        return "uninformative"
    if cohort_role == "target":
        return "support" if status == "traverses" else "absence"
    else:  # off_target
        return "contradiction" if status == "traverses" else "absence"


def _segment_length_class(length: int) -> str:
    """Return a descriptive size class for a coordinate-tagged graph segment."""
    if length <= 0:
        return "unknown"
    if length == 1:
        return "snp_like"
    if length < 50:
        return "small_indel_like"
    if length < 1000:
        return "sv_like"
    return "large_sv_like"


def _graph_segment_interpretation(pattern: AllelePattern) -> str:
    """Explain what the GFA segment call means without VCF-style overclaiming."""
    if pattern.offtarget_missing_n == pattern.offtarget_total_n:
        return (
            "Targets traverse this graph segment; off-target coordinate-overlapping "
            "graph coverage was not observed, so this is private-node evidence "
            "rather than confirmed alternate-path absence."
        )
    if pattern.offtarget_missing_n == 0 and pattern.offtarget_support_n == 0:
        return (
            "Targets traverse this graph segment; off-targets have "
            "coordinate-overlapping graph coverage but do not traverse this same "
            "segment."
        )
    return (
        "Targets traverse this graph segment; off-target same-segment traversal "
        "and coordinate coverage are mixed or incomplete."
    )


def _parse_region(region: str) -> tuple[str, int | None, int | None]:
    """Parse a ``contig:start-end`` region string (1-based samtools format).

    Returns:
        (contig, 0-based start, 0-based exclusive end)
    """
    parts = region.split(":")
    ctg = parts[0]
    if len(parts) < 2:
        return ctg, None, None
    coord_parts = parts[1].split("-")
    if len(coord_parts) != 2:
        return ctg, None, None
    try:
        start_1based = int(coord_parts[0])
        end_1based = int(coord_parts[1])
        return ctg, start_1based - 1, end_1based
    except ValueError:
        return ctg, None, None
