"""VCF-first private-allele discovery backend.

Implements the complete ``privy scan`` workflow for VCF inputs:

    1.  Validate inputs (VCF index, sample overlap with cohort).
    2.  Open VCF and read header metadata.
    3.  Warn about cohort samples absent from the VCF.
    4.  Stream through VCF by contig — one record at a time.
    5.  Apply FILTER/QUAL/multiallelic filters.
    6.  Enumerate alternate alleles; count cohort support via
        :func:`~privy.io.vcf.extract_cohort_counts`.
    7.  Classify each allele with :func:`~privy.core.patterns.classify_strictness`.
    8.  Accumulate passing :class:`HitRecord` objects.
    9.  Score and rank hits.
    10. Merge passing loci into candidate regions.
    11. Write all output files.

Design note — in-memory accumulation:
    Hits are accumulated in a list before scoring and writing.  This allows
    ranking (which requires all final scores to exist) and region merging
    (which requires all passing loci to exist).  The list is much smaller
    than the VCF itself — typically tens of thousands of rows at most — so
    this is acceptable for plant pangenome datasets.

    The VCF scan loop itself is fully streaming: only one
    ``pysam.VariantRecord`` is live in memory at a time.
"""

from __future__ import annotations

import logging
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
from privy.io.jsonio import write_run_json
from privy.io.tsv import (
    EVIDENCE_COLUMNS,
    HITS_COLUMNS,
    QC_COLUMNS,
    REGIONS_COLUMNS,
    SAMPLE_SUPPORT_COLUMNS,
    TsvWriter,
)
from privy.core.evidence import EvidenceRecord
from privy.io.vcf import (
    Genotype,
    classify_variant_type,
    extract_cohort_counts,
    format_allele_key,
    get_vcf_contigs,
    get_vcf_samples,
    stream_vcf_records,
    validate_vcf_index,
)
from privy.utils.metrics import ScanStats
from privy.utils.misc import now_iso

log = logging.getLogger("privy.backends.vcf_scan")


# ---------------------------------------------------------------------------
# Internal data structures
# ---------------------------------------------------------------------------

@dataclass
class HitRecord:
    """Internal record for one passing allele, produced during the scan loop.

    Bundles the locus geometry, allele pattern, and per-sample genotypes
    before scoring.  Not exposed outside this module.
    """

    locus_id: str
    locus: Locus
    allele_key: str
    variant_type: str
    variant_qual: float | None
    pattern: AllelePattern
    # sample_id → GT tuple (pysam allele indices), captured for sample_support.tsv
    sample_genotypes: dict[str, Genotype]
    # Raw alleles from VCF, needed by the BAM support layer for SNP pileup
    ref_allele: str = ""
    alt_allele: str = ""
    # Scores populated after the scan loop
    discovery_score: float = 0.0
    support_score: float = 0.0
    penalty_score: float = 0.0
    final_score: float = 0.0
    rank: int = 0


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def run_vcf_scan(
    vcf: Path | None,
    cohort: CohortDefinition,
    cfg: PrivyConfig,
    outdir: Path,
    mode: str = "private_allele",
    bam: list[Path] | None = None,
    bam_manifest: Path | None = None,
    gfa: Path | None = None,
    xmfa: Path | None = None,
    region: str | None = None,
    contig: str | None = None,
    write_hits: bool = True,
    write_regions: bool = True,
    write_evidence: bool = True,
    write_sample_support: bool = True,
    write_qc: bool = True,
    write_run_json: bool = True,
    threads: int = 1,
) -> None:
    """Run the VCF-first private-allele scan.

    Args:
        vcf: Path to indexed multisample VCF.  Required unless *xmfa* is
             provided as an alternate primary input.
        cohort: Validated :class:`~privy.core.cohort.CohortDefinition`.
        cfg: Resolved :class:`~privy.core.config.PrivyConfig`.
        outdir: Output directory (must already exist).
        mode: Discovery mode.  Only ``"private_allele"`` is implemented
              in Phase 2.
        bam: Optional list of BAM files (ignored in Phase 2).
        bam_manifest: Optional BAM manifest (ignored in Phase 2).
        gfa: Optional GFA file (ignored in Phase 2).
        xmfa: Optional XMFA file (ignored in Phase 2).
        region: Optional region string ``contig:start-end``.
        contig: Optional single-contig restriction.
        write_hits: Write ``hits.tsv``.
        write_regions: Write ``regions.tsv``.
        write_evidence: Write ``evidence.tsv``.
        write_sample_support: Write ``sample_support.tsv``.
        write_qc: Write ``qc.tsv``.
        write_run_json: Write ``run.json``.
        threads: Worker threads (serial only in Phase 2).
    """
    start_time = now_iso()

    if mode != "private_allele":
        raise NotImplementedError(
            f"Discovery mode {mode!r} is not yet implemented.  "
            "Only 'private_allele' is available in Phase 2."
        )

    if vcf is None:
        raise ValueError("vcf path is required for private_allele mode.")

    # ── Step 1: Validate inputs ──────────────────────────────────────────────
    _validate_vcf_inputs(vcf, cohort, cfg)

    # ── Step 2: Read VCF header ──────────────────────────────────────────────
    vcf_samples = get_vcf_samples(vcf)
    vcf_contigs = get_vcf_contigs(vcf)

    # ── Step 3: Validate cohort vs. VCF samples ──────────────────────────────
    vcf_sample_set = set(vcf_samples)
    active_targets = [s for s in cohort.targets if s in vcf_sample_set]
    active_offtargets = [s for s in cohort.off_targets if s in vcf_sample_set]

    missing_targets = [s for s in cohort.targets if s not in vcf_sample_set]
    missing_offtargets = [s for s in cohort.off_targets if s not in vcf_sample_set]

    if missing_targets:
        log.warning(
            "Target samples not found in VCF (will count as missing): %s",
            ", ".join(missing_targets),
        )
    if missing_offtargets:
        log.warning(
            "Off-target samples not found in VCF (will count as missing): %s",
            ", ".join(missing_offtargets),
        )

    if not active_targets:
        raise ValueError(
            "No target samples from the cohort definition were found in the VCF header.  "
            f"Cohort targets: {list(cohort.targets)}.  "
            f"VCF samples (first 10): {vcf_samples[:10]}."
        )
    if not active_offtargets:
        raise ValueError(
            "No off-target samples from the cohort definition were found in the VCF header.  "
            f"Cohort off-targets: {list(cohort.off_targets)}."
        )

    log.info(
        "Cohort validated | targets=%d/%d in VCF | off-targets=%d/%d in VCF",
        len(active_targets), cohort.n_targets,
        len(active_offtargets), cohort.n_off_targets,
    )

    # ── Step 4: Determine contigs to scan ────────────────────────────────────
    if contig is not None:
        contigs_to_scan = [contig]
    elif region is not None:
        contigs_to_scan = [_parse_region_contig(region)]
    else:
        contigs_to_scan = vcf_contigs if vcf_contigs else [None]  # type: ignore[list-item]

    # ── Steps 5–8: Scan VCF and accumulate hits ───────────────────────────────
    stats = ScanStats(
        n_target_samples=len(active_targets),
        n_offtarget_samples=len(active_offtargets),
    )

    # All samples from the cohort (present or absent from VCF) are used in
    # classify_strictness so that cohort totals are stable.
    all_targets = list(cohort.targets)
    all_offtargets = list(cohort.off_targets)

    log.info(
        "Starting scan | contigs=%d | mode=%s | targets=%d | off-targets=%d",
        len(contigs_to_scan), mode, cohort.n_targets, cohort.n_off_targets,
    )

    hit_records: list[HitRecord] = []
    locus_counter = 0

    for ctg in contigs_to_scan:
        ctg_start: int | None = None
        ctg_end: int | None = None

        if region is not None and ctg == _parse_region_contig(region):
            ctg_start, ctg_end = _parse_region_coords(region)

        ctg_hits = _scan_contig(
            vcf_path=vcf,
            contig=ctg,
            start=ctg_start,
            end=ctg_end,
            all_targets=all_targets,
            all_offtargets=all_offtargets,
            cfg=cfg,
            stats=stats,
            locus_counter_start=locus_counter,
        )
        hit_records.extend(ctg_hits)
        locus_counter += len(ctg_hits)
        stats.n_contigs_scanned += 1
        if ctg_hits:
            log.info("  %s: %d hits", ctg, len(ctg_hits))

    stats.loci_emitted = len(hit_records)
    log.info(
        "Scan complete | hits=%d | contradicted=%d | skipped_filter=%d | skipped_qual=%d",
        len(hit_records), stats.alleles_contradicted,
        stats.records_skipped_filter, stats.records_skipped_qual,
    )

    # ── Optional BAM support layer ────────────────────────────────────────────
    bam_result = None
    if bam or bam_manifest:
        from privy.backends.bam_support import (  # noqa: PLC0415
            HitLocusInfo,
            annotate_loci_with_bam,
            resolve_bam_sample_pairs,
        )
        bam_sample_pairs = resolve_bam_sample_pairs(bam, bam_manifest, cohort)
        if bam_sample_pairs and hit_records:
            loci_info = [
                HitLocusInfo(
                    locus_id=hr.locus_id,
                    contig=hr.locus.contig,
                    start=hr.locus.start,
                    end=hr.locus.end,
                    variant_type=hr.variant_type,
                    ref_allele=hr.ref_allele,
                    alt_allele=hr.alt_allele,
                )
                for hr in hit_records
            ]
            bam_result = annotate_loci_with_bam(
                loci_info, bam_sample_pairs, cohort, cfg.bam
            )

    # ── Step 9: Score hits ───────────────────────────────────────────────────
    _score_hit_records(
        hit_records, cfg,
        support_scores=bam_result.support_score_by_locus if bam_result else None,
    )

    # ── Step 10: Merge to regions ─────────────────────────────────────────────
    hit_loci = [hr.locus for hr in hit_records]
    hit_by_id: dict[str, HitRecord] = {hr.locus_id: hr for hr in hit_records}

    region_loci: list[Locus] = []
    if cfg.scan.merge_distance >= 0:
        region_loci = merge_loci_to_regions(
            hit_loci,
            merge_distance=cfg.scan.merge_distance,
            same_variant_class_only=cfg.scan.same_variant_class_only,
            region_id_prefix="REGION",
        )

    stats.regions_emitted = len(region_loci)
    log.info("Region merging complete | regions=%d (merge_distance=%d bp)",
             len(region_loci), cfg.scan.merge_distance)

    # ── Step 11: Write outputs ────────────────────────────────────────────────
    if write_hits:
        _write_hits_tsv(hit_records, outdir)

    if write_regions:
        _write_regions_tsv(region_loci, hit_by_id, outdir)

    if write_evidence:
        _write_evidence_tsv(
            hit_records, outdir,
            bam_records=bam_result.evidence_records if bam_result else None,
        )

    if write_sample_support:
        _write_sample_support_tsv(
            hit_records, cohort, outdir,
            bam_metrics=bam_result.bam_metrics if bam_result else None,
        )

    if write_qc:
        _write_qc_tsv(stats, outdir)

    if write_run_json:
        _write_run_json(
            cfg=cfg,
            cohort=cohort,
            stats=stats,
            vcf_path=vcf,
            outdir=outdir,
            start_time=start_time,
            end_time=now_iso(),
        )

    log.info("Outputs written to %s", outdir)


# ---------------------------------------------------------------------------
# Scan loop
# ---------------------------------------------------------------------------

def _scan_contig(
    vcf_path: Path,
    contig: str | None,
    start: int | None,
    end: int | None,
    all_targets: list[str],
    all_offtargets: list[str],
    cfg: PrivyConfig,
    stats: ScanStats,
    locus_counter_start: int,
) -> list[HitRecord]:
    """Stream through one contig and return passing HitRecord objects."""
    hits: list[HitRecord] = []
    locus_n = locus_counter_start

    for typed_record in stream_vcf_records(vcf_path, contig=contig, start=start, end=end):
        stats.records_evaluated += 1

        # ── FILTER check ─────────────────────────────────────────────────────
        if cfg.scan.pass_only:
            filters = list(typed_record.filter)
            if filters and filters != ["PASS"] and "PASS" not in filters:
                stats.records_skipped_filter += 1
                continue

        # ── QUAL check ───────────────────────────────────────────────────────
        qual = typed_record.qual
        if cfg.scan.min_qual is not None and qual is not None:
            if qual < cfg.scan.min_qual:
                stats.records_skipped_qual += 1
                continue

        # ── Multiallelic check ───────────────────────────────────────────────
        alts = typed_record.alts
        if alts is None:
            continue  # monomorphic reference — skip
        if not cfg.scan.allow_multiallelic and len(alts) > 1:
            stats.records_skipped_multiallelic += 1
            continue

        ref = typed_record.ref
        pos = typed_record.pos  # 1-based VCF POS
        chrom = typed_record.chrom

        # ── Per-alt allele evaluation ─────────────────────────────────────────
        for alt_index, alt in enumerate(alts):
            if alt is None:
                continue  # symbolic padding allele

            stats.alleles_evaluated += 1

            # Count cohort support
            (ts_n, tt_n, os_n, ot_n, tm_n, om_n) = extract_cohort_counts(
                record=typed_record,
                target_samples=all_targets,
                offtarget_samples=all_offtargets,
                alt_index=alt_index,
            )

            # Classify
            pattern = build_allele_pattern(
                allele_key=format_allele_key(chrom, pos, ref, alt),
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

            if not pattern.pattern_pass:
                continue

            stats.alleles_passed += 1

            # 0-based half-open coordinates for Locus
            locus_start = pos - 1
            locus_end = pos - 1 + len(ref)

            locus_n += 1
            locus_id = f"PPX{locus_n:08d}"

            locus = Locus(
                locus_id=locus_id,
                contig=chrom,
                start=locus_start,
                end=locus_end,
                locus_type=LocusType(classify_variant_type(ref, alt)),
                primary_source=PrimarySource.VCF,
                source_ids=[locus_id],  # used for region reconstruction
            )

            # Capture sample genotypes for sample_support.tsv
            sample_genotypes: dict[str, Genotype] = {}
            for s in all_targets + all_offtargets:
                try:
                    sample_genotypes[s] = typed_record.samples[s]["GT"] or (None,)
                except (KeyError, TypeError):
                    sample_genotypes[s] = (None,)

            hits.append(HitRecord(
                locus_id=locus_id,
                locus=locus,
                allele_key=pattern.allele_key,
                variant_type=classify_variant_type(ref, alt),
                variant_qual=qual,
                pattern=pattern,
                sample_genotypes=sample_genotypes,
                ref_allele=ref,
                alt_allele=alt,
            ))

    return hits


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------

def _score_hit_records(
    hit_records: list[HitRecord],
    cfg: PrivyConfig,
    support_scores: dict[str, float] | None = None,
) -> None:
    """Score all hits in-place using the configured scoring weights."""
    scored: list[ScoredHit] = []

    for hr in hit_records:
        ds = compute_discovery_score(
            hr.pattern,
            variant_qual=hr.variant_qual,
            discovery_weight=cfg.scoring.discovery_weight,
        )
        ps = compute_penalty_score(hr.pattern, penalty_weight=cfg.scoring.penalty_weight)
        raw_ss = support_scores.get(hr.locus_id, 0.0) if support_scores else 0.0
        ss = round(raw_ss * cfg.scoring.support_weight, 6)
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

    # Rank in-place
    ranked = rank_scored_hits(scored)
    rank_map = {s.locus_id: s.rank for s in ranked}
    for hr in hit_records:
        hr.rank = rank_map.get(hr.locus_id, 0)


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def _write_hits_tsv(hit_records: list[HitRecord], outdir: Path) -> None:
    """Write ``hits.tsv`` — one row per passing locus, ranked by final_score."""
    path = outdir / "hits.tsv"
    sorted_hits = sorted(hit_records, key=lambda h: h.rank)
    with TsvWriter(path, HITS_COLUMNS) as w:
        for hr in sorted_hits:
            w.write_row({
                "locus_id":          hr.locus_id,
                "contig":            hr.locus.contig,
                "start":             hr.locus.start,
                "end":               hr.locus.end,
                "variant_type":      hr.variant_type,
                "allele_key":        hr.allele_key,
                "target_support_n":  hr.pattern.target_support_n,
                "target_total_n":    hr.pattern.target_total_n,
                "offtarget_support_n": hr.pattern.offtarget_support_n,
                "offtarget_total_n": hr.pattern.offtarget_total_n,
                "target_missing_n":  hr.pattern.target_missing_n,
                "offtarget_missing_n": hr.pattern.offtarget_missing_n,
                "strictness_class":  hr.pattern.strictness_class.value,
                "discovery_score":   hr.discovery_score,
                "support_score":     hr.support_score,
                "penalty_score":     hr.penalty_score,
                "final_score":       hr.final_score,
            })
    log.info("Wrote %s (%d rows)", path, len(hit_records))


def _write_regions_tsv(
    region_loci: list[Locus],
    hit_by_id: dict[str, HitRecord],
    outdir: Path,
) -> None:
    """Write ``regions.tsv`` — one row per merged candidate region."""
    path = outdir / "regions.tsv"
    with TsvWriter(path, REGIONS_COLUMNS) as w:
        for region in region_loci:
            constituent = [hit_by_id[sid] for sid in region.source_ids if sid in hit_by_id]
            if not constituent:
                continue

            variant_types = sorted(set(c.variant_type for c in constituent))
            strictness_counter: Counter[str] = Counter(
                c.pattern.strictness_class.value for c in constituent
            )
            dominant_strictness = strictness_counter.most_common(1)[0][0]

            n_target_complete = sum(
                1 for c in constituent if c.pattern.target_missing_n == 0
            )
            target_consistency = n_target_complete / len(constituent)

            n_offtarget_excluded = sum(
                1 for c in constituent if c.pattern.offtarget_support_n == 0
            )
            offtarget_exclusion = n_offtarget_excluded / len(constituent)

            mean_score = sum(c.final_score for c in constituent) / len(constituent)

            w.write_row({
                "region_id":               region.locus_id,
                "contig":                  region.contig,
                "start":                   region.start,
                "end":                     region.end,
                "n_loci":                  len(constituent),
                "variant_types":           ",".join(variant_types),
                "dominant_strictness_class": dominant_strictness,
                "target_consistency":      round(target_consistency, 4),
                "offtarget_exclusion":     round(offtarget_exclusion, 4),
                "final_score":             round(mean_score, 6),
            })
    log.info("Wrote %s (%d rows)", path, len(region_loci))


def _write_evidence_tsv(
    hit_records: list[HitRecord],
    outdir: Path,
    bam_records: list[EvidenceRecord] | None = None,
) -> None:
    """Write ``evidence.tsv`` — VCF evidence plus any BAM evidence records."""
    path = outdir / "evidence.tsv"
    vcf_count = len(hit_records)
    bam_count = len(bam_records) if bam_records else 0
    with TsvWriter(path, EVIDENCE_COLUMNS) as w:
        for hr in hit_records:
            w.write_row({
                "locus_id":       hr.locus_id,
                "source_type":    "vcf",
                "sample_id":      "",
                "evidence_class": "support" if hr.pattern.pattern_pass else "contradiction",
                "metric_name":    "allele_pattern",
                "metric_value":   hr.final_score,
                "details":        hr.pattern.pattern_reason,
            })
        if bam_records:
            for er in bam_records:
                w.write_row({
                    "locus_id":       er.locus_id,
                    "source_type":    er.source_type.value,
                    "sample_id":      er.sample_id or "",
                    "evidence_class": er.evidence_class.value,
                    "metric_name":    er.metric_name,
                    "metric_value":   er.metric_value,
                    "details":        er.provenance,
                })
    log.info("Wrote %s (%d VCF + %d BAM rows)", path, vcf_count, bam_count)


def _write_sample_support_tsv(
    hit_records: list[HitRecord],
    cohort: CohortDefinition,
    outdir: Path,
    bam_metrics: dict[tuple[str, str], dict[str, str]] | None = None,
) -> None:
    """Write ``sample_support.tsv`` — one row per sample per passing locus."""
    path = outdir / "sample_support.tsv"
    with TsvWriter(path, SAMPLE_SUPPORT_COLUMNS) as w:
        for hr in sorted(hit_records, key=lambda h: h.rank):
            for sample, gt in hr.sample_genotypes.items():
                if cohort.is_ignored(sample):
                    continue
                cohort_role = "target" if cohort.is_target(sample) else "off_target"
                gt_str = _format_gt(gt)
                alt_supported = _gt_supports_alt(gt)
                evidence_class = _gt_to_evidence_class(gt, alt_supported, cohort_role)
                metrics = bam_metrics.get((hr.locus_id, sample), {}) if bam_metrics else {}
                w.write_row({
                    "locus_id":         hr.locus_id,
                    "sample_id":        sample,
                    "cohort_role":      cohort_role,
                    "genotype":         gt_str,
                    "allele_supported": str(alt_supported).lower(),
                    "depth":            metrics.get("depth", "NA"),
                    "allele_fraction":  metrics.get("allele_fraction", "NA"),
                    "evidence_class":   evidence_class,
                })
    n_rows = len(hit_records) * len(hit_records[0].sample_genotypes) if hit_records else 0
    log.info("Wrote %s (~%d rows)", path, n_rows)


def _write_qc_tsv(stats: ScanStats, outdir: Path) -> None:
    """Write ``qc.tsv``."""
    path = outdir / "qc.tsv"
    with TsvWriter(path, QC_COLUMNS) as w:
        w.write_rows(stats.as_qc_rows())
    log.info("Wrote %s", path)


def _write_run_json(
    cfg: PrivyConfig,
    cohort: CohortDefinition,
    stats: ScanStats,
    vcf_path: Path | None,
    outdir: Path,
    start_time: str,
    end_time: str,
) -> None:
    """Write ``run.json`` — full resolved configuration and run metadata."""
    from privy import __version__  # noqa: PLC0415

    data = {
        "privy_version": __version__,
        "start_time": start_time,
        "end_time": end_time,
        "project_name": cfg.project_name,
        "config": cfg.as_run_dict(),
        "cohort": {
            "targets": list(cohort.targets),
            "off_targets": list(cohort.off_targets),
            "ignored_samples": list(cohort.ignored_samples),
            "n_targets": cohort.n_targets,
            "n_off_targets": cohort.n_off_targets,
        },
        "inputs": {
            "vcf": str(vcf_path) if vcf_path else None,
        },
        "scan_stats": stats.as_summary_dict(),
        "output_dir": str(outdir),
    }
    write_run_json(outdir / "run.json", data)
    log.info("Wrote %s/run.json", outdir)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _validate_vcf_inputs(vcf: Path, cohort: CohortDefinition, cfg: PrivyConfig) -> None:
    """Validate that the VCF file exists and is indexed."""
    if not vcf.exists():
        raise FileNotFoundError(f"VCF file not found: {vcf}")
    validate_vcf_index(vcf)


def _parse_region_contig(region: str) -> str:
    """Extract contig name from a ``contig:start-end`` string."""
    return region.split(":")[0]


def _parse_region_coords(region: str) -> tuple[int | None, int | None]:
    """Extract 0-based start/end from a ``contig:start-end`` string.

    The input is assumed to be 1-based coordinates (samtools convention).
    Returns 0-based half-open coordinates for pysam.
    """
    parts = region.split(":")
    if len(parts) < 2:
        return None, None
    coord_str = parts[1]
    coord_parts = coord_str.split("-")
    if len(coord_parts) != 2:
        return None, None
    try:
        start_1based = int(coord_parts[0])
        end_1based = int(coord_parts[1])
        return start_1based - 1, end_1based  # convert to 0-based half-open
    except ValueError:
        return None, None


def _format_gt(gt: Genotype | None) -> str:
    """Format a pysam GT tuple as a human-readable string, e.g. ``0/1``."""
    if gt is None:
        return "./."
    return "/".join("." if a is None else str(a) for a in gt)


def _gt_supports_alt(gt: Genotype | None) -> bool:
    """Return True if the GT contains any non-zero, non-None allele."""
    if gt is None:
        return False
    return any(a is not None and a > 0 for a in gt)


def _gt_to_evidence_class(
    gt: Genotype | None,
    alt_supported: bool,
    cohort_role: str,
) -> str:
    """Map GT and cohort role to an evidence class string for sample_support.tsv."""
    if gt is None or any(a is None for a in gt) if gt else True:
        return "uninformative"
    if cohort_role == "target":
        return "support" if alt_supported else "absence"
    else:  # off_target
        return "contradiction" if alt_supported else "absence"
