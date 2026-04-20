"""BAM support layer — read-level evidence at candidate loci.

BAM is not a discovery caller.  It queries depth and allele counts at loci
already identified by the VCF or GFA backend and classifies the observations
as SUPPORT, ABSENCE, CONTRADICTION, AMBIGUOUS, or UNINFORMATIVE evidence.

Design notes:
  - :class:`HitLocusInfo` is a lightweight dataclass that decouples this
    module from ``vcf_scan``'s internal ``HitRecord`` type.
  - :class:`BamAnnotationResult` bundles all outputs so the VCF backend
    can integrate them in a single step.
  - UNINFORMATIVE values (low depth or non-SNP loci) are excluded from
    the per-locus support-score mean to avoid penalising loci for
    low-coverage samples.
  - Negative values (CONTRADICTION → -1.0) are passed through to
    :func:`~privy.core.scoring.compute_support_score` unchanged, so
    contradicted loci receive a penalty in the final score.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path

from privy.core.cohort import CohortDefinition
from privy.core.config import BamConfig
from privy.core.evidence import EvidenceClass, EvidenceRecord, SourceType
from privy.core.scoring import compute_support_score
from privy.io.bam import (
    get_bam_sample_name,
    load_bam_manifest,
    query_allele_counts_at_locus,
    query_position_depth,
    validate_bam_index,
)

log = logging.getLogger("privy.backends.bam_support")


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class HitLocusInfo:
    """Lightweight descriptor for a candidate locus passed to the BAM layer."""

    locus_id: str
    contig: str
    start: int       # 0-based
    end: int         # 0-based half-open
    variant_type: str
    ref_allele: str
    alt_allele: str


@dataclass
class BamAnnotationResult:
    """Aggregated outputs from BAM annotation of candidate loci."""

    evidence_records: list[EvidenceRecord]
    support_score_by_locus: dict[str, float]
    # (locus_id, sample_id) → {"depth": str, "allele_fraction": str}
    bam_metrics: dict[tuple[str, str], dict[str, str]] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def resolve_bam_sample_pairs(
    bam_paths: list[Path] | None,
    bam_manifest: Path | None,
    cohort: CohortDefinition,
) -> list[tuple[Path, str]]:
    """Return ``(bam_path, sample_id)`` pairs for all provided BAM inputs.

    Priority: manifest > explicit paths.  For explicit paths the sample name
    is read from the @RG SM header tag; if absent, the filename stem is used
    and a warning is logged.

    Args:
        bam_paths: Optional list of individual BAM paths.
        bam_manifest: Optional TSV manifest (takes precedence over paths).
        cohort: Used only for logging — unknown samples are included but noted.

    Returns:
        List of ``(Path, str)`` pairs, one per BAM/sample entry.
    """
    pairs: list[tuple[Path, str]] = []

    if bam_manifest is not None:
        for row in load_bam_manifest(bam_manifest):
            pairs.append((Path(row["bam_path"]), row["sample_id"]))
        log.info("Loaded %d BAM(s) from manifest %s", len(pairs), bam_manifest)
        return pairs

    if bam_paths:
        for bam_p in bam_paths:
            sample = get_bam_sample_name(bam_p)
            if sample is None:
                sample = bam_p.stem
                log.warning(
                    "BAM %s has no @RG SM tag; using filename stem as sample: %s",
                    bam_p, sample,
                )
            pairs.append((bam_p, sample))
        log.info("Using %d explicitly provided BAM(s)", len(pairs))

    return pairs


def annotate_loci_with_bam(
    loci_info: list[HitLocusInfo],
    bam_sample_pairs: list[tuple[Path, str]],
    cohort: CohortDefinition,
    cfg: BamConfig,
) -> BamAnnotationResult:
    """Query BAM files at candidate loci and classify the evidence.

    For each (locus, BAM) pair:
      - Queries allele counts (SNPs) or depth (indels/SVs) via pysam.
      - Classifies the observation against the sample's cohort role.
      - Excludes UNINFORMATIVE observations from the per-locus score mean.

    Args:
        loci_info: Candidate loci from the prior scan step.
        bam_sample_pairs: ``(bam_path, sample_id)`` pairs to query.
        cohort: Cohort definition (target / off-target roles).
        cfg: BAM configuration (thresholds and quality filters).

    Returns:
        :class:`BamAnnotationResult` containing all evidence records,
        per-locus support scores, and per-(locus, sample) depth/AF metrics.
    """
    evidence_records: list[EvidenceRecord] = []
    actionable_by_locus: dict[str, list[float]] = {li.locus_id: [] for li in loci_info}
    bam_metrics: dict[tuple[str, str], dict[str, str]] = {}

    for bam_path, sample_id in bam_sample_pairs:
        try:
            validate_bam_index(bam_path)
        except FileNotFoundError as exc:
            log.warning("Skipping BAM %s: %s", bam_path, exc)
            continue

        if cohort.is_ignored(sample_id):
            continue

        cohort_role = "target" if cohort.is_target(sample_id) else "off_target"

        for locus in loci_info:
            is_snp = len(locus.ref_allele) == 1 and len(locus.alt_allele) == 1

            if is_snp:
                ref_count, alt_count, other_count = query_allele_counts_at_locus(
                    bam_path,
                    locus.contig,
                    locus.start,
                    locus.ref_allele,
                    locus.alt_allele,
                    min_mapq=cfg.min_mapq,
                    min_baseq=cfg.min_baseq,
                )
                depth = ref_count + alt_count + other_count
                allele_fraction = round(alt_count / depth, 4) if depth > 0 else 0.0
                depth_str = str(depth)
                af_str = str(allele_fraction)
                metric_name = "allele_fraction"
                metric_value = allele_fraction
            else:
                depths = query_position_depth(
                    bam_path, locus.contig, locus.start, locus.end,
                    min_mapq=cfg.min_mapq,
                )
                depth = int(sum(depths) / len(depths)) if depths else 0
                ref_count, alt_count, other_count = depth, 0, 0
                allele_fraction = 0.0
                depth_str = str(depth)
                af_str = "NA"
                metric_name = "depth"
                metric_value = float(depth)

            bam_metrics[(locus.locus_id, sample_id)] = {
                "depth": depth_str,
                "allele_fraction": af_str,
            }

            evidence_class, ev_value = _classify_bam_evidence(
                ref_count=ref_count,
                alt_count=alt_count,
                depth=depth,
                allele_fraction=allele_fraction if is_snp else None,
                cohort_role=cohort_role,
                cfg=cfg,
                is_snp=is_snp,
            )

            qualifiers: dict[str, object] = {
                "depth": depth,
                "ref_count": ref_count,
                "alt_count": alt_count,
                "other_count": other_count,
                "min_mapq": cfg.min_mapq,
            }
            if is_snp:
                qualifiers["allele_fraction"] = allele_fraction
                qualifiers["min_baseq"] = cfg.min_baseq

            evidence_records.append(EvidenceRecord(
                locus_id=locus.locus_id,
                source_type=SourceType.BAM,
                evidence_class=evidence_class,
                metric_name=metric_name,
                metric_value=metric_value,
                sample_id=sample_id,
                qualifiers=qualifiers,
                provenance=f"BAM pileup at {locus.contig}:{locus.start}",
            ))

            if evidence_class != EvidenceClass.UNINFORMATIVE:
                actionable_by_locus[locus.locus_id].append(ev_value)

    support_score_by_locus: dict[str, float] = {
        locus.locus_id: compute_support_score(
            actionable_by_locus[locus.locus_id],
            support_weight=1.0,
        )
        for locus in loci_info
    }

    log.info(
        "BAM annotation complete | loci=%d | bam_samples=%d | evidence_records=%d",
        len(loci_info), len(bam_sample_pairs), len(evidence_records),
    )

    return BamAnnotationResult(
        evidence_records=evidence_records,
        support_score_by_locus=support_score_by_locus,
        bam_metrics=bam_metrics,
    )


# ---------------------------------------------------------------------------
# Evidence classification
# ---------------------------------------------------------------------------

def _classify_bam_evidence(
    ref_count: int,
    alt_count: int,
    depth: int,
    allele_fraction: float | None,
    cohort_role: str,
    cfg: BamConfig,
    is_snp: bool,
) -> tuple[EvidenceClass, float]:
    """Classify a BAM observation into an evidence class and normalised value.

    Classification logic:

    - Depth below ``cfg.min_depth`` → UNINFORMATIVE (0.0), regardless of role.
    - Non-SNP locus → UNINFORMATIVE (0.0); depth alone cannot distinguish alleles.
    - Target sample:
        - alt present (count ≥ min_alt_count **and** AF ≥ allele_fraction_min)
          → SUPPORT (1.0)
        - Otherwise (VCF says present, BAM cannot confirm) → AMBIGUOUS (0.3)
    - Off-target sample:
        - alt present → CONTRADICTION (−1.0)
        - alt absent → ABSENCE (1.0)

    Returns:
        ``(EvidenceClass, normalised_float)`` tuple.
    """
    if depth < cfg.min_depth:
        return EvidenceClass.UNINFORMATIVE, 0.0

    if not is_snp:
        return EvidenceClass.UNINFORMATIVE, 0.0

    af = allele_fraction if allele_fraction is not None else 0.0

    if cohort_role == "target":
        if alt_count >= cfg.min_alt_count and af >= cfg.allele_fraction_min:
            return EvidenceClass.SUPPORT, 1.0
        return EvidenceClass.AMBIGUOUS, 0.3

    # off_target
    if alt_count >= cfg.min_alt_count and af >= cfg.allele_fraction_min:
        return EvidenceClass.CONTRADICTION, -1.0
    return EvidenceClass.ABSENCE, 1.0
