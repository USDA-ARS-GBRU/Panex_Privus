#!/usr/bin/env python3
"""Find group-private variant positions from a multi-sample VCF.

Designed for whole-genome assembly alignments (e.g. minimap2 asm-to-ref BAMs
converted to VCF with bcftools mpileup + call). Coordinates are taken directly
from the VCF POS column and match what IGV displays.

A position qualifies when:
  1. Every target sample with a called, homozygous genotype shares the same allele
  2. At least one non-target sample has a called genotype
  3. Not every non-target carries the target allele (site must be informative)

Heterozygous target genotypes are treated as "partial" (reduces strictness).
Missing genotypes (./.) are skipped for both targets and non-targets.

Strictness levels (per interval: worst site in the interval determines the level):
  1  strict                          all targets hom + agree; no non-target carries it
  2  nontarget_shared                all targets hom + agree; ≥1 non-target also carries it
  3  target_partial                  ≥1 target het/missing; hom targets agree; no non-target carries it
  4  target_partial_nontarget_shared ≥1 target het/missing AND ≥1 non-target carries the allele

Typical workflow:
  # 1. Call variants per genome (disable read-level quality filters for assemblies)
  for bam in *.bam; do
      sample=$(basename "$bam" .bam)
      bcftools mpileup -Ou -f ref.fasta --no-BAQ -Q 0 -q 0 "$bam" \\
          | bcftools call -mv -Oz -o "${sample}.vcf.gz"
      bcftools index "${sample}.vcf.gz"
  done

  # 2. Merge into a single multi-sample VCF
  bcftools merge -Oz -o all_genomes.vcf.gz *.vcf.gz
  bcftools index all_genomes.vcf.gz

  # 3. Run this script
  python find_group_private_vcf_regions.py \\
      --vcf all_genomes.vcf.gz \\
      --targets Minsoy MiyakoWhite Kingawa PI407303 \\
      --region Gm15 \\
      --out-prefix results/group01_private_positions
"""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

import pysam

STRICTNESS_LEVELS: Dict[int, str] = {
    1: "strict",
    2: "nontarget_shared",
    3: "target_partial",
    4: "target_partial_nontarget_shared",
}


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class SiteHit:
    chrom: str
    ref_pos: int                    # 1-based, matches VCF POS and IGV display
    ref_allele: str
    target_allele: str
    classification: str             # snp or indel
    n_targets_observed: int         # hom targets agreeing at this site
    n_nontarget_observed: int       # non-targets with any called genotype
    n_nontarget_carries: int        # non-targets that carry the target allele
    n_nontarget_distinct: int       # distinct non-target alleles (excl. target allele)
    nontarget_alleles: Tuple[str, ...]
    strictness: str
    strictness_level: int


@dataclass
class IntervalHit:
    interval_id: str
    chrom: str
    strictness: str
    strictness_level: int
    start: int
    end: int
    span_bp: int
    n_sites: int
    classification: str             # snp, indel, or mixed
    targets: str
    target_alleles: str             # unique target alleles observed across sites
    ref_alleles: str                # unique ref alleles across sites
    nontarget_alleles: str          # unique non-target alleles (excl. target allele)
    n_targets_expected: int
    n_targets_observed: int         # minimum hom-target count across sites in interval
    min_nontarget_observed: int
    max_nontarget_observed: int
    max_nontarget_carries: int      # worst-case non-target carry count (0 for strict)


@dataclass
class ScanStats:
    total_sites: int = 0
    sites_in_window: int = 0
    sites_skipped_no_hom_target: int = 0
    sites_skipped_targets_disagree: int = 0
    sites_skipped_no_nontarget: int = 0
    sites_skipped_all_nontarget_carry: int = 0
    qualifying_sites: int = 0
    qualifying_intervals: int = 0


# ---------------------------------------------------------------------------
# Core helpers
# ---------------------------------------------------------------------------

def _resolve_allele(gt_index: int, ref: str, alts: Tuple[str, ...]) -> str:
    if gt_index == 0:
        return ref
    return alts[gt_index - 1]


def _get_sample_call(sample_rec, ref: str, alts: Tuple[str, ...]):
    """Parse a VCF sample genotype field.

    Returns:
        (allele_str, True)       homozygous call
        (frozenset[str], False)  heterozygous call (set of allele strings)
        None                     missing or uncallable
    """
    gt = sample_rec["GT"]
    if gt is None:
        return None
    indices = [i for i in gt if i is not None]
    if not indices:
        return None
    alleles = [_resolve_allele(i, ref, alts) for i in indices]
    unique = set(alleles)
    if len(unique) == 1:
        return (alleles[0], True)
    return (frozenset(unique), False)


def _classify(ref_allele: str, target_allele: str) -> str:
    if len(ref_allele) == 1 and len(target_allele) == 1:
        return "snp"
    return "indel"


# ---------------------------------------------------------------------------
# VCF scanning
# ---------------------------------------------------------------------------

def scan_vcf(
    vcf_path: str,
    targets: Sequence[str],
    region: Optional[str],
    start: Optional[int],
    end: Optional[int],
) -> Tuple[List[SiteHit], ScanStats, List[str]]:

    stats = ScanStats()
    hits: List[SiteHit] = []

    vcf = pysam.VariantFile(vcf_path)
    all_samples = list(vcf.header.samples)

    missing = [s for s in targets if s not in all_samples]
    if missing:
        raise SystemExit(
            f"ERROR: target samples not found in VCF: {', '.join(missing)}\n"
            f"Samples in VCF: {', '.join(all_samples)}"
        )

    other_samples = [s for s in all_samples if s not in targets]

    records = vcf.fetch(region=region) if region else vcf.fetch()

    for rec in records:
        pos = rec.pos          # 1-based, matches IGV display
        chrom = rec.chrom
        ref = rec.ref
        alts = rec.alts or ()

        stats.total_sites += 1

        if start is not None and pos < start:
            continue
        if end is not None and pos > end:
            continue
        stats.sites_in_window += 1

        # ---- target genotypes ----
        hom_alleles: List[str] = []
        for sid in targets:
            call = _get_sample_call(rec.samples[sid], ref, alts)
            if call is None:
                continue
            if call[1]:                     # homozygous
                hom_alleles.append(call[0])
            # het targets: not counted as hom, will lower n_targets_observed

        if not hom_alleles:
            stats.sites_skipped_no_hom_target += 1
            continue

        if len(set(hom_alleles)) != 1:
            stats.sites_skipped_targets_disagree += 1
            continue

        target_allele = hom_alleles[0]
        n_targets_observed = len(hom_alleles)
        all_targets_observed = (n_targets_observed == len(targets))

        # ---- non-target genotypes ----
        nontarget_observed: List[str] = []
        nontarget_carries: List[str] = []
        nontarget_allele_set: Set[str] = set()

        for sid in other_samples:
            call = _get_sample_call(rec.samples[sid], ref, alts)
            if call is None:
                continue
            nontarget_observed.append(sid)
            if call[1]:                     # hom
                nontarget_allele_set.add(call[0])
                if call[0] == target_allele:
                    nontarget_carries.append(sid)
            else:                           # het — carries if target allele is one of the two
                nontarget_allele_set.update(call[0])
                if target_allele in call[0]:
                    nontarget_carries.append(sid)

        if not nontarget_observed:
            stats.sites_skipped_no_nontarget += 1
            continue

        # skip positions where every called non-target carries the target allele
        if len(nontarget_carries) == len(nontarget_observed):
            stats.sites_skipped_all_nontarget_carry += 1
            continue

        # ---- strictness ----
        strict_nontarget = len(nontarget_carries) == 0
        if all_targets_observed and strict_nontarget:
            level = 1
        elif all_targets_observed and not strict_nontarget:
            level = 2
        elif not all_targets_observed and strict_nontarget:
            level = 3
        else:
            level = 4

        nontarget_other = nontarget_allele_set - {target_allele}

        hits.append(SiteHit(
            chrom=chrom,
            ref_pos=pos,
            ref_allele=ref,
            target_allele=target_allele,
            classification=_classify(ref, target_allele),
            n_targets_observed=n_targets_observed,
            n_nontarget_observed=len(nontarget_observed),
            n_nontarget_carries=len(nontarget_carries),
            n_nontarget_distinct=len(nontarget_other),
            nontarget_alleles=tuple(sorted(nontarget_other)),
            strictness=STRICTNESS_LEVELS[level],
            strictness_level=level,
        ))
        stats.qualifying_sites += 1

    return hits, stats, all_samples


# ---------------------------------------------------------------------------
# Interval merging
# ---------------------------------------------------------------------------

def merge_sites_to_intervals(
    sites: Sequence[SiteHit],
    n_targets_expected: int,
    target_str: str,
) -> List[IntervalHit]:
    if not sites:
        return []

    sorted_sites = sorted(sites, key=lambda s: (s.chrom, s.ref_pos))
    intervals: List[IntervalHit] = []
    current: List[SiteHit] = [sorted_sites[0]]

    def finalize(chunk: Sequence[SiteHit], num: int) -> IntervalHit:
        level = max(s.strictness_level for s in chunk)
        site_classes = sorted({s.classification for s in chunk})
        cls = site_classes[0] if len(site_classes) == 1 else "mixed"
        return IntervalHit(
            interval_id=f"region_{num:06d}",
            chrom=chunk[0].chrom,
            strictness=STRICTNESS_LEVELS[level],
            strictness_level=level,
            start=chunk[0].ref_pos,
            end=chunk[-1].ref_pos,
            span_bp=chunk[-1].ref_pos - chunk[0].ref_pos + 1,
            n_sites=len(chunk),
            classification=cls,
            targets=target_str,
            target_alleles=",".join(sorted({s.target_allele for s in chunk})),
            ref_alleles=",".join(sorted({s.ref_allele for s in chunk})),
            nontarget_alleles=",".join(sorted({a for s in chunk for a in s.nontarget_alleles})),
            n_targets_expected=n_targets_expected,
            n_targets_observed=min(s.n_targets_observed for s in chunk),
            min_nontarget_observed=min(s.n_nontarget_observed for s in chunk),
            max_nontarget_observed=max(s.n_nontarget_observed for s in chunk),
            max_nontarget_carries=max(s.n_nontarget_carries for s in chunk),
        )

    num = 1
    for site in sorted_sites[1:]:
        prev = current[-1]
        if site.chrom == prev.chrom and site.ref_pos == prev.ref_pos + 1:
            current.append(site)
        else:
            intervals.append(finalize(current, num))
            num += 1
            current = [site]
    intervals.append(finalize(current, num))
    return intervals


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def write_intervals_tsv(
    out_path: Path,
    intervals: Sequence[IntervalHit],
) -> None:
    header = [
        "interval_id", "chrom", "strictness", "strictness_level",
        "start", "end", "span_bp", "n_sites", "classification",
        "targets", "target_alleles", "ref_alleles", "nontarget_alleles",
        "n_targets_expected", "n_targets_observed",
        "min_nontarget_observed", "max_nontarget_observed", "max_nontarget_carries",
    ]
    with out_path.open("w") as fh:
        fh.write("\t".join(header) + "\n")
        for row in intervals:
            fh.write("\t".join([
                row.interval_id, row.chrom, row.strictness, str(row.strictness_level),
                str(row.start), str(row.end), str(row.span_bp), str(row.n_sites),
                row.classification, row.targets, row.target_alleles, row.ref_alleles,
                row.nontarget_alleles, str(row.n_targets_expected),
                str(row.n_targets_observed), str(row.min_nontarget_observed),
                str(row.max_nontarget_observed), str(row.max_nontarget_carries),
            ]) + "\n")


def write_summary_tsv(
    out_path: Path,
    stats: ScanStats,
    vcf_path: str,
    targets: Sequence[str],
    all_samples: Sequence[str],
    intervals: Sequence[IntervalHit],
    region: Optional[str],
    start: Optional[int],
    end: Optional[int],
) -> None:
    level_counts = Counter(i.strictness_level for i in intervals)
    class_counts = Counter(i.classification for i in intervals)
    rows = [
        ("vcf",                                 vcf_path),
        ("targets",                             ",".join(targets)),
        ("all_samples",                         ",".join(all_samples)),
        ("n_targets",                           str(len(targets))),
        ("n_nontargets",                        str(len(all_samples) - len(targets))),
        ("region",                              region or ""),
        ("window_start",                        "" if start is None else str(start)),
        ("window_end",                          "" if end is None else str(end)),
        ("total_sites",                         str(stats.total_sites)),
        ("sites_in_window",                     str(stats.sites_in_window)),
        ("sites_skipped_no_hom_target",         str(stats.sites_skipped_no_hom_target)),
        ("sites_skipped_targets_disagree",      str(stats.sites_skipped_targets_disagree)),
        ("sites_skipped_no_nontarget",          str(stats.sites_skipped_no_nontarget)),
        ("sites_skipped_all_nontarget_carry",   str(stats.sites_skipped_all_nontarget_carry)),
        ("qualifying_sites",                    str(stats.qualifying_sites)),
        ("qualifying_intervals",                str(len(intervals))),
        ("qualifying_interval_bp_total",        str(sum(i.span_bp for i in intervals))),
        ("level_1_strict",                      str(level_counts.get(1, 0))),
        ("level_2_nontarget_shared",            str(level_counts.get(2, 0))),
        ("level_3_target_partial",              str(level_counts.get(3, 0))),
        ("level_4_target_partial_nontarget_shared", str(level_counts.get(4, 0))),
        ("class_snp",                           str(class_counts.get("snp", 0))),
        ("class_indel",                         str(class_counts.get("indel", 0))),
        ("class_mixed",                         str(class_counts.get("mixed", 0))),
    ]
    with out_path.open("w") as fh:
        fh.write("metric\tvalue\n")
        for k, v in rows:
            fh.write(f"{k}\t{v}\n")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Scan a multi-sample VCF for group-private variant positions. "
            "Coordinates match IGV when VCF is built from the same reference."
        )
    )
    p.add_argument(
        "--vcf", required=True,
        help="Multi-sample VCF or BCF (plain, gzipped, or bgzipped+indexed)",
    )
    p.add_argument(
        "--targets", required=True, nargs="+",
        help="Sample IDs of the target group (must match VCF sample names exactly)",
    )
    p.add_argument(
        "--out-prefix", required=True,
        help="Output path prefix; writes <prefix>.intervals.tsv and <prefix>.summary.tsv",
    )
    p.add_argument(
        "--region", default=None,
        help=(
            "Genomic region to scan in samtools format, e.g. 'Gm15' or "
            "'Gm15:1000000-2000000'. Requires the VCF to be bgzipped and indexed."
        ),
    )
    p.add_argument(
        "--start", type=int, default=None,
        help="Additional 1-based start position filter (applied after --region)",
    )
    p.add_argument(
        "--end", type=int, default=None,
        help="Additional 1-based end position filter (applied after --region)",
    )
    return p


def main() -> None:
    args = build_argparser().parse_args()

    if args.start is not None and args.end is not None and args.start > args.end:
        raise SystemExit("ERROR: --start must be <= --end")
    if len(set(args.targets)) != len(args.targets):
        raise SystemExit("ERROR: --targets contains duplicate sample IDs")

    hits, stats, all_samples = scan_vcf(
        vcf_path=args.vcf,
        targets=args.targets,
        region=args.region,
        start=args.start,
        end=args.end,
    )

    target_str = ",".join(args.targets)
    intervals = merge_sites_to_intervals(hits, len(args.targets), target_str)
    stats.qualifying_intervals = len(intervals)

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    intervals_path = Path(f"{args.out_prefix}.intervals.tsv")
    summary_path = Path(f"{args.out_prefix}.summary.tsv")

    write_intervals_tsv(intervals_path, intervals)
    write_summary_tsv(
        summary_path, stats, args.vcf, args.targets, all_samples,
        intervals, args.region, args.start, args.end,
    )

    print(f"Samples in VCF : {', '.join(all_samples)}")
    print(f"Targets        : {', '.join(args.targets)}")
    print(f"Qualifying sites    : {stats.qualifying_sites}")
    print(f"Qualifying intervals: {len(intervals)}")
    print(f"Wrote: {intervals_path}")
    print(f"Wrote: {summary_path}")


if __name__ == "__main__":
    main()
