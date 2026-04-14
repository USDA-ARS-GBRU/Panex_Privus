#!/usr/bin/env python3
"""Find group-private variant positions by scanning assembly BAM files directly.

Designed for whole-genome assembly alignments (minimap2 asm-to-ref) where
contigs can be tens of megabases long. Unlike bcftools mpileup, this script
processes one alignment record at a time so memory usage scales with the number
of variant sites found, not the length of the contigs.

Memory profile (approximate):
  - Per genome: ~50 bytes × number of SNP sites (typically < 500 MB for a full genome)
  - Coverage intervals: negligible (one entry per contig alignment)
  Total for 10 assemblies across a full ~1 Gbp genome: typically 2–5 GB

Requires:
  - pysam
  - BAM files sorted and indexed (samtools sort / samtools index)
  - Reference FASTA indexed with samtools faidx

NOTE: Insertions relative to the reference are not reported. SNPs and
deletions are detected. Insertions do not have a unique reference coordinate
and are rarely used as group-private markers.

Usage:
  python find_group_private_bam_regions.py \\
      --samples samples.tsv \\
      --ref ref/Gmax_880_v6.0.fa \\
      --targets Minsoy MiyakoWhite Kingawa PI407303 \\
      --region Gm15 \\
      --out-prefix results/group01_private_positions

  samples.tsv is a two-column tab-separated file (no header):
      Benning    01_bams/Benning_reheader.bam
      Clark      01_bams/Clark_reheader.bam
      ...
"""

from __future__ import annotations

import argparse
import bisect
from collections import defaultdict
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

# Sentinel for a deletion allele (query has no base at this ref position)
DELETION = "*"


# ---------------------------------------------------------------------------
# Data classes  (identical structure to find_group_private_vcf_regions.py)
# ---------------------------------------------------------------------------

@dataclass
class SiteHit:
    chrom: str
    ref_pos: int
    ref_allele: str
    target_allele: str
    classification: str             # snp or indel
    n_targets_observed: int
    n_nontarget_observed: int
    n_nontarget_carries: int
    n_nontarget_distinct: int
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
    classification: str
    targets: str
    target_alleles: str
    ref_alleles: str
    nontarget_alleles: str
    n_targets_expected: int
    n_targets_observed: int
    min_nontarget_observed: int
    max_nontarget_observed: int
    max_nontarget_carries: int


@dataclass
class ScanStats:
    total_variant_positions: int = 0
    sites_skipped_no_target_coverage: int = 0
    sites_skipped_targets_disagree: int = 0
    sites_skipped_no_nontarget: int = 0
    sites_skipped_all_nontarget_carry: int = 0
    qualifying_sites: int = 0
    qualifying_intervals: int = 0


# ---------------------------------------------------------------------------
# BAM parsing
# ---------------------------------------------------------------------------

def _collect_variants_and_coverage(
    bam_path: str,
    ref_fasta: pysam.FastaFile,
    chrom: Optional[str],
    start: Optional[int],
    end: Optional[int],
) -> Tuple[Dict[Tuple[str, int], Tuple[str, str]], List[Tuple[str, int, int]]]:
    """Single pass through one BAM.

    Returns:
        variants  – {(chrom, pos_1based): (ref_allele, query_allele)}
                    query_allele is DELETION for deleted bases
        intervals – [(chrom, start_0based, end_0based)] one per primary alignment,
                    sorted by (chrom, start)
    """
    variants: Dict[Tuple[str, int], Tuple[str, str]] = {}
    intervals: List[Tuple[str, int, int]] = []

    bam = pysam.AlignmentFile(bam_path, "rb")

    fetch_kwargs: dict = {}
    if chrom is not None:
        fetch_kwargs["contig"] = chrom
        if start is not None:
            fetch_kwargs["start"] = start - 1   # pysam fetch uses 0-based
        if end is not None:
            fetch_kwargs["stop"] = end

    for read in bam.fetch(**fetch_kwargs):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        ref_name = read.reference_name
        ref_start = read.reference_start    # 0-based inclusive
        ref_end = read.reference_end        # 0-based exclusive

        intervals.append((ref_name, ref_start, ref_end))

        query_seq = read.query_sequence
        if query_seq is None:
            continue

        # Fetch the reference slice once for this alignment
        try:
            ref_seq = ref_fasta.fetch(ref_name, ref_start, ref_end).upper()
        except (KeyError, ValueError):
            continue

        # Iterate over all aligned pairs (including gaps)
        for qpos, rpos in read.get_aligned_pairs(matches_only=False):
            if rpos is None:
                continue    # insertion — no reference coordinate, skip

            pos_1based = rpos + 1

            if qpos is None:
                # Deletion: query has no base at this reference position
                ref_base = ref_seq[rpos - ref_start] if rpos - ref_start < len(ref_seq) else "N"
                if ref_base not in "ACGT":
                    continue
                key = (ref_name, pos_1based)
                if key not in variants:
                    variants[key] = (ref_base, DELETION)
            else:
                # Match or mismatch
                ref_base = ref_seq[rpos - ref_start] if rpos - ref_start < len(ref_seq) else "N"
                query_base = query_seq[qpos].upper()
                if ref_base not in "ACGT" or query_base not in "ACGT":
                    continue
                if query_base != ref_base:
                    key = (ref_name, pos_1based)
                    if key not in variants:
                        variants[key] = (ref_base, query_base)

    intervals.sort()
    return variants, intervals


def _build_coverage_index(
    intervals: List[Tuple[str, int, int]],
) -> Dict[str, Tuple[List[int], List[int]]]:
    """Build per-chromosome sorted arrays of (starts, ends) for fast lookup."""
    by_chrom: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    for chrom, s, e in intervals:
        by_chrom[chrom].append((s, e))
    index = {}
    for chrom, ivs in by_chrom.items():
        ivs.sort()
        index[chrom] = ([iv[0] for iv in ivs], [iv[1] for iv in ivs])
    return index


def _is_covered(
    index: Dict[str, Tuple[List[int], List[int]]],
    chrom: str,
    pos_1based: int,
) -> bool:
    """O(log n) coverage check using binary search on sorted start positions."""
    if chrom not in index:
        return False
    pos_0based = pos_1based - 1
    starts, ends = index[chrom]
    i = bisect.bisect_right(starts, pos_0based) - 1
    return i >= 0 and ends[i] > pos_0based


# ---------------------------------------------------------------------------
# Main scan logic
# ---------------------------------------------------------------------------

def scan_bams(
    sample_bams: Dict[str, str],
    targets: Sequence[str],
    ref_path: str,
    chrom: Optional[str],
    start: Optional[int],
    end: Optional[int],
) -> Tuple[List[SiteHit], ScanStats, List[str]]:

    all_samples = list(sample_bams.keys())
    other_samples = [s for s in all_samples if s not in targets]
    ref_fasta = pysam.FastaFile(ref_path)

    # ------------------------------------------------------------------
    # Phase 1 – one pass per BAM: collect variants and coverage intervals
    # ------------------------------------------------------------------
    print("Phase 1: extracting variants from BAM files...")
    sample_variants: Dict[str, Dict[Tuple[str, int], Tuple[str, str]]] = {}
    sample_cov_index: Dict[str, Dict[str, Tuple[List[int], List[int]]]] = {}

    for sample, bam_path in sample_bams.items():
        print(f"  {sample} ({bam_path})")
        variants, intervals = _collect_variants_and_coverage(
            bam_path, ref_fasta, chrom, start, end
        )
        sample_variants[sample] = variants
        sample_cov_index[sample] = _build_coverage_index(intervals)

    # ------------------------------------------------------------------
    # Phase 2 – collect union of all variant positions
    # ------------------------------------------------------------------
    all_positions: Set[Tuple[str, int]] = set()
    for sv in sample_variants.values():
        all_positions.update(sv.keys())

    print(f"Phase 2: testing {len(all_positions):,} variant positions...")

    # ------------------------------------------------------------------
    # Phase 3 – evaluate each position
    # ------------------------------------------------------------------
    stats = ScanStats()
    hits: List[SiteHit] = []

    for pos_chrom, pos in sorted(all_positions):
        stats.total_variant_positions += 1

        # ---- target alleles ----
        target_alleles: List[str] = []
        ref_allele: Optional[str] = None

        for sid in targets:
            if (pos_chrom, pos) in sample_variants[sid]:
                ref_b, alt_b = sample_variants[sid][(pos_chrom, pos)]
                target_alleles.append(alt_b)
                if ref_allele is None:
                    ref_allele = ref_b
            elif _is_covered(sample_cov_index[sid], pos_chrom, pos):
                if ref_allele is None:
                    ref_allele = ref_fasta.fetch(pos_chrom, pos - 1, pos).upper()
                target_alleles.append(ref_allele)
            # else: not covered — missing, excluded from count

        if not target_alleles:
            stats.sites_skipped_no_target_coverage += 1
            continue

        if len(set(target_alleles)) != 1:
            stats.sites_skipped_targets_disagree += 1
            continue

        target_allele = target_alleles[0]
        n_targets_observed = len(target_alleles)
        all_targets_observed = (n_targets_observed == len(targets))

        if ref_allele is None:
            ref_allele = ref_fasta.fetch(pos_chrom, pos - 1, pos).upper()

        # ---- non-target alleles ----
        nontarget_observed: List[str] = []
        nontarget_carries: List[str] = []
        nontarget_allele_set: Set[str] = set()

        for sid in other_samples:
            if (pos_chrom, pos) in sample_variants[sid]:
                _, alt_b = sample_variants[sid][(pos_chrom, pos)]
                nontarget_observed.append(sid)
                nontarget_allele_set.add(alt_b)
                if alt_b == target_allele:
                    nontarget_carries.append(sid)
            elif _is_covered(sample_cov_index[sid], pos_chrom, pos):
                nontarget_observed.append(sid)
                nontarget_allele_set.add(ref_allele)
                if ref_allele == target_allele:
                    nontarget_carries.append(sid)

        if not nontarget_observed:
            stats.sites_skipped_no_nontarget += 1
            continue

        if len(nontarget_carries) == len(nontarget_observed):
            stats.sites_skipped_all_nontarget_carry += 1
            continue

        # ---- classify strictness ----
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
        classification = "indel" if target_allele == DELETION or ref_allele != target_allele and len(ref_allele) != len(target_allele) else "snp"

        hits.append(SiteHit(
            chrom=pos_chrom,
            ref_pos=pos,
            ref_allele=ref_allele,
            target_allele=target_allele,
            classification=classification,
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
# Interval merging  (identical to VCF script)
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
# Output  (identical column layout to VCF script)
# ---------------------------------------------------------------------------

def write_intervals_tsv(out_path: Path, intervals: Sequence[IntervalHit]) -> None:
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
    samples_file: str,
    ref_path: str,
    targets: Sequence[str],
    all_samples: Sequence[str],
    intervals: Sequence[IntervalHit],
    chrom: Optional[str],
    start: Optional[int],
    end: Optional[int],
) -> None:
    from collections import Counter
    level_counts = Counter(i.strictness_level for i in intervals)
    class_counts = Counter(i.classification for i in intervals)
    rows = [
        ("samples_file",                        samples_file),
        ("ref",                                 ref_path),
        ("targets",                             ",".join(targets)),
        ("all_samples",                         ",".join(all_samples)),
        ("n_targets",                           str(len(targets))),
        ("n_nontargets",                        str(len(all_samples) - len(targets))),
        ("region_chrom",                        chrom or ""),
        ("window_start",                        "" if start is None else str(start)),
        ("window_end",                          "" if end is None else str(end)),
        ("total_variant_positions",             str(stats.total_variant_positions)),
        ("sites_skipped_no_target_coverage",    str(stats.sites_skipped_no_target_coverage)),
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

def load_samples_file(path: str) -> Dict[str, str]:
    """Read a two-column TSV: sample_name <tab> bam_path"""
    samples: Dict[str, str] = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                raise SystemExit(
                    f"ERROR: samples file must have two tab-separated columns "
                    f"(sample_name, bam_path). Bad line: {line!r}"
                )
            sample, bam = parts[0].strip(), parts[1].strip()
            if sample in samples:
                raise SystemExit(f"ERROR: duplicate sample name in samples file: {sample}")
            samples[sample] = bam
    return samples


def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Scan assembly BAM files directly for group-private variant positions. "
            "Memory-efficient alternative to bcftools mpileup for long contigs."
        )
    )
    p.add_argument(
        "--samples", required=True,
        help=(
            "Two-column tab-separated file (no header): "
            "sample_name<TAB>bam_path. One sample per line."
        ),
    )
    p.add_argument(
        "--ref", required=True,
        help="Reference FASTA indexed with samtools faidx",
    )
    p.add_argument(
        "--targets", required=True, nargs="+",
        help="Sample names of the target group (must match first column of --samples)",
    )
    p.add_argument(
        "--out-prefix", required=True,
        help="Output path prefix; writes <prefix>.intervals.tsv and <prefix>.summary.tsv",
    )
    p.add_argument(
        "--region", default=None,
        help="Chromosome/contig name to restrict scanning, e.g. Gm15",
    )
    p.add_argument(
        "--start", type=int, default=None,
        help="1-based start position filter (used with --region)",
    )
    p.add_argument(
        "--end", type=int, default=None,
        help="1-based end position filter (used with --region)",
    )
    return p


def main() -> None:
    args = build_argparser().parse_args()

    if args.start is not None and args.end is not None and args.start > args.end:
        raise SystemExit("ERROR: --start must be <= --end")
    if len(set(args.targets)) != len(args.targets):
        raise SystemExit("ERROR: --targets contains duplicate sample names")

    sample_bams = load_samples_file(args.samples)

    missing = [s for s in args.targets if s not in sample_bams]
    if missing:
        raise SystemExit(
            f"ERROR: target samples not found in samples file: {', '.join(missing)}\n"
            f"Available: {', '.join(sample_bams)}"
        )

    hits, stats, all_samples = scan_bams(
        sample_bams=sample_bams,
        targets=args.targets,
        ref_path=args.ref,
        chrom=args.region,
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
        summary_path, stats, args.samples, args.ref,
        args.targets, all_samples, intervals,
        args.region, args.start, args.end,
    )

    print(f"Samples       : {', '.join(all_samples)}")
    print(f"Targets       : {', '.join(args.targets)}")
    print(f"Variant sites : {stats.total_variant_positions:,}")
    print(f"Qualifying sites    : {stats.qualifying_sites:,}")
    print(f"Qualifying intervals: {len(intervals):,}")
    print(f"Wrote: {intervals_path}")
    print(f"Wrote: {summary_path}")


if __name__ == "__main__":
    main()
