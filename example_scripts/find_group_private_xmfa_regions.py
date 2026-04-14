#!/usr/bin/env python3
"""Scan an XMFA for strict private-state regions for a user-defined target group.

A position qualifies when:
  1. all target genomes present in the block share the exact same state
  2. the target state is one of A/C/G/T/-
  3. no non-target genome with an unambiguous observed state matches the target state
  4. at least one non-target genome with an unambiguous observed state is present

Ambiguous states (e.g. N or other IUPAC codes) are ignored.
Shared gaps are allowed and treated as indel-like private states.

Coordinates are reported in the reference genome supplied by --ref-sid.
Only alignment columns anchored to a reference base are reportable. This means
columns where the reference itself has a gap are not emitted because they do not
have a reference coordinate.
"""

from __future__ import annotations

import argparse, gzip, re, os
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Set, Tuple


VALID_STATES = {"A", "C", "G", "T", "-"}
AMBIGUOUS_STATES = {"N", "?", "X"}

# Strictness levels (1 = most strict, higher = more relaxed)
# 1: all targets present, no non-target shares the allele
# 2: all targets present, but ≥1 non-target also carries the allele
# 3: ≥1 target absent from block, but no non-target shares the allele
# 4: ≥1 target absent from block AND ≥1 non-target shares the allele
STRICTNESS_LEVELS: Dict[int, str] = {
    1: "strict",
    2: "nontarget_shared",
    3: "target_partial",
    4: "target_partial_nontarget_shared",
}
XMFA_HEADER_RE = re.compile(
    r"^>\s*(?P<idx>\S+):(?P<start>\d+)-(?P<end>\d+)\s+(?P<strand>[+-])\s+(?P<label>\S+)$"
)


def open_text(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def normalize_sid(label: str) -> str:
    return os.path.basename(label).split(".")[0]


@dataclass
class SiteHit:
    ref_pos: int
    block_index: int
    target_allele: str
    classification: str  # snp or indel
    n_nontarget_observed: int
    n_nontarget_distinct: int
    nontarget_alleles: Tuple[str, ...]
    n_targets_observed: int
    strictness: str
    strictness_level: int


@dataclass
class IntervalHit:
    interval_id: str
    start: int
    end: int
    span_bp: int
    n_sites: int
    n_blocks: int
    classification: str
    target_alleles: str
    nontarget_alleles: str
    n_targets_expected: int
    n_targets_observed: int
    min_nontarget_observed: int
    max_nontarget_observed: int
    block_ids: str
    strictness: str
    strictness_level: int


@dataclass
class ScanStats:
    total_blocks: int = 0
    blocks_with_ref: int = 0
    blocks_scanned: int = 0
    blocks_skipped_no_ref: int = 0
    blocks_partial_targets: int = 0
    blocks_outside_window: int = 0
    alignment_columns_seen: int = 0
    ref_anchored_columns_seen: int = 0
    ref_anchored_columns_in_window: int = 0
    columns_with_all_targets_present: int = 0
    columns_with_valid_target_state: int = 0
    informative_columns_tested: int = 0
    qualifying_columns: int = 0
    qualifying_intervals: int = 0


def read_xmfa_blocks(xmfa_path: str) -> Iterator[Dict[str, dict]]:
    block: Dict[str, dict] = {}
    cur_sid: Optional[str] = None

    def flush() -> Optional[Dict[str, dict]]:
        nonlocal block
        if not block:
            return None
        out: Dict[str, dict] = {}
        for sid, rec in block.items():
            out[sid] = {
                "start": rec["start"],
                "end": rec["end"],
                "strand": rec["strand"],
                "label": rec["label"],
                "seq": "".join(rec["seq_chunks"]),
            }
        block = {}
        return out

    with open_text(xmfa_path) as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue

            if line == "=":
                emitted = flush()
                if emitted is not None:
                    yield emitted
                cur_sid = None
                continue

            if line.startswith(">"):
                match = XMFA_HEADER_RE.match(line)
                if not match:
                    raise ValueError(f"Unrecognized XMFA header: {line}")
                label = match.group("label")
                sid = normalize_sid(label)
                block[sid] = {
                    "start": int(match.group("start")),
                    "end": int(match.group("end")),
                    "strand": match.group("strand"),
                    "label": label,
                    "seq_chunks": [],
                }
                cur_sid = sid
                continue

            if cur_sid is not None:
                block[cur_sid]["seq_chunks"].append(line.strip())

    emitted = flush()
    if emitted is not None:
        yield emitted


def discover_sids(xmfa_path: str) -> List[str]:
    sids: Set[str] = set()
    for block in read_xmfa_blocks(xmfa_path):
        sids.update(block.keys())
    return sorted(sids)


def normalize_state(base: str) -> str:
    b = base.upper()
    if b == ".":
        return "-"
    return b


def classify_site(target_allele: str, nontarget_alleles: Set[str]) -> str:
    if target_allele == "-" or "-" in nontarget_alleles:
        return "indel"
    return "snp"


def scan_block(
    block: Dict[str, dict],
    block_index: int,
    ref_sid: str,
    targets: Sequence[str],
    start: Optional[int],
    end: Optional[int],
    stats: ScanStats,
    ref_offset: int = 0,
) -> List[SiteHit]:
    hits: List[SiteHit] = []
    stats.total_blocks += 1

    if ref_sid not in block:
        stats.blocks_skipped_no_ref += 1
        return hits
    stats.blocks_with_ref += 1

    targets_in_block = [sid for sid in targets if sid in block]
    if not targets_in_block:
        return hits
    if len(targets_in_block) < len(targets):
        stats.blocks_partial_targets += 1

    ref_seq = block[ref_sid]["seq"]
    block_len = len(ref_seq)
    stats.blocks_scanned += 1
    stats.alignment_columns_seen += block_len

    seqs = {sid: rec["seq"] for sid, rec in block.items()}
    seq_lens = {len(seq) for seq in seqs.values()}
    if len(seq_lens) != 1:
        raise ValueError(
            f"Block {block_index} has inconsistent sequence lengths: {sorted(seq_lens)}"
        )

    ref_start = block[ref_sid]["start"]
    ref_end = block[ref_sid]["end"]
    ref_lo, ref_hi = (ref_start, ref_end) if ref_start <= ref_end else (ref_end, ref_start)
    if start is not None and ref_hi < start:
        stats.blocks_outside_window += 1
        return hits
    if end is not None and ref_lo > end:
        stats.blocks_outside_window += 1
        return hits

    other_sids = [sid for sid in block.keys() if sid not in targets]
    ref_strand = block[ref_sid]["strand"]
    ref_pos = ref_lo if ref_strand == "+" else ref_hi
    ref_step = 1 if ref_strand == "+" else -1

    for col in range(block_len):
        ref_state = normalize_state(ref_seq[col])
        if ref_state not in {"A", "C", "G", "T"}:
            continue

        stats.ref_anchored_columns_seen += 1
        cur_pos = ref_pos + ref_offset
        ref_pos += ref_step

        if start is not None and cur_pos < start:
            continue
        if end is not None and cur_pos > end:
            continue
        stats.ref_anchored_columns_in_window += 1

        target_states = [normalize_state(seqs[sid][col]) for sid in targets_in_block]
        target_valid = [s for s in target_states if s in VALID_STATES]
        if len(target_valid) == 0:
            continue
        if len(set(target_valid)) != 1:
            continue

        n_targets_observed = len(target_valid)
        all_targets_observed = (n_targets_observed == len(targets))
        target_allele = target_valid[0]
        stats.columns_with_valid_target_state += 1
        if all_targets_observed:
            stats.columns_with_all_targets_present += 1

        nontarget_valid_states: List[str] = []
        for sid in other_sids:
            state = normalize_state(seqs[sid][col])
            if state in VALID_STATES:
                nontarget_valid_states.append(state)

        if not nontarget_valid_states:
            continue

        nontarget_unique = set(nontarget_valid_states)
        # Skip positions where every non-target also shares the target allele (not private at all)
        if nontarget_unique == {target_allele}:
            continue

        stats.informative_columns_tested += 1
        strict_nontarget = target_allele not in nontarget_unique
        if all_targets_observed and strict_nontarget:
            strictness_level = 1
        elif all_targets_observed and not strict_nontarget:
            strictness_level = 2
        elif not all_targets_observed and strict_nontarget:
            strictness_level = 3
        else:
            strictness_level = 4
        strictness = STRICTNESS_LEVELS[strictness_level]

        classification = classify_site(target_allele, nontarget_unique)
        hits.append(
            SiteHit(
                ref_pos=cur_pos,
                block_index=block_index,
                target_allele=target_allele,
                classification=classification,
                n_nontarget_observed=len(nontarget_valid_states),
                n_nontarget_distinct=len(nontarget_unique),
                nontarget_alleles=tuple(sorted(nontarget_unique)),
                n_targets_observed=n_targets_observed,
                strictness=strictness,
                strictness_level=strictness_level,
            )
        )
        stats.qualifying_columns += 1

    return hits


def merge_sites_to_intervals(
    sites: Sequence[SiteHit],
    n_targets_expected: int,
) -> List[IntervalHit]:
    if not sites:
        return []

    sorted_sites = sorted(sites, key=lambda x: x.ref_pos)
    intervals: List[IntervalHit] = []
    current: List[SiteHit] = [sorted_sites[0]]

    def finalize(chunk: Sequence[SiteHit], interval_num: int) -> IntervalHit:
        start = chunk[0].ref_pos
        end = chunk[-1].ref_pos
        site_classes = sorted({site.classification for site in chunk})
        if len(site_classes) == 1:
            interval_class = site_classes[0]
        else:
            interval_class = "mixed"
        target_alleles = ",".join(sorted({site.target_allele for site in chunk}))
        nontarget_alleles = ",".join(
            sorted({allele for site in chunk for allele in site.nontarget_alleles})
        )
        block_ids = sorted({str(site.block_index) for site in chunk}, key=int)
        interval_strictness_level = max(s.strictness_level for s in chunk)
        interval_strictness = STRICTNESS_LEVELS[interval_strictness_level]
        return IntervalHit(
            interval_id=f"region_{interval_num:06d}",
            start=start,
            end=end,
            span_bp=end - start + 1,
            n_sites=len(chunk),
            n_blocks=len(block_ids),
            classification=interval_class,
            target_alleles=target_alleles,
            nontarget_alleles=nontarget_alleles,
            n_targets_expected=n_targets_expected,
            n_targets_observed=min(s.n_targets_observed for s in chunk),
            min_nontarget_observed=min(site.n_nontarget_observed for site in chunk),
            max_nontarget_observed=max(site.n_nontarget_observed for site in chunk),
            block_ids=",".join(block_ids),
            strictness=interval_strictness,
            strictness_level=interval_strictness_level,
        )

    interval_num = 1
    for site in sorted_sites[1:]:
        if site.ref_pos == current[-1].ref_pos + 1:
            current.append(site)
        else:
            intervals.append(finalize(current, interval_num))
            interval_num += 1
            current = [site]
    intervals.append(finalize(current, interval_num))
    return intervals


def write_intervals_tsv(
    out_path: Path,
    intervals: Sequence[IntervalHit],
    targets: Sequence[str],
) -> None:
    with out_path.open("w") as out:
        header = [
            "interval_id",
            "strictness",
            "strictness_level",
            "start",
            "end",
            "span_bp",
            "n_sites",
            "n_blocks",
            "classification",
            "targets",
            "target_alleles",
            "nontarget_alleles",
            "n_targets_expected",
            "n_targets_observed",
            "min_nontarget_observed",
            "max_nontarget_observed",
            "block_ids",
        ]
        out.write("\t".join(header) + "\n")
        target_str = ",".join(targets)
        for row in intervals:
            vals = [
                row.interval_id,
                row.strictness,
                str(row.strictness_level),
                str(row.start),
                str(row.end),
                str(row.span_bp),
                str(row.n_sites),
                str(row.n_blocks),
                row.classification,
                target_str,
                row.target_alleles,
                row.nontarget_alleles,
                str(row.n_targets_expected),
                str(row.n_targets_observed),
                str(row.min_nontarget_observed),
                str(row.max_nontarget_observed),
                row.block_ids,
            ]
            out.write("\t".join(vals) + "\n")


def write_summary_tsv(
    out_path: Path,
    stats: ScanStats,
    xmfa_path: str,
    ref_sid: str,
    targets: Sequence[str],
    all_sids: Sequence[str],
    intervals: Sequence[IntervalHit],
    start: Optional[int],
    end: Optional[int],
    ref_offset: int = 0,
) -> None:
    class_counts = Counter(interval.classification for interval in intervals)
    total_interval_bp = sum(interval.span_bp for interval in intervals)
    summary_rows = [
        ("xmfa", xmfa_path),
        ("ref_sid", ref_sid),
        ("targets", ",".join(targets)),
        ("all_genomes_discovered", ",".join(all_sids)),
        ("ref_offset", str(ref_offset)),
        ("window_start", "" if start is None else str(start)),
        ("window_end", "" if end is None else str(end)),
        ("total_blocks", str(stats.total_blocks)),
        ("blocks_with_ref", str(stats.blocks_with_ref)),
        ("blocks_scanned", str(stats.blocks_scanned)),
        ("blocks_skipped_no_ref", str(stats.blocks_skipped_no_ref)),
        ("blocks_partial_targets", str(stats.blocks_partial_targets)),
        ("blocks_outside_window", str(stats.blocks_outside_window)),
        ("alignment_columns_seen", str(stats.alignment_columns_seen)),
        ("ref_anchored_columns_seen", str(stats.ref_anchored_columns_seen)),
        ("ref_anchored_columns_in_window", str(stats.ref_anchored_columns_in_window)),
        ("columns_with_all_targets_present", str(stats.columns_with_all_targets_present)),
        ("columns_with_valid_target_state", str(stats.columns_with_valid_target_state)),
        ("informative_columns_tested", str(stats.informative_columns_tested)),
        ("qualifying_columns", str(stats.qualifying_columns)),
        ("qualifying_intervals", str(len(intervals))),
        ("qualifying_interval_bp_total", str(total_interval_bp)),
        ("interval_class_snp", str(class_counts.get("snp", 0))),
        ("interval_class_indel", str(class_counts.get("indel", 0))),
        ("interval_class_mixed", str(class_counts.get("mixed", 0))),
    ]
    with out_path.open("w") as out:
        out.write("metric\tvalue\n")
        for metric, value in summary_rows:
            out.write(f"{metric}\t{value}\n")


def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Scan an XMFA for strict private-state regions where all target genomes share "
            "the same exact allele and every observed non-target genome differs. "
            "Shared gaps are allowed. Ambiguous states are ignored."
        )
    )
    p.add_argument("--xmfa", required=True, help="Input XMFA file (.xmfa or .xmfa.gz)")
    p.add_argument("--ref-sid", required=True, help="Reference genome ID after XMFA normalization")
    p.add_argument(
        "--targets",
        required=True,
        nargs="+",
        help="One or more target genomes to test as a strict private-state group",
    )
    p.add_argument(
        "--out-prefix",
        required=True,
        help="Output prefix. Writes <prefix>.intervals.tsv and <prefix>.summary.tsv",
    )
    p.add_argument(
        "--start",
        type=int,
        default=None,
        help="Optional reference start coordinate for scanning a sub-window",
    )
    p.add_argument(
        "--end",
        type=int,
        default=None,
        help="Optional reference end coordinate for scanning a sub-window",
    )
    p.add_argument(
        "--ref-offset",
        type=int,
        default=0,
        help=(
            "Offset added to every reported reference coordinate. Use this when the XMFA "
            "was built from an extracted sub-region of the chromosome so that output "
            "positions are in chromosomal coordinates (e.g. --ref-offset 10000000 if the "
            "extract started at chr15:10000001)."
        ),
    )
    return p


def main() -> None:
    args = build_argparser().parse_args()
    if args.start is not None and args.end is not None and args.start > args.end:
        raise SystemExit("ERROR: --start must be <= --end")

    all_sids = discover_sids(args.xmfa)
    missing_targets = [sid for sid in args.targets if sid not in all_sids]
    if args.ref_sid not in all_sids:
        raise SystemExit(
            "ERROR: --ref-sid not found in XMFA after normalization. "
            f"Found: {', '.join(all_sids)}"
        )
    if missing_targets:
        raise SystemExit(
            "ERROR: One or more --targets were not found in the XMFA after normalization: "
            f"{', '.join(missing_targets)}\nFound: {', '.join(all_sids)}"
        )
    if len(set(args.targets)) != len(args.targets):
        raise SystemExit("ERROR: --targets contains duplicate genome IDs")

    stats = ScanStats()
    all_hits: List[SiteHit] = []
    for block_index, block in enumerate(read_xmfa_blocks(args.xmfa), start=1):
        all_hits.extend(
            scan_block(
                block=block,
                block_index=block_index,
                ref_sid=args.ref_sid,
                targets=args.targets,
                start=args.start,
                end=args.end,
                stats=stats,
                ref_offset=args.ref_offset,
            )
        )

    intervals = merge_sites_to_intervals(all_hits, n_targets_expected=len(args.targets))
    stats.qualifying_intervals = len(intervals)

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    intervals_path = Path(f"{args.out_prefix}.intervals.tsv")
    summary_path = Path(f"{args.out_prefix}.summary.tsv")

    write_intervals_tsv(intervals_path, intervals, args.targets)
    write_summary_tsv(
        summary_path,
        stats,
        args.xmfa,
        args.ref_sid,
        args.targets,
        all_sids,
        intervals,
        args.start,
        args.end,
        args.ref_offset,
    )

    print(f"Discovered genomes: {', '.join(all_sids)}")
    print(f"Targets: {', '.join(args.targets)}")
    print(f"Qualifying columns: {stats.qualifying_columns}")
    print(f"Qualifying intervals: {len(intervals)}")
    print(f"Wrote: {intervals_path}")
    print(f"Wrote: {summary_path}")


if __name__ == "__main__":
    main()
