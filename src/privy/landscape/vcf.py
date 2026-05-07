"""VCF sliding-window landscape analysis.

This module computes Panex-native window metrics from multisample VCF records.
It is deliberately focused on target/off-target-aware interpretation rather
than replacing mature population-genetic tools.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from statistics import mean, median
from typing import Literal

from privy.io.vcf import (
    Genotype,
    VariantRecordLike,
    get_vcf_samples,
    has_alt_allele,
    is_missing_genotype,
    stream_vcf_records,
)
from privy.pangenome.model import PangenomeGroups, resolve_pangenome_groups

WindowMode = Literal["records", "bp"]

LANDSCAPE_SAMPLE_WINDOW_COLUMNS = [
    "window_id",
    "contig",
    "window_index",
    "start",
    "end",
    "midpoint",
    "window_mode",
    "n_variants",
    "sample",
    "cohort_role",
    "called_n",
    "missing_n",
    "missing_rate",
    "het_n",
    "het_rate",
    "nonref_n",
    "nonref_rate",
    "minor_genotype_n",
    "minor_genotype_rate",
    "rare_alt_n",
    "rare_alt_rate",
    "private_alt_n",
    "private_alt_rate",
    "median_call_freq",
    "nearest_background",
    "nearest_background_role",
    "nearest_similarity",
    "similarity_compared_variants",
]

LANDSCAPE_WINDOW_COLUMNS = [
    "window_id",
    "contig",
    "window_index",
    "start",
    "end",
    "midpoint",
    "window_mode",
    "n_variants",
    "span_bp",
    "density_variants_per_kb",
    "target_mean_missing_rate",
    "offtarget_mean_missing_rate",
    "target_mean_nonref_rate",
    "offtarget_mean_nonref_rate",
    "target_private_alt_n",
    "offtarget_private_alt_n",
    "target_private_alt_rate",
    "offtarget_private_alt_rate",
    "top_nearest_background",
    "top_nearest_background_n",
]

BACKGROUND_BLOCK_COLUMNS = [
    "block_id",
    "sample",
    "cohort_role",
    "contig",
    "start",
    "end",
    "n_windows",
    "nearest_background",
    "nearest_background_role",
    "mean_similarity",
]

CANDIDATE_INTROGRESSION_BLOCK_COLUMNS = [
    "block_id",
    "sample",
    "contig",
    "start",
    "end",
    "n_windows",
    "candidate_donor",
    "candidate_donor_role",
    "mean_donor_similarity",
    "mean_nearest_target_similarity",
    "mean_similarity_delta",
    "max_missing_rate",
    "mean_private_alt_rate",
    "mean_nonref_rate",
    "evidence_class",
    "interpretation",
]

LANDSCAPE_SIMILARITY_COLUMNS = [
    "window_id",
    "contig",
    "window_index",
    "start",
    "end",
    "sample_a",
    "sample_b",
    "similarity",
    "compared_variants",
]


@dataclass(frozen=True)
class LandscapeResult:
    """Rows produced by one VCF landscape run."""

    samples: tuple[str, ...]
    groups: PangenomeGroups
    window_mode: WindowMode
    sample_rows: list[dict[str, object]]
    window_rows: list[dict[str, object]]
    background_block_rows: list[dict[str, object]]
    candidate_introgression_rows: list[dict[str, object]]
    similarity_rows: list[dict[str, object]]


@dataclass(frozen=True)
class _VariantSnapshot:
    contig: str
    start: int
    end: int
    genotypes: tuple[Genotype | None, ...]
    alt_carriers: tuple[frozenset[int], ...]


@dataclass(frozen=True)
class _Window:
    window_id: str
    contig: str
    window_index: int
    start: int
    end: int
    mode: WindowMode
    variants: tuple[_VariantSnapshot, ...]


@dataclass(frozen=True)
class _PairwiseStats:
    similarity: float | None
    compared: int


def run_vcf_landscape(
    vcf_path: Path,
    targets: list[str],
    off_targets: list[str] | None = None,
    ignored_samples: list[str] | None = None,
    window_records: int = 200,
    step_records: int = 50,
    window_bp: int | None = None,
    step_bp: int | None = None,
    pass_only: bool = True,
    min_qual: float | None = None,
    rare_max_count: int = 1,
    rare_max_freq: float = 0.05,
    min_called_for_freq: int = 10,
    min_freq_values: int = 10,
    min_background_similarity: float = 0.65,
    min_introgression_similarity: float | None = None,
    min_introgression_delta: float = 0.0,
    max_introgression_missing_rate: float = 0.5,
    min_introgression_windows: int = 1,
) -> LandscapeResult:
    """Run sliding-window VCF landscape analysis.

    Coordinates in output rows are 0-based, half-open.
    """
    if not vcf_path.exists():
        raise FileNotFoundError(f"VCF file not found: {vcf_path}")
    if window_bp is not None and window_bp <= 0:
        raise ValueError("window_bp must be positive.")
    if step_bp is not None and step_bp <= 0:
        raise ValueError("step_bp must be positive.")
    if window_records <= 0:
        raise ValueError("window_records must be positive.")
    if step_records <= 0:
        raise ValueError("step_records must be positive.")
    if rare_max_count < 0:
        raise ValueError("rare_max_count must be non-negative.")
    if rare_max_freq < 0 or rare_max_freq > 1:
        raise ValueError("rare_max_freq must be between 0 and 1.")
    if min_called_for_freq < 0:
        raise ValueError("min_called_for_freq must be non-negative.")
    if min_freq_values < 0:
        raise ValueError("min_freq_values must be non-negative.")
    if min_background_similarity < 0 or min_background_similarity > 1:
        raise ValueError("min_background_similarity must be between 0 and 1.")
    if (
        min_introgression_similarity is not None
        and (min_introgression_similarity < 0 or min_introgression_similarity > 1)
    ):
        raise ValueError("min_introgression_similarity must be between 0 and 1.")
    if min_introgression_delta < 0 or min_introgression_delta > 1:
        raise ValueError("min_introgression_delta must be between 0 and 1.")
    if max_introgression_missing_rate < 0 or max_introgression_missing_rate > 1:
        raise ValueError("max_introgression_missing_rate must be between 0 and 1.")
    if min_introgression_windows <= 0:
        raise ValueError("min_introgression_windows must be positive.")

    samples = tuple(get_vcf_samples(vcf_path))
    groups = resolve_pangenome_groups(
        all_samples=samples,
        targets=targets,
        off_targets=off_targets,
        ignored_samples=ignored_samples,
    )
    active_sample_indices = tuple(i for i, sample in enumerate(samples) if sample in groups.full)
    sample_to_role = _sample_roles(groups)

    mode: WindowMode = "bp" if window_bp is not None else "records"
    effective_step_bp = step_bp or window_bp

    sample_rows: list[dict[str, object]] = []
    window_rows: list[dict[str, object]] = []
    similarity_rows: list[dict[str, object]] = []
    background_rows_input: list[dict[str, object]] = []

    contig_records: list[_VariantSnapshot] = []
    current_contig: str | None = None
    window_counter = 0

    for record in stream_vcf_records(vcf_path):
        if record.alts is None:
            continue
        if pass_only and not _record_passes(record_filter_values=list(record.filter)):
            continue
        if min_qual is not None and (record.qual is None or record.qual < min_qual):
            continue

        if current_contig is not None and record.chrom != current_contig:
            window_counter = _analyze_contig_windows(
                records=tuple(contig_records),
                mode=mode,
                window_counter_start=window_counter,
                samples=samples,
                active_sample_indices=active_sample_indices,
                sample_to_role=sample_to_role,
                groups=groups,
                sample_rows=sample_rows,
                window_rows=window_rows,
                similarity_rows=similarity_rows,
                background_rows_input=background_rows_input,
                window_records=window_records,
                step_records=step_records,
                window_bp=window_bp,
                step_bp=effective_step_bp,
                rare_max_count=rare_max_count,
                rare_max_freq=rare_max_freq,
                min_called_for_freq=min_called_for_freq,
                min_freq_values=min_freq_values,
            )
            contig_records = []

        current_contig = record.chrom
        contig_records.append(_snapshot_record(record, samples))

    if contig_records:
        window_counter = _analyze_contig_windows(
            records=tuple(contig_records),
            mode=mode,
            window_counter_start=window_counter,
            samples=samples,
            active_sample_indices=active_sample_indices,
            sample_to_role=sample_to_role,
            groups=groups,
            sample_rows=sample_rows,
            window_rows=window_rows,
            similarity_rows=similarity_rows,
            background_rows_input=background_rows_input,
            window_records=window_records,
            step_records=step_records,
            window_bp=window_bp,
            step_bp=effective_step_bp,
            rare_max_count=rare_max_count,
            rare_max_freq=rare_max_freq,
            min_called_for_freq=min_called_for_freq,
            min_freq_values=min_freq_values,
        )

    del window_counter
    block_rows = build_background_blocks(
        background_rows_input,
        min_background_similarity=min_background_similarity,
    )
    introgression_rows = build_candidate_introgression_blocks(
        sample_window_rows=background_rows_input,
        similarity_rows=similarity_rows,
        target_samples=groups.target,
        min_introgression_similarity=(
            min_background_similarity
            if min_introgression_similarity is None
            else min_introgression_similarity
        ),
        min_introgression_delta=min_introgression_delta,
        max_introgression_missing_rate=max_introgression_missing_rate,
        min_introgression_windows=min_introgression_windows,
    )
    return LandscapeResult(
        samples=samples,
        groups=groups,
        window_mode=mode,
        sample_rows=sample_rows,
        window_rows=window_rows,
        background_block_rows=block_rows,
        candidate_introgression_rows=introgression_rows,
        similarity_rows=similarity_rows,
    )


def build_background_blocks(
    sample_window_rows: list[dict[str, object]],
    min_background_similarity: float = 0.65,
) -> list[dict[str, object]]:
    """Merge adjacent sample windows with the same nearest local background."""
    grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in sample_window_rows:
        grouped[(str(row["sample"]), str(row["contig"]))].append(row)

    block_rows: list[dict[str, object]] = []
    block_index = 1
    for sample_contig in sorted(grouped):
        rows = sorted(grouped[sample_contig], key=lambda r: int(str(r["window_index"])))
        current: dict[str, object] | None = None
        similarity_values: list[float] = []

        for row in rows:
            nearest = str(row["nearest_background"])
            nearest_role = str(row["nearest_background_role"])
            similarity = _to_optional_float(row["nearest_similarity"])
            if nearest == "NA" or similarity is None or similarity < min_background_similarity:
                nearest = "unassigned"
                nearest_role = "unassigned"
                similarity = None

            should_extend = (
                current is not None
                and current["nearest_background"] == nearest
                and current["nearest_background_role"] == nearest_role
                and current["contig"] == row["contig"]
            )
            if not should_extend:
                if current is not None:
                    current["mean_similarity"] = _fmt_optional_mean(similarity_values)
                    block_rows.append(current)
                    block_index += 1
                current = {
                    "block_id": f"LB{block_index:08d}",
                    "sample": row["sample"],
                    "cohort_role": row["cohort_role"],
                    "contig": row["contig"],
                    "start": row["start"],
                    "end": row["end"],
                    "n_windows": 1,
                    "nearest_background": nearest,
                    "nearest_background_role": nearest_role,
                    "mean_similarity": "NA",
                }
                similarity_values = []
            else:
                assert current is not None
                current["end"] = row["end"]
                current["n_windows"] = int(str(current["n_windows"])) + 1

            if similarity is not None:
                similarity_values.append(similarity)

        if current is not None:
            current["mean_similarity"] = _fmt_optional_mean(similarity_values)
            block_rows.append(current)
            block_index += 1

    return block_rows


def build_candidate_introgression_blocks(
    sample_window_rows: list[dict[str, object]],
    similarity_rows: list[dict[str, object]],
    target_samples: tuple[str, ...],
    min_introgression_similarity: float = 0.65,
    min_introgression_delta: float = 0.0,
    max_introgression_missing_rate: float = 0.5,
    min_introgression_windows: int = 1,
) -> list[dict[str, object]]:
    """Merge adjacent donor-like target windows into candidate introgression blocks.

    These are exploratory local-background calls, not formal introgression tests.
    A target window is eligible when its nearest sample is an off-target sample,
    the nearest similarity passes the configured threshold, missingness is low
    enough, and the off-target similarity exceeds the best target-to-target
    similarity by at least ``min_introgression_delta`` when such a comparison is
    available.
    """
    target_set = set(target_samples)
    nearest_target_similarity = _nearest_target_similarity_by_window(
        similarity_rows,
        target_set,
    )

    grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in sample_window_rows:
        sample = str(row["sample"])
        if sample not in target_set or str(row["cohort_role"]) != "target":
            continue

        candidate_donor = str(row["nearest_background"])
        candidate_donor_role = str(row["nearest_background_role"])
        donor_similarity = _to_optional_float(row["nearest_similarity"])
        missing_rate = _to_optional_float(row["missing_rate"])
        if (
            candidate_donor == "NA"
            or candidate_donor_role != "off_target"
            or donor_similarity is None
            or donor_similarity < min_introgression_similarity
            or missing_rate is None
            or missing_rate > max_introgression_missing_rate
        ):
            continue

        target_similarity = nearest_target_similarity.get((str(row["window_id"]), sample))
        similarity_delta = (
            None if target_similarity is None else donor_similarity - target_similarity
        )
        if (
            similarity_delta is not None
            and similarity_delta < min_introgression_delta
        ):
            continue

        candidate = dict(row)
        candidate["candidate_donor"] = candidate_donor
        candidate["candidate_donor_role"] = candidate_donor_role
        candidate["donor_similarity"] = donor_similarity
        candidate["nearest_target_similarity"] = target_similarity
        candidate["similarity_delta"] = similarity_delta
        candidate["missing_rate_float"] = missing_rate
        grouped[(sample, str(row["contig"]))].append(candidate)

    block_rows: list[dict[str, object]] = []
    block_index = 1
    for sample_contig in sorted(grouped):
        rows = sorted(grouped[sample_contig], key=lambda r: int(str(r["window_index"])))
        current: dict[str, object] | None = None
        donor_values: list[float] = []
        target_values: list[float] = []
        delta_values: list[float] = []
        missing_values: list[float] = []
        private_values: list[float] = []
        nonref_values: list[float] = []

        for row in rows:
            window_index = int(str(row["window_index"]))
            should_extend = (
                current is not None
                and current["candidate_donor"] == row["candidate_donor"]
                and current["contig"] == row["contig"]
                and int(str(current["_last_window_index"])) + 1 == window_index
            )
            if not should_extend:
                if current is not None:
                    _finalize_candidate_introgression_block(
                        current=current,
                        donor_values=donor_values,
                        target_values=target_values,
                        delta_values=delta_values,
                        missing_values=missing_values,
                        private_values=private_values,
                        nonref_values=nonref_values,
                        min_introgression_windows=min_introgression_windows,
                        block_rows=block_rows,
                    )
                    if int(str(current["n_windows"])) >= min_introgression_windows:
                        block_index += 1
                current = {
                    "block_id": f"IB{block_index:08d}",
                    "sample": row["sample"],
                    "contig": row["contig"],
                    "start": row["start"],
                    "end": row["end"],
                    "n_windows": 1,
                    "candidate_donor": row["candidate_donor"],
                    "candidate_donor_role": row["candidate_donor_role"],
                    "mean_donor_similarity": "NA",
                    "mean_nearest_target_similarity": "NA",
                    "mean_similarity_delta": "NA",
                    "max_missing_rate": "NA",
                    "mean_private_alt_rate": "NA",
                    "mean_nonref_rate": "NA",
                    "evidence_class": "candidate_introgression",
                    "interpretation": (
                        "Target sample is locally closest to an off-target sample; "
                        "review as a candidate donor-like or introgressed block."
                    ),
                    "_last_window_index": window_index,
                }
                donor_values = []
                target_values = []
                delta_values = []
                missing_values = []
                private_values = []
                nonref_values = []
            else:
                assert current is not None
                current["end"] = row["end"]
                current["n_windows"] = int(str(current["n_windows"])) + 1
                current["_last_window_index"] = window_index

            donor_values.append(float(row["donor_similarity"]))
            if row["nearest_target_similarity"] is not None:
                target_values.append(float(row["nearest_target_similarity"]))
            if row["similarity_delta"] is not None:
                delta_values.append(float(row["similarity_delta"]))
            missing_values.append(float(row["missing_rate_float"]))
            private_rate = _to_optional_float(row["private_alt_rate"])
            nonref_rate = _to_optional_float(row["nonref_rate"])
            if private_rate is not None:
                private_values.append(private_rate)
            if nonref_rate is not None:
                nonref_values.append(nonref_rate)

        if current is not None:
            _finalize_candidate_introgression_block(
                current=current,
                donor_values=donor_values,
                target_values=target_values,
                delta_values=delta_values,
                missing_values=missing_values,
                private_values=private_values,
                nonref_values=nonref_values,
                min_introgression_windows=min_introgression_windows,
                block_rows=block_rows,
            )
            if int(str(current["n_windows"])) >= min_introgression_windows:
                block_index += 1

    return block_rows


def _nearest_target_similarity_by_window(
    similarity_rows: list[dict[str, object]],
    target_samples: set[str],
) -> dict[tuple[str, str], float]:
    nearest: dict[tuple[str, str], float] = {}
    for row in similarity_rows:
        sample_a = str(row["sample_a"])
        sample_b = str(row["sample_b"])
        if sample_a not in target_samples or sample_b not in target_samples:
            continue
        similarity = _to_optional_float(row["similarity"])
        if similarity is None:
            continue
        window_id = str(row["window_id"])
        for sample in (sample_a, sample_b):
            key = (window_id, sample)
            if key not in nearest or similarity > nearest[key]:
                nearest[key] = similarity
    return nearest


def _finalize_candidate_introgression_block(
    current: dict[str, object],
    donor_values: list[float],
    target_values: list[float],
    delta_values: list[float],
    missing_values: list[float],
    private_values: list[float],
    nonref_values: list[float],
    min_introgression_windows: int,
    block_rows: list[dict[str, object]],
) -> None:
    current["mean_donor_similarity"] = _fmt_optional_mean(donor_values)
    current["mean_nearest_target_similarity"] = _fmt_optional_mean(target_values)
    current["mean_similarity_delta"] = _fmt_optional_mean(delta_values)
    current["max_missing_rate"] = (
        "NA" if not missing_values else _fmt_float(max(missing_values))
    )
    current["mean_private_alt_rate"] = _fmt_optional_mean(private_values)
    current["mean_nonref_rate"] = _fmt_optional_mean(nonref_values)

    if not delta_values:
        current["evidence_class"] = "donor_like_single_target"
    elif mean(delta_values) > 0:
        current["evidence_class"] = "offtarget_closer_than_target"

    current.pop("_last_window_index", None)
    if int(str(current["n_windows"])) >= min_introgression_windows:
        block_rows.append(current)


def _analyze_contig_windows(
    records: tuple[_VariantSnapshot, ...],
    mode: WindowMode,
    window_counter_start: int,
    samples: tuple[str, ...],
    active_sample_indices: tuple[int, ...],
    sample_to_role: dict[str, str],
    groups: PangenomeGroups,
    sample_rows: list[dict[str, object]],
    window_rows: list[dict[str, object]],
    similarity_rows: list[dict[str, object]],
    background_rows_input: list[dict[str, object]],
    window_records: int,
    step_records: int,
    window_bp: int | None,
    step_bp: int | None,
    rare_max_count: int,
    rare_max_freq: float,
    min_called_for_freq: int,
    min_freq_values: int,
) -> int:
    if not records:
        return window_counter_start

    windows = _make_record_windows(records, window_records, step_records)
    if mode == "bp":
        if window_bp is None or step_bp is None:
            raise ValueError("bp window mode requires window_bp and step_bp.")
        windows = _make_bp_windows(records, window_bp, step_bp)

    window_counter = window_counter_start
    for local_index, window_variants in enumerate(windows, 1):
        window_counter += 1
        window = _make_window(
            variants=window_variants,
            mode=mode,
            global_index=window_counter,
            local_index=local_index,
        )
        rows = _analyze_window(
            window=window,
            samples=samples,
            active_sample_indices=active_sample_indices,
            sample_to_role=sample_to_role,
            groups=groups,
            rare_max_count=rare_max_count,
            rare_max_freq=rare_max_freq,
            min_called_for_freq=min_called_for_freq,
            min_freq_values=min_freq_values,
        )
        sample_rows.extend(rows.sample_rows)
        window_rows.append(rows.window_row)
        similarity_rows.extend(rows.similarity_rows)
        background_rows_input.extend(rows.sample_rows)

    return window_counter


@dataclass(frozen=True)
class _AnalyzedWindow:
    sample_rows: list[dict[str, object]]
    window_row: dict[str, object]
    similarity_rows: list[dict[str, object]]


def _analyze_window(
    window: _Window,
    samples: tuple[str, ...],
    active_sample_indices: tuple[int, ...],
    sample_to_role: dict[str, str],
    groups: PangenomeGroups,
    rare_max_count: int,
    rare_max_freq: float,
    min_called_for_freq: int,
    min_freq_values: int,
) -> _AnalyzedWindow:
    n_samples = len(samples)
    n_variants = len(window.variants)
    called_n = [0] * n_samples
    missing_n = [0] * n_samples
    het_n = [0] * n_samples
    nonref_n = [0] * n_samples
    minor_gt_n = [0] * n_samples
    rare_alt_n = [0] * n_samples
    private_alt_n = [0] * n_samples
    call_freqs: list[list[float]] = [[] for _ in samples]
    private_alt_events_target = 0
    private_alt_events_offtarget = 0

    target_indices = {samples.index(sample) for sample in groups.target}
    offtarget_indices = {samples.index(sample) for sample in groups.off_target}

    pair_match_counts: dict[tuple[int, int], int] = defaultdict(int)
    pair_compared_counts: dict[tuple[int, int], int] = defaultdict(int)

    for variant in window.variants:
        genotype_counts: Counter[Genotype] = Counter()
        called_indices: list[int] = []
        for sample_index in active_sample_indices:
            gt = variant.genotypes[sample_index]
            if is_missing_genotype(gt):
                missing_n[sample_index] += 1
                continue
            if gt is None:
                missing_n[sample_index] += 1
                continue
            normalized_gt = _normalize_genotype(gt)
            genotype_counts[normalized_gt] += 1
            called_n[sample_index] += 1
            called_indices.append(sample_index)
            if _is_het(gt):
                het_n[sample_index] += 1
            if _is_nonref(gt):
                nonref_n[sample_index] += 1

        if len(genotype_counts) >= 2:
            min_count = min(genotype_counts.values())
            minor_gts = {gt for gt, count in genotype_counts.items() if count == min_count}
            for sample_index in called_indices:
                gt = variant.genotypes[sample_index]
                if gt is not None and _normalize_genotype(gt) in minor_gts:
                    minor_gt_n[sample_index] += 1

        called_denominator = len(called_indices)
        if called_denominator > min_called_for_freq:
            for sample_index in called_indices:
                gt = variant.genotypes[sample_index]
                if gt is not None:
                    call_freqs[sample_index].append(
                        genotype_counts[_normalize_genotype(gt)] / called_denominator
                    )

        for carriers in variant.alt_carriers:
            active_carriers = set(carriers) & set(active_sample_indices)
            if not active_carriers:
                continue
            alt_freq = len(active_carriers) / len(active_sample_indices)
            is_rare = len(active_carriers) <= rare_max_count or alt_freq <= rare_max_freq
            target_carriers = active_carriers & target_indices
            offtarget_carriers = active_carriers & offtarget_indices
            is_target_private = bool(target_carriers) and not offtarget_carriers
            is_offtarget_private = bool(offtarget_carriers) and not target_carriers
            if is_target_private:
                private_alt_events_target += 1
            if is_offtarget_private:
                private_alt_events_offtarget += 1
            for sample_index in active_carriers:
                if is_rare:
                    rare_alt_n[sample_index] += 1
                if (
                    (sample_index in target_indices and is_target_private)
                    or (sample_index in offtarget_indices and is_offtarget_private)
                ):
                    private_alt_n[sample_index] += 1

        for left_pos, left_index in enumerate(called_indices):
            left_gt = variant.genotypes[left_index]
            if left_gt is None:
                continue
            for right_index in called_indices[left_pos + 1:]:
                right_gt = variant.genotypes[right_index]
                if right_gt is None:
                    continue
                pair = (left_index, right_index)
                pair_compared_counts[pair] += 1
                if _normalize_genotype(left_gt) == _normalize_genotype(right_gt):
                    pair_match_counts[pair] += 1

    pairwise = _pairwise_similarity(samples, active_sample_indices, pair_match_counts,
                                    pair_compared_counts)
    nearest = _nearest_backgrounds(samples, active_sample_indices, pairwise)

    midpoint = window.start + ((window.end - window.start) // 2)
    sample_rows: list[dict[str, object]] = []
    for sample_index in active_sample_indices:
        sample = samples[sample_index]
        nearest_index, nearest_stats = nearest[sample_index]
        nearest_sample = "NA" if nearest_index is None else samples[nearest_index]
        nearest_role = (
            "NA" if nearest_index is None else sample_to_role.get(nearest_sample, "other")
        )
        median_freq: str | float = "NA"
        if len(call_freqs[sample_index]) > min_freq_values:
            median_freq = _fmt_float(median(call_freqs[sample_index]))
        row = {
            "window_id": window.window_id,
            "contig": window.contig,
            "window_index": window.window_index,
            "start": window.start,
            "end": window.end,
            "midpoint": midpoint,
            "window_mode": window.mode,
            "n_variants": n_variants,
            "sample": sample,
            "cohort_role": sample_to_role.get(sample, "other"),
            "called_n": called_n[sample_index],
            "missing_n": missing_n[sample_index],
            "missing_rate": _fmt_ratio(missing_n[sample_index], n_variants),
            "het_n": het_n[sample_index],
            "het_rate": _fmt_ratio(het_n[sample_index], called_n[sample_index]),
            "nonref_n": nonref_n[sample_index],
            "nonref_rate": _fmt_ratio(nonref_n[sample_index], called_n[sample_index]),
            "minor_genotype_n": minor_gt_n[sample_index],
            "minor_genotype_rate": _fmt_ratio(minor_gt_n[sample_index], called_n[sample_index]),
            "rare_alt_n": rare_alt_n[sample_index],
            "rare_alt_rate": _fmt_ratio(rare_alt_n[sample_index], called_n[sample_index]),
            "private_alt_n": private_alt_n[sample_index],
            "private_alt_rate": _fmt_ratio(private_alt_n[sample_index], called_n[sample_index]),
            "median_call_freq": median_freq,
            "nearest_background": nearest_sample,
            "nearest_background_role": nearest_role,
            "nearest_similarity": (
                "NA" if nearest_stats.similarity is None else _fmt_float(nearest_stats.similarity)
            ),
            "similarity_compared_variants": nearest_stats.compared,
        }
        sample_rows.append(row)

    similarity_rows = _similarity_rows(window, samples, active_sample_indices, pairwise)
    target_rows = [r for r in sample_rows if r["cohort_role"] == "target"]
    offtarget_rows = [r for r in sample_rows if r["cohort_role"] == "off_target"]
    span_bp = max(window.end - window.start, 1)
    nearest_counter = Counter(str(r["nearest_background"]) for r in sample_rows
                              if r["nearest_background"] != "NA")
    top_nearest = nearest_counter.most_common(1)

    window_row = {
        "window_id": window.window_id,
        "contig": window.contig,
        "window_index": window.window_index,
        "start": window.start,
        "end": window.end,
        "midpoint": midpoint,
        "window_mode": window.mode,
        "n_variants": n_variants,
        "span_bp": span_bp,
        "density_variants_per_kb": _fmt_float(n_variants / (span_bp / 1000)),
        "target_mean_missing_rate": _mean_row_float(target_rows, "missing_rate"),
        "offtarget_mean_missing_rate": _mean_row_float(offtarget_rows, "missing_rate"),
        "target_mean_nonref_rate": _mean_row_float(target_rows, "nonref_rate"),
        "offtarget_mean_nonref_rate": _mean_row_float(offtarget_rows, "nonref_rate"),
        "target_private_alt_n": private_alt_events_target,
        "offtarget_private_alt_n": private_alt_events_offtarget,
        "target_private_alt_rate": _fmt_ratio(private_alt_events_target, n_variants),
        "offtarget_private_alt_rate": _fmt_ratio(private_alt_events_offtarget, n_variants),
        "top_nearest_background": top_nearest[0][0] if top_nearest else "NA",
        "top_nearest_background_n": top_nearest[0][1] if top_nearest else 0,
    }
    return _AnalyzedWindow(
        sample_rows=sample_rows,
        window_row=window_row,
        similarity_rows=similarity_rows,
    )


def _snapshot_record(record: VariantRecordLike, samples: tuple[str, ...]) -> _VariantSnapshot:
    alts = record.alts
    genotypes: list[Genotype | None] = []
    for sample in samples:
        gt: Genotype | None
        try:
            gt = record.samples[sample]["GT"]
        except (KeyError, TypeError):
            gt = None
        genotypes.append(gt)

    alt_carriers: list[frozenset[int]] = []
    for alt_index, _alt in enumerate(alts or ()):
        carriers = {
            sample_index
            for sample_index, gt in enumerate(genotypes)
            if gt is not None and not is_missing_genotype(gt) and has_alt_allele(gt, alt_index)
        }
        alt_carriers.append(frozenset(carriers))

    start = record.pos - 1
    end = start + max(len(record.ref), 1)
    return _VariantSnapshot(
        contig=record.chrom,
        start=start,
        end=end,
        genotypes=tuple(genotypes),
        alt_carriers=tuple(alt_carriers),
    )


def _make_record_windows(
    records: tuple[_VariantSnapshot, ...],
    window_records: int,
    step_records: int,
) -> list[tuple[_VariantSnapshot, ...]]:
    windows: list[tuple[_VariantSnapshot, ...]] = []
    start = 0
    n_records = len(records)
    while start < n_records:
        end = min(start + window_records, n_records)
        windows.append(records[start:end])
        if end == n_records:
            break
        start += step_records
    return windows


def _make_bp_windows(
    records: tuple[_VariantSnapshot, ...],
    window_bp: int,
    step_bp: int,
) -> list[tuple[_VariantSnapshot, ...]]:
    if not records:
        return []
    contig_start = min(record.start for record in records)
    contig_end = max(record.end for record in records)
    windows: list[tuple[_VariantSnapshot, ...]] = []
    start = contig_start
    while start < contig_end:
        end = start + window_bp
        variants = tuple(record for record in records if record.start < end and record.end > start)
        if variants:
            windows.append(variants)
        if end >= contig_end:
            break
        start += step_bp
    return windows


def _make_window(
    variants: tuple[_VariantSnapshot, ...],
    mode: WindowMode,
    global_index: int,
    local_index: int,
) -> _Window:
    if not variants:
        raise ValueError("Cannot create a window without variants.")
    contig = variants[0].contig
    start = min(variant.start for variant in variants)
    end = max(variant.end for variant in variants)
    return _Window(
        window_id=f"LW{global_index:08d}",
        contig=contig,
        window_index=local_index,
        start=start,
        end=end,
        mode=mode,
        variants=variants,
    )


def _pairwise_similarity(
    samples: tuple[str, ...],
    active_sample_indices: tuple[int, ...],
    pair_match_counts: dict[tuple[int, int], int],
    pair_compared_counts: dict[tuple[int, int], int],
) -> dict[tuple[int, int], _PairwiseStats]:
    del samples
    pairwise: dict[tuple[int, int], _PairwiseStats] = {}
    for left_pos, left_index in enumerate(active_sample_indices):
        for right_index in active_sample_indices[left_pos + 1:]:
            pair = (left_index, right_index)
            compared = pair_compared_counts.get(pair, 0)
            similarity = None if compared == 0 else pair_match_counts.get(pair, 0) / compared
            pairwise[pair] = _PairwiseStats(similarity=similarity, compared=compared)
    return pairwise


def _nearest_backgrounds(
    samples: tuple[str, ...],
    active_sample_indices: tuple[int, ...],
    pairwise: dict[tuple[int, int], _PairwiseStats],
) -> dict[int, tuple[int | None, _PairwiseStats]]:
    nearest: dict[int, tuple[int | None, _PairwiseStats]] = {}
    for sample_index in active_sample_indices:
        best_index: int | None = None
        best_stats = _PairwiseStats(similarity=None, compared=0)
        for other_index in active_sample_indices:
            if sample_index == other_index:
                continue
            pair = (
                min(sample_index, other_index),
                max(sample_index, other_index),
            )
            stats = pairwise.get(pair, _PairwiseStats(similarity=None, compared=0))
            if stats.similarity is None:
                continue
            if best_stats.similarity is None or stats.similarity > best_stats.similarity:
                best_index = other_index
                best_stats = stats
            elif stats.similarity == best_stats.similarity and best_index is not None:
                if samples[other_index] < samples[best_index]:
                    best_index = other_index
                    best_stats = stats
        nearest[sample_index] = (best_index, best_stats)
    return nearest


def _similarity_rows(
    window: _Window,
    samples: tuple[str, ...],
    active_sample_indices: tuple[int, ...],
    pairwise: dict[tuple[int, int], _PairwiseStats],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for left_pos, left_index in enumerate(active_sample_indices):
        for right_index in active_sample_indices[left_pos + 1:]:
            pair = (left_index, right_index)
            stats = pairwise[pair]
            rows.append({
                "window_id": window.window_id,
                "contig": window.contig,
                "window_index": window.window_index,
                "start": window.start,
                "end": window.end,
                "sample_a": samples[left_index],
                "sample_b": samples[right_index],
                "similarity": "NA" if stats.similarity is None else _fmt_float(stats.similarity),
                "compared_variants": stats.compared,
            })
    return rows


def _sample_roles(groups: PangenomeGroups) -> dict[str, str]:
    roles = {sample: "target" for sample in groups.target}
    roles.update({sample: "off_target" for sample in groups.off_target})
    return roles


def _record_passes(record_filter_values: list[str]) -> bool:
    return (
        not record_filter_values
        or record_filter_values == ["PASS"]
        or "PASS" in record_filter_values
    )


def _normalize_genotype(gt: Genotype) -> Genotype:
    return tuple(sorted(gt, key=lambda allele: -1 if allele is None else allele))


def _is_het(gt: Genotype) -> bool:
    called = [allele for allele in gt if allele is not None]
    return len(set(called)) > 1


def _is_nonref(gt: Genotype) -> bool:
    return any(allele is not None and allele > 0 for allele in gt)


def _fmt_ratio(numerator: int, denominator: int) -> str:
    if denominator == 0:
        return "NA"
    return _fmt_float(numerator / denominator)


def _fmt_float(value: float) -> str:
    return f"{value:.6f}"


def _mean_row_float(rows: list[dict[str, object]], key: str) -> str:
    values = [_to_optional_float(row[key]) for row in rows]
    usable = [value for value in values if value is not None]
    return _fmt_optional_mean(usable)


def _fmt_optional_mean(values: list[float]) -> str:
    if not values:
        return "NA"
    return _fmt_float(mean(values))


def _to_optional_float(value: object) -> float | None:
    if value == "NA" or value is None:
        return None
    try:
        return float(str(value))
    except (TypeError, ValueError):
        return None
