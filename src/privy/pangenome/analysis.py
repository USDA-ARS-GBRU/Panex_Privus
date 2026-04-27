"""Shared pangenome analyses for feature matrices."""

from __future__ import annotations

import random
from collections import Counter
from collections.abc import Iterable

from privy.pangenome.model import FeatureMatrix, FeatureRecord, PangenomeGroups


def build_feature_summary_rows(
    matrix: FeatureMatrix,
    groups: PangenomeGroups,
) -> list[dict[str, object]]:
    """Summarize each feature across full, target, and off-target groups."""
    rows: list[dict[str, object]] = []
    target_set = set(groups.target)
    off_target_set = set(groups.off_target)

    for feature in matrix.features:
        samples = matrix.samples_for_feature(feature.feature_id)
        target_count = len(samples & target_set)
        off_target_count = len(samples & off_target_set)
        total_count = target_count + off_target_count
        rows.append({
            "feature_id": feature.feature_id,
            "source_type": feature.source_type,
            "feature_type": feature.feature_type,
            "contig": feature.contig or "",
            "start": "" if feature.start is None else feature.start,
            "end": "" if feature.end is None else feature.end,
            "length": feature.length,
            "total_present_n": total_count,
            "target_present_n": target_count,
            "target_total_n": len(groups.target),
            "offtarget_present_n": off_target_count,
            "offtarget_total_n": len(groups.off_target),
            "full_category": _category(total_count, len(groups.full)),
            "target_category": _category(target_count, len(groups.target)),
            "offtarget_category": _category(off_target_count, len(groups.off_target)),
            "target_private": target_count > 0 and off_target_count == 0,
            "offtarget_private": off_target_count > 0 and target_count == 0,
        })

    return rows


def build_coverage_histogram_rows(
    matrix: FeatureMatrix,
    groups: PangenomeGroups,
) -> list[dict[str, object]]:
    """Build coverage histograms by group."""
    rows: list[dict[str, object]] = []
    for group_name, samples in groups.as_dict().items():
        sample_set = set(samples)
        feature_hist: Counter[int] = Counter()
        bp_hist: Counter[int] = Counter()
        for feature in matrix.features:
            coverage = len(matrix.samples_for_feature(feature.feature_id) & sample_set)
            feature_hist[coverage] += 1
            bp_hist[coverage] += feature.length

        for coverage in range(len(samples) + 1):
            rows.append({
                "group": group_name,
                "coverage": coverage,
                "n_features": feature_hist[coverage],
                "n_bp": bp_hist[coverage],
            })
    return rows


def build_composition_rows(
    matrix: FeatureMatrix,
    groups: PangenomeGroups,
) -> list[dict[str, object]]:
    """Count absent/private/accessory/core features by group."""
    rows: list[dict[str, object]] = []
    for group_name, samples in groups.as_dict().items():
        sample_set = set(samples)
        category_features: Counter[str] = Counter()
        category_bp: Counter[str] = Counter()
        for feature in matrix.features:
            coverage = len(matrix.samples_for_feature(feature.feature_id) & sample_set)
            category = _category(coverage, len(samples))
            category_features[category] += 1
            category_bp[category] += feature.length

        for category in ("absent", "private", "accessory", "core"):
            rows.append({
                "group": group_name,
                "category": category,
                "n_features": category_features[category],
                "n_bp": category_bp[category],
            })
    return rows


def build_growth_curve_rows(
    matrix: FeatureMatrix,
    groups: PangenomeGroups,
    permutations: int = 100,
    seed: int = 42,
) -> list[dict[str, object]]:
    """Build permutation-based pangenome growth curves for each built-in group."""
    if permutations < 1:
        raise ValueError("permutations must be at least 1.")

    rows: list[dict[str, object]] = []
    for group_name, samples in groups.as_dict().items():
        rows.extend(_growth_rows_for_group(matrix.features, matrix.presence, group_name, samples,
                                           permutations, seed))
    return rows


def _growth_rows_for_group(
    features: Iterable[FeatureRecord],
    presence: dict[str, frozenset[str]],
    group_name: str,
    samples: tuple[str, ...],
    permutations: int,
    seed: int,
) -> list[dict[str, object]]:
    feature_list = tuple(features)
    rows: list[dict[str, object]] = []
    if not samples:
        return rows

    for trial in range(1, permutations + 1):
        order = list(samples)
        random.Random(seed + trial - 1).shuffle(order)
        selected: set[str] = set()
        prev_features = 0
        prev_bp = 0

        for n, sample in enumerate(order, 1):
            selected.add(sample)
            present_features = 0
            present_bp = 0
            singleton_features = 0
            singleton_bp = 0

            for feature in feature_list:
                count = len(presence.get(feature.feature_id, frozenset()) & selected)
                if count > 0:
                    present_features += 1
                    present_bp += feature.length
                if count == 1:
                    singleton_features += 1
                    singleton_bp += feature.length

            rows.append({
                "group": group_name,
                "trial": trial,
                "n": n,
                "sample_added": sample,
                "features": present_features,
                "bp": present_bp,
                "new_features": present_features - prev_features,
                "new_bp": present_bp - prev_bp,
                "singleton_features": singleton_features,
                "singleton_bp": singleton_bp,
            })
            prev_features = present_features
            prev_bp = present_bp

    return rows


def _category(present_n: int, group_n: int) -> str:
    if present_n == 0:
        return "absent"
    if present_n == group_n:
        return "core"
    if present_n == 1:
        return "private"
    return "accessory"
