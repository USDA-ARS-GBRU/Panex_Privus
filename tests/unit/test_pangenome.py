"""Unit tests for shared pangenome analysis primitives."""

from __future__ import annotations

from pathlib import Path

import pytest

from privy.io.gfa import parse_gfa
from privy.pangenome import (
    build_composition_rows,
    build_coverage_histogram_rows,
    build_feature_summary_rows,
    build_gfa_feature_matrix,
    build_growth_curve_rows,
    build_vcf_feature_matrix,
    resolve_pangenome_groups,
)

GFA_PATH = Path(__file__).parent.parent / "data" / "small_cohort.gfa"


def test_resolve_groups_infers_off_targets_from_remaining_samples() -> None:
    groups = resolve_pangenome_groups(
        all_samples=["T1", "T2", "O1", "O2", "O3"],
        targets=["T1", "T2"],
    )

    assert groups.target == ("T1", "T2")
    assert groups.off_target == ("O1", "O2", "O3")
    assert groups.full == ("T1", "T2", "O1", "O2", "O3")


def test_resolve_groups_rejects_unknown_target() -> None:
    with pytest.raises(ValueError, match="Target samples were not found"):
        resolve_pangenome_groups(
            all_samples=["T1", "O1"],
            targets=["T1", "T2"],
        )


def test_build_gfa_feature_matrix_uses_segments_as_features() -> None:
    graph = parse_gfa(GFA_PATH)
    matrix = build_gfa_feature_matrix(graph)

    assert len(matrix.features) == 7
    assert matrix.source_type == "gfa"
    assert set(matrix.samples) == {"T1", "T2", "O1", "O2", "O3"}
    assert matrix.samples_for_feature("s2_target") == frozenset({"T1", "T2"})
    assert matrix.samples_for_feature("s2_offt") == frozenset({"O1", "O2", "O3"})


def test_build_vcf_feature_matrix_uses_alt_alleles_as_features(indexed_vcf: Path) -> None:
    matrix = build_vcf_feature_matrix(indexed_vcf)

    assert matrix.source_type == "vcf"
    assert len(matrix.features) == 10
    assert set(matrix.samples) == {"T1", "T2", "O1", "O2", "O3"}
    assert matrix.samples_for_feature("chr1:100:A:T") == frozenset({"T1", "T2"})
    assert matrix.samples_for_feature("chr1:500:A:T") == frozenset({"T1", "T2", "O1"})
    assert matrix.samples_for_feature("chr1:800:A:G") == frozenset()


def test_vcf_feature_summary_marks_target_private_alleles(indexed_vcf: Path) -> None:
    matrix = build_vcf_feature_matrix(indexed_vcf)
    groups = resolve_pangenome_groups(matrix.samples, targets=["T1", "T2"])
    rows = build_feature_summary_rows(matrix, groups)

    by_id = {str(row["feature_id"]): row for row in rows}
    assert by_id["chr1:100:A:T"]["target_private"] is True
    assert by_id["chr1:500:A:T"]["target_private"] is False
    assert by_id["chr1:800:A:G"]["full_category"] == "absent"


def test_feature_summary_marks_target_private_segments() -> None:
    matrix = build_gfa_feature_matrix(parse_gfa(GFA_PATH))
    groups = resolve_pangenome_groups(matrix.samples, targets=["T1", "T2"])
    rows = build_feature_summary_rows(matrix, groups)

    by_id = {str(row["feature_id"]): row for row in rows}
    assert by_id["s2_target"]["target_private"] is True
    assert by_id["s4_target"]["target_private"] is True
    assert by_id["s1"]["full_category"] == "core"
    assert by_id["s2_offt"]["offtarget_category"] == "core"


def test_coverage_histogram_counts_full_graph_coverage() -> None:
    matrix = build_gfa_feature_matrix(parse_gfa(GFA_PATH))
    groups = resolve_pangenome_groups(matrix.samples, targets=["T1", "T2"])
    rows = build_coverage_histogram_rows(matrix, groups)

    full = {int(row["coverage"]): row for row in rows if row["group"] == "full"}
    assert full[1]["n_features"] == 1  # s4_target
    assert full[2]["n_features"] == 1  # s2_target
    assert full[3]["n_features"] == 2  # s2_offt, s4_offt
    assert full[4]["n_features"] == 1  # s5, because T2 is missing at bubble 2
    assert full[5]["n_features"] == 2  # s1, s3 backbone in all samples


def test_composition_rows_are_emitted_for_each_group() -> None:
    matrix = build_gfa_feature_matrix(parse_gfa(GFA_PATH))
    groups = resolve_pangenome_groups(matrix.samples, targets=["T1", "T2"])
    rows = build_composition_rows(matrix, groups)

    assert len(rows) == 12
    target = {
        str(row["category"]): row
        for row in rows
        if row["group"] == "target"
    }
    assert target["core"]["n_features"] == 3
    assert target["private"]["n_features"] == 2
    assert target["absent"]["n_features"] == 2


def test_growth_curve_rows_are_deterministic() -> None:
    matrix = build_gfa_feature_matrix(parse_gfa(GFA_PATH))
    groups = resolve_pangenome_groups(matrix.samples, targets=["T1", "T2"])

    first = build_growth_curve_rows(matrix, groups, permutations=2, seed=7)
    second = build_growth_curve_rows(matrix, groups, permutations=2, seed=7)

    assert first == second
    assert len(first) == 20  # 2 permutations * (5 full + 2 target + 3 off-target)
