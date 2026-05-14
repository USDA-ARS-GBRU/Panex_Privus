"""Tests for shared CLI cohort parsing helpers."""

from __future__ import annotations

from pathlib import Path

from privy.cli.cohort_args import (
    collect_sample_values,
    load_cohort_file,
    parse_grouped_cohort_args,
)


def test_parse_grouped_cohort_args_accepts_repeated_and_grouped_values() -> None:
    parsed = parse_grouped_cohort_args(
        [
            "--targets",
            "T1",
            "T2",
            "--off-targets",
            "O1",
            "--off-targets",
            "O2",
            "--ignore-samples=LowQual",
        ]
    )

    assert parsed["targets"] == ["T1", "T2"]
    assert parsed["off_targets"] == ["O1", "O2"]
    assert parsed["ignore_samples"] == ["LowQual"]


def test_collect_sample_values_merges_flags_and_file(tmp_path: Path) -> None:
    samples_file = tmp_path / "targets.txt"
    samples_file.write_text("# comment\nT2\nT1\nT3 extra-column\n", encoding="utf-8")

    samples = collect_sample_values(["T1"], samples_file)

    assert samples == ["T1", "T2", "T3"]


def test_load_cohort_file_reads_yaml_with_nested_cohorts(tmp_path: Path) -> None:
    cohort_yaml = tmp_path / "cohort.yaml"
    cohort_yaml.write_text(
        "cohorts:\n"
        "  targets: [T1, T2]\n"
        "  off_targets: [O1, O2]\n"
        "  ignored_samples: [LowQual]\n",
        encoding="utf-8",
    )

    spec = load_cohort_file(cohort_yaml)

    assert spec.targets == ("T1", "T2")
    assert spec.off_targets == ("O1", "O2")
    assert spec.ignored_samples == ("LowQual",)


def test_load_cohort_file_reads_tsv_roles(tmp_path: Path) -> None:
    cohort_tsv = tmp_path / "cohort.tsv"
    cohort_tsv.write_text(
        "sample_id\tcohort_role\n"
        "T1\ttarget\n"
        "T2\ttarget\n"
        "O1\toff-target\n"
        "O2\tbackground\n"
        "LowQual\tignored\n",
        encoding="utf-8",
    )

    spec = load_cohort_file(cohort_tsv)

    assert spec.targets == ("T1", "T2")
    assert spec.off_targets == ("O1", "O2")
    assert spec.ignored_samples == ("LowQual",)
