"""Unit tests for VCF landscape primitives."""

from __future__ import annotations

from pathlib import Path

from privy.landscape.vcf import run_vcf_landscape


def test_vcf_landscape_record_windows_emit_sample_and_window_rows(indexed_vcf: Path) -> None:
    result = run_vcf_landscape(
        indexed_vcf,
        targets=["T1", "T2"],
        window_records=3,
        step_records=3,
        min_called_for_freq=0,
        min_freq_values=0,
    )

    assert result.window_mode == "records"
    assert len(result.window_rows) == 3
    assert len(result.sample_rows) == 15
    first_window = result.window_rows[0]
    assert first_window["contig"] == "chr1"
    assert first_window["n_variants"] == 3
    assert first_window["target_private_alt_n"] == 3

    t1_first = next(
        row for row in result.sample_rows
        if row["window_id"] == first_window["window_id"] and row["sample"] == "T1"
    )
    assert t1_first["missing_rate"] == "0.000000"
    assert t1_first["private_alt_n"] == 3
    assert t1_first["nearest_background"] == "T2"
    assert t1_first["nearest_similarity"] == "1.000000"


def test_vcf_landscape_bp_windows_are_user_adjustable(indexed_vcf: Path) -> None:
    result = run_vcf_landscape(
        indexed_vcf,
        targets=["T1", "T2"],
        window_bp=250,
        step_bp=250,
        min_called_for_freq=0,
        min_freq_values=0,
    )

    assert result.window_mode == "bp"
    assert result.window_rows
    assert all(row["window_mode"] == "bp" for row in result.window_rows)
    assert result.window_rows[0]["span_bp"] <= 250


def test_vcf_landscape_background_blocks_merge_nearest_assignments(indexed_vcf: Path) -> None:
    result = run_vcf_landscape(
        indexed_vcf,
        targets=["T1", "T2"],
        window_records=3,
        step_records=3,
        min_called_for_freq=0,
        min_freq_values=0,
        min_background_similarity=0.5,
    )

    t1_blocks = [row for row in result.background_block_rows if row["sample"] == "T1"]
    assert t1_blocks
    assert any(row["nearest_background"] == "T2" for row in t1_blocks)
