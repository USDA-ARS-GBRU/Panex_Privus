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


def test_vcf_landscape_reports_candidate_introgression_blocks(tmp_path: Path) -> None:
    vcf = _write_indexed_vcf(
        tmp_path,
        "\n".join([
            "##fileformat=VCFv4.2",
            '##FILTER=<ID=PASS,Description="All filters passed">',
            "##contig=<ID=chr1,length=10000>",
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tT1\tT2\tO1\tO2",
            "chr1\t100\t.\tA\tT\t50\tPASS\t.\tGT\t0/1\t1/1\t0/1\t0/0",
            "chr1\t200\t.\tA\tT\t50\tPASS\t.\tGT\t0/1\t1/1\t0/1\t0/0",
            "chr1\t300\t.\tA\tT\t50\tPASS\t.\tGT\t0/1\t1/1\t0/1\t0/0",
            "chr1\t400\t.\tA\tT\t50\tPASS\t.\tGT\t0/1\t1/1\t0/1\t0/0",
            "chr1\t500\t.\tA\tT\t50\tPASS\t.\tGT\t0/1\t1/1\t0/1\t0/0",
            "chr1\t600\t.\tA\tT\t50\tPASS\t.\tGT\t0/1\t1/1\t0/1\t0/0",
            "",
        ]),
    )

    result = run_vcf_landscape(
        vcf,
        targets=["T1", "T2"],
        off_targets=["O1", "O2"],
        window_records=3,
        step_records=3,
        min_called_for_freq=0,
        min_freq_values=0,
        min_introgression_similarity=0.9,
        min_introgression_delta=0.5,
        min_introgression_windows=2,
    )

    assert len(result.candidate_introgression_rows) == 1
    block = result.candidate_introgression_rows[0]
    assert block["sample"] == "T1"
    assert block["candidate_donor"] == "O1"
    assert block["n_windows"] == 2
    assert block["mean_donor_similarity"] == "1.000000"
    assert block["mean_similarity_delta"] == "1.000000"


def _write_indexed_vcf(tmp_path: Path, text: str) -> Path:
    import pysam  # noqa: PLC0415

    plain = tmp_path / "introgression.vcf"
    plain.write_text(text, encoding="utf-8")
    gz_path = str(tmp_path / "introgression.vcf.gz")
    pysam.tabix_compress(str(plain), gz_path, force=True)
    pysam.tabix_index(gz_path, preset="vcf", force=True)
    return Path(gz_path)
