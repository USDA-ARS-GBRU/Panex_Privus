"""End-to-end integration tests for the VCF scan backend.

Runs :func:`~privy.backends.vcf_scan.run_vcf_scan` on the small synthetic
VCF defined in ``conftest.py`` and asserts on the written output files.

Expected behaviour with default config (pass_only=True, min_qual=None,
allow_multiallelic=True):

    pos=100  strict_complete   (QUAL 50)   → emitted  (PPX00000001)
    pos=200  strict_target_missing          → emitted  (PPX00000002)
    pos=300  strict_offtarget_missing       → emitted  (PPX00000003)
    pos=400  strict_both_missing            → emitted  (PPX00000004)
    pos=500  contradicted                   → NOT emitted
    pos=600  FILTER=FAIL                    → skipped (records_skipped_filter)
    pos=700  strict_complete   (QUAL  5)   → emitted  (PPX00000005)
    pos=800  ALT=T  strict_complete        → emitted  (PPX00000006)
    pos=800  ALT=G  RELAXED_THRESHOLD      → NOT emitted (target support 0/2)
    pos=900  strict_complete (indel)        → emitted  (PPX00000007)

Total: 7 hits, 7 regions (merge_distance=0 default).
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from privy.backends.vcf_scan import run_vcf_scan


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _qc(outdir: Path) -> dict[str, str]:
    """Return qc.tsv as a {metric: value} dict."""
    rows = _read_tsv(outdir / "qc.tsv")
    return {r["metric"]: r["value"] for r in rows}


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def outdir(tmp_path: Path) -> Path:
    d = tmp_path / "out"
    d.mkdir()
    return d


# ---------------------------------------------------------------------------
# Output files
# ---------------------------------------------------------------------------


class TestOutputFiles:
    def test_all_six_output_files_written(
        self, indexed_vcf, small_cohort, default_cfg, outdir
    ) -> None:
        run_vcf_scan(vcf=indexed_vcf, cohort=small_cohort, cfg=default_cfg, outdir=outdir)
        for fname in (
            "hits.tsv",
            "regions.tsv",
            "evidence.tsv",
            "sample_support.tsv",
            "qc.tsv",
            "run.json",
        ):
            assert (outdir / fname).exists(), f"{fname} was not written"

    def test_run_json_is_parseable(
        self, indexed_vcf, small_cohort, default_cfg, outdir
    ) -> None:
        run_vcf_scan(vcf=indexed_vcf, cohort=small_cohort, cfg=default_cfg, outdir=outdir)
        data = json.loads((outdir / "run.json").read_text())
        assert "privy_version" in data
        assert "scan_stats" in data
        assert "config" in data

    def test_run_json_cohort_sizes(
        self, indexed_vcf, small_cohort, default_cfg, outdir
    ) -> None:
        run_vcf_scan(vcf=indexed_vcf, cohort=small_cohort, cfg=default_cfg, outdir=outdir)
        data = json.loads((outdir / "run.json").read_text())
        assert data["scan_stats"]["n_target_samples"] == 2
        assert data["scan_stats"]["n_offtarget_samples"] == 3


# ---------------------------------------------------------------------------
# hits.tsv content
# ---------------------------------------------------------------------------


class TestHitsTsv:
    def test_correct_column_names(
        self, indexed_vcf, small_cohort, default_cfg, outdir
    ) -> None:
        run_vcf_scan(vcf=indexed_vcf, cohort=small_cohort, cfg=default_cfg, outdir=outdir)
        rows = _read_tsv(outdir / "hits.tsv")
        assert rows, "hits.tsv should not be empty"
        expected = {
            "locus_id", "contig", "start", "end", "variant_type", "allele_key",
            "target_support_n", "target_total_n", "offtarget_support_n",
            "offtarget_total_n", "target_missing_n", "offtarget_missing_n",
            "strictness_class", "discovery_score", "support_score",
            "penalty_score", "final_score",
        }
        assert set(rows[0].keys()) == expected

    def test_default_hit_count_is_seven(
        self, indexed_vcf, small_cohort, default_cfg, outdir
    ) -> None:
        run_vcf_scan(vcf=indexed_vcf, cohort=small_cohort, cfg=default_cfg, outdir=outdir)
        rows = _read_tsv(outdir / "hits.tsv")
        assert len(rows) == 7

    def test_strictness_class_distribution(
        self, indexed_vcf, small_cohort, default_cfg, outdir
    ) -> None:
        run_vcf_scan(vcf=indexed_vcf, cohort=small_cohort, cfg=default_cfg, outdir=outdir)
        rows = _read_tsv(outdir / "hits.tsv")
        classes = [r["strictness_class"] for r in rows]
        assert classes.count("strict_complete") == 4        # pos=100, 700, 800/T, 900
        assert classes.count("strict_target_missing") == 1  # pos=200
        assert classes.count("strict_offtarget_missing") == 1  # pos=300
        assert classes.count("strict_both_missing") == 1    # pos=400

    def test_hits_sorted_by_descending_score(
        self, indexed_vcf, small_cohort, default_cfg, outdir
    ) -> None:
        run_vcf_scan(vcf=indexed_vcf, cohort=small_cohort, cfg=default_cfg, outdir=outdir)
        rows = _read_tsv(outdir / "hits.tsv")
        scores = [float(r["final_score"]) for r in rows]
        assert scores == sorted(scores, reverse=True), (
            "hits.tsv rows must be ordered by descending final_score (rank order)"
        )

    def test_locus_ids_are_sequential_ppx(
        self, indexed_vcf, small_cohort, default_cfg, outdir
    ) -> None:
        run_vcf_scan(vcf=indexed_vcf, cohort=small_cohort, cfg=default_cfg, outdir=outdir)
        rows = _read_tsv(outdir / "hits.tsv")
        ids = [r["locus_id"] for r in rows]
        assert all(lid.startswith("PPX") for lid in ids)
        assert len(ids) == len(set(ids)), "locus_ids must be unique"

    def test_indel_variant_type_at_pos_900(
        self, indexed_vcf, small_cohort, default_cfg, outdir
    ) -> None:
        run_vcf_scan(vcf=indexed_vcf, cohort=small_cohort, cfg=default_cfg, outdir=outdir)
        rows = _read_tsv(outdir / "hits.tsv")
        types_by_key = {r["allele_key"]: r["variant_type"] for r in rows}
        indel_key = next((k for k in types_by_key if "900" in k), None)
        assert indel_key is not None, "Expected a hit at pos=900 (indel)"
        assert types_by_key[indel_key] == "indel"

    def test_allele_key_uses_1based_pos(
        self, indexed_vcf, small_cohort, default_cfg, outdir
    ) -> None:
        run_vcf_scan(vcf=indexed_vcf, cohort=small_cohort, cfg=default_cfg, outdir=outdir)
        rows = _read_tsv(outdir / "hits.tsv")
        keys = [r["allele_key"] for r in rows]
        # allele_key format is contig:pos:ref:alt (1-based VCF POS)
        assert any(k.startswith("chr1:100:") for k in keys)

    def test_locus_coordinates_are_0based_halfopen(
        self, indexed_vcf, small_cohort, default_cfg, outdir
    ) -> None:
        run_vcf_scan(vcf=indexed_vcf, cohort=small_cohort, cfg=default_cfg, outdir=outdir)
        rows = _read_tsv(outdir / "hits.tsv")
        # pos=100 SNP A>T → locus_start = 99, locus_end = 100
        hit_100 = next(r for r in rows if "chr1:100:" in r["allele_key"])
        assert int(hit_100["start"]) == 99
        assert int(hit_100["end"]) == 100


# ---------------------------------------------------------------------------
# QC metrics
# ---------------------------------------------------------------------------


class TestQcMetrics:
    def test_records_evaluated_count(
        self, indexed_vcf, small_cohort, default_cfg, outdir
    ) -> None:
        run_vcf_scan(vcf=indexed_vcf, cohort=small_cohort, cfg=default_cfg, outdir=outdir)
        # 9 records total (pos=100 through 900)
        assert _qc(outdir)["records_evaluated"] == "9"

    def test_records_skipped_filter(
        self, indexed_vcf, small_cohort, default_cfg, outdir
    ) -> None:
        run_vcf_scan(vcf=indexed_vcf, cohort=small_cohort, cfg=default_cfg, outdir=outdir)
        # pos=600 has FILTER=FAIL
        assert _qc(outdir)["records_skipped_filter"] == "1"

    def test_alleles_contradicted_count(
        self, indexed_vcf, small_cohort, default_cfg, outdir
    ) -> None:
        run_vcf_scan(vcf=indexed_vcf, cohort=small_cohort, cfg=default_cfg, outdir=outdir)
        # pos=500: O1 carries the allele → contradicted
        assert _qc(outdir)["alleles_contradicted"] == "1"

    def test_alleles_passed_count(
        self, indexed_vcf, small_cohort, default_cfg, outdir
    ) -> None:
        run_vcf_scan(vcf=indexed_vcf, cohort=small_cohort, cfg=default_cfg, outdir=outdir)
        assert _qc(outdir)["alleles_passed"] == "7"

    def test_n_contigs_scanned(
        self, indexed_vcf, small_cohort, default_cfg, outdir
    ) -> None:
        run_vcf_scan(vcf=indexed_vcf, cohort=small_cohort, cfg=default_cfg, outdir=outdir)
        assert _qc(outdir)["n_contigs_scanned"] == "1"

    def test_n_target_and_offtarget_samples(
        self, indexed_vcf, small_cohort, default_cfg, outdir
    ) -> None:
        run_vcf_scan(vcf=indexed_vcf, cohort=small_cohort, cfg=default_cfg, outdir=outdir)
        qc = _qc(outdir)
        assert qc["n_target_samples"] == "2"
        assert qc["n_offtarget_samples"] == "3"

    def test_min_qual_skipping(
        self, indexed_vcf, small_cohort, default_cfg, outdir
    ) -> None:
        """With min_qual=10, pos=700 (QUAL=5) is skipped → 6 hits."""
        cfg = default_cfg.model_copy(
            update={"scan": default_cfg.scan.model_copy(update={"min_qual": 10.0})}
        )
        run_vcf_scan(vcf=indexed_vcf, cohort=small_cohort, cfg=cfg, outdir=outdir)
        assert _qc(outdir)["records_skipped_qual"] == "1"
        assert len(_read_tsv(outdir / "hits.tsv")) == 6

    def test_no_allow_multiallelic_skipping(
        self, indexed_vcf, small_cohort, default_cfg, outdir
    ) -> None:
        """With allow_multiallelic=False, pos=800 (2 ALTs) is skipped → 6 hits."""
        cfg = default_cfg.model_copy(
            update={
                "scan": default_cfg.scan.model_copy(
                    update={"allow_multiallelic": False}
                )
            }
        )
        run_vcf_scan(vcf=indexed_vcf, cohort=small_cohort, cfg=cfg, outdir=outdir)
        assert _qc(outdir)["records_skipped_multiallelic"] == "1"
        assert len(_read_tsv(outdir / "hits.tsv")) == 6


# ---------------------------------------------------------------------------
# regions.tsv
# ---------------------------------------------------------------------------


class TestRegionsTsv:
    def test_regions_written_with_default_merge_distance(
        self, indexed_vcf, small_cohort, default_cfg, outdir
    ) -> None:
        """merge_distance=0 (default): each locus becomes its own region."""
        run_vcf_scan(vcf=indexed_vcf, cohort=small_cohort, cfg=default_cfg, outdir=outdir)
        rows = _read_tsv(outdir / "regions.tsv")
        assert len(rows) == 7

    def test_regions_columns(
        self, indexed_vcf, small_cohort, default_cfg, outdir
    ) -> None:
        run_vcf_scan(vcf=indexed_vcf, cohort=small_cohort, cfg=default_cfg, outdir=outdir)
        rows = _read_tsv(outdir / "regions.tsv")
        assert rows
        expected = {
            "region_id", "contig", "start", "end", "n_loci",
            "variant_types", "dominant_strictness_class",
            "target_consistency", "offtarget_exclusion", "final_score",
        }
        assert set(rows[0].keys()) == expected

    def test_merge_distance_collapses_nearby_loci(
        self, indexed_vcf, small_cohort, default_cfg, outdir
    ) -> None:
        """merge_distance=500 should merge all 7 hits on chr1 into one region."""
        cfg = default_cfg.model_copy(
            update={"scan": default_cfg.scan.model_copy(update={"merge_distance": 500})}
        )
        run_vcf_scan(vcf=indexed_vcf, cohort=small_cohort, cfg=cfg, outdir=outdir)
        rows = _read_tsv(outdir / "regions.tsv")
        assert len(rows) == 1
        assert int(rows[0]["n_loci"]) == 7


# ---------------------------------------------------------------------------
# sample_support.tsv
# ---------------------------------------------------------------------------


class TestSampleSupportTsv:
    def test_sample_support_has_rows(
        self, indexed_vcf, small_cohort, default_cfg, outdir
    ) -> None:
        run_vcf_scan(vcf=indexed_vcf, cohort=small_cohort, cfg=default_cfg, outdir=outdir)
        rows = _read_tsv(outdir / "sample_support.tsv")
        # 7 hits × 5 samples = 35 rows
        assert len(rows) == 35

    def test_depth_and_allele_fraction_are_na(
        self, indexed_vcf, small_cohort, default_cfg, outdir
    ) -> None:
        """BAM support is not implemented in Phase 2; columns must be 'NA'."""
        run_vcf_scan(vcf=indexed_vcf, cohort=small_cohort, cfg=default_cfg, outdir=outdir)
        rows = _read_tsv(outdir / "sample_support.tsv")
        assert all(r["depth"] == "NA" for r in rows)
        assert all(r["allele_fraction"] == "NA" for r in rows)


# ---------------------------------------------------------------------------
# Error handling
# ---------------------------------------------------------------------------


class TestErrorHandling:
    def test_missing_vcf_raises_file_not_found(
        self, small_cohort, default_cfg, tmp_path
    ) -> None:
        with pytest.raises(FileNotFoundError):
            run_vcf_scan(
                vcf=tmp_path / "nonexistent.vcf.gz",
                cohort=small_cohort,
                cfg=default_cfg,
                outdir=tmp_path,
            )

    def test_no_target_samples_in_vcf_raises_value_error(
        self, indexed_vcf, default_cfg, tmp_path
    ) -> None:
        from privy.core.cohort import CohortDefinition  # noqa: PLC0415

        ghost_cohort = CohortDefinition.from_lists(
            targets=["GHOST1", "GHOST2"],
            off_targets=["GHOST3", "GHOST4"],
        )
        with pytest.raises(ValueError, match="No target samples"):
            run_vcf_scan(
                vcf=indexed_vcf,
                cohort=ghost_cohort,
                cfg=default_cfg,
                outdir=tmp_path,
            )

    def test_unsupported_mode_raises_not_implemented(
        self, indexed_vcf, small_cohort, default_cfg, tmp_path
    ) -> None:
        with pytest.raises(NotImplementedError):
            run_vcf_scan(
                vcf=indexed_vcf,
                cohort=small_cohort,
                cfg=default_cfg,
                outdir=tmp_path,
                mode="private_sv_state",
            )
