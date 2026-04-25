"""Integration tests for the standalone GFA scan backend.

Exercises run_gfa_scan() end-to-end against the small_cohort.gfa fixture and
verifies output file schema, hit counts, strictness classification, region
merging, and error handling.

Expected results with default config (min_target_support=1.0,
max_off_target_support=0.0, min_segment_length=1):
    - 2 hits: GPX00000001 (s2_target, strict_complete)
               GPX00000002 (s4_target, strict_target_missing)
    - 2 regions (merge_distance=0 → each hit is its own region)
"""

from __future__ import annotations

from pathlib import Path

import pytest

from privy.backends.gfa_scan import run_gfa_scan
from privy.core.cohort import CohortDefinition
from privy.core.config import default_config
from privy.io.tsv import read_tsv

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

GFA_PATH = Path(__file__).parent.parent / "data" / "small_cohort.gfa"

TARGETS = ["T1", "T2"]
OFF_TARGETS = ["O1", "O2", "O3"]


@pytest.fixture
def cohort() -> CohortDefinition:
    return CohortDefinition.from_lists(targets=TARGETS, off_targets=OFF_TARGETS)


@pytest.fixture
def cfg():
    return default_config()


@pytest.fixture
def scan_outdir(tmp_path: Path, cohort: CohortDefinition, cfg) -> Path:
    """Run a default scan and return the output directory."""
    run_gfa_scan(
        gfa=GFA_PATH,
        cohort=cohort,
        cfg=cfg,
        outdir=tmp_path,
    )
    return tmp_path


# ---------------------------------------------------------------------------
# TestOutputFiles
# ---------------------------------------------------------------------------


class TestOutputFiles:
    def test_hits_tsv_created(self, scan_outdir: Path) -> None:
        assert (scan_outdir / "hits.tsv").exists()

    def test_regions_tsv_created(self, scan_outdir: Path) -> None:
        assert (scan_outdir / "regions.tsv").exists()

    def test_evidence_tsv_created(self, scan_outdir: Path) -> None:
        assert (scan_outdir / "evidence.tsv").exists()

    def test_sample_support_tsv_created(self, scan_outdir: Path) -> None:
        assert (scan_outdir / "sample_support.tsv").exists()

    def test_qc_tsv_created(self, scan_outdir: Path) -> None:
        assert (scan_outdir / "qc.tsv").exists()

    def test_run_json_created(self, scan_outdir: Path) -> None:
        import json
        p = scan_outdir / "run.json"
        assert p.exists()
        data = json.loads(p.read_text())
        assert "inputs" in data
        assert "gfa" in data["inputs"]


# ---------------------------------------------------------------------------
# TestHitsTsv
# ---------------------------------------------------------------------------


class TestHitsTsv:
    def test_default_hit_count(self, scan_outdir: Path) -> None:
        rows = read_tsv(scan_outdir / "hits.tsv")
        assert len(rows) == 2

    def test_hit_locus_ids(self, scan_outdir: Path) -> None:
        rows = read_tsv(scan_outdir / "hits.tsv")
        ids = {r["locus_id"] for r in rows}
        assert ids == {"GPX00000001", "GPX00000002"}

    def test_s2_target_is_strict_complete(self, scan_outdir: Path) -> None:
        rows = read_tsv(scan_outdir / "hits.tsv")
        # s2_target sits at chr1:8-18 → start=8
        s2_row = next(r for r in rows if r["start"] == "8")
        assert s2_row["strictness_class"] == "strict_complete"
        assert s2_row["target_support_n"] == "2"
        assert s2_row["offtarget_support_n"] == "0"
        assert s2_row["target_missing_n"] == "0"

    def test_s4_target_is_strict_target_missing(self, scan_outdir: Path) -> None:
        rows = read_tsv(scan_outdir / "hits.tsv")
        # s4_target sits at chr1:60-67 → start=60
        s4_row = next(r for r in rows if r["start"] == "60")
        assert s4_row["strictness_class"] == "strict_target_missing"
        assert s4_row["target_support_n"] == "1"
        assert s4_row["target_missing_n"] == "1"
        assert s4_row["offtarget_support_n"] == "0"

    def test_variant_type_is_graph_region(self, scan_outdir: Path) -> None:
        rows = read_tsv(scan_outdir / "hits.tsv")
        for r in rows:
            assert r["variant_type"] == "graph_region"

    def test_allele_key_contains_seg_prefix(self, scan_outdir: Path) -> None:
        rows = read_tsv(scan_outdir / "hits.tsv")
        for r in rows:
            assert "SEG:" in r["allele_key"]

    def test_scores_are_numbers(self, scan_outdir: Path) -> None:
        rows = read_tsv(scan_outdir / "hits.tsv")
        for r in rows:
            float(r["discovery_score"])
            float(r["final_score"])

    def test_strict_complete_has_higher_score(self, scan_outdir: Path) -> None:
        rows = read_tsv(scan_outdir / "hits.tsv")
        s2_row = next(r for r in rows if r["start"] == "8")
        s4_row = next(r for r in rows if r["start"] == "60")
        assert float(s2_row["final_score"]) >= float(s4_row["final_score"])

    def test_locus_coordinates_0based(self, scan_outdir: Path) -> None:
        rows = read_tsv(scan_outdir / "hits.tsv")
        s2_row = next(r for r in rows if r["start"] == "8")
        # s2_target: SO=8, LN=10 → start=8, end=18
        assert s2_row["start"] == "8"
        assert s2_row["end"] == "18"

    def test_backbone_segments_not_in_hits(self, scan_outdir: Path) -> None:
        rows = read_tsv(scan_outdir / "hits.tsv")
        allele_keys = {r["allele_key"] for r in rows}
        for seg in ("s1", "s3", "s5", "s2_offt", "s4_offt"):
            assert not any(seg in k for k in allele_keys)


# ---------------------------------------------------------------------------
# TestRegionsTsv
# ---------------------------------------------------------------------------


class TestRegionsTsv:
    def test_default_region_count(self, scan_outdir: Path) -> None:
        rows = read_tsv(scan_outdir / "regions.tsv")
        # merge_distance=0, two hits at different positions → 2 regions
        assert len(rows) == 2

    def test_region_variant_type(self, scan_outdir: Path) -> None:
        rows = read_tsv(scan_outdir / "regions.tsv")
        for r in rows:
            assert r["variant_types"] == "graph_region"

    def test_merge_distance_collapses_hits(
        self, tmp_path: Path, cohort: CohortDefinition, cfg
    ) -> None:
        merged_dir = tmp_path / "merged"
        merged_dir.mkdir()
        merged_cfg = cfg.model_copy(
            update={"scan": cfg.scan.model_copy(update={"merge_distance": 10000})}
        )
        run_gfa_scan(
            gfa=GFA_PATH,
            cohort=cohort,
            cfg=merged_cfg,
            outdir=merged_dir,
        )
        rows = read_tsv(merged_dir / "regions.tsv")
        assert len(rows) == 1


# ---------------------------------------------------------------------------
# TestQcMetrics
# ---------------------------------------------------------------------------


class TestQcMetrics:
    def test_qc_has_required_metrics(self, scan_outdir: Path) -> None:
        rows = read_tsv(scan_outdir / "qc.tsv")
        metrics = {r["metric"] for r in rows}
        assert "alleles_passed" in metrics
        assert "alleles_contradicted" in metrics
        assert "loci_emitted" in metrics
        assert "regions_emitted" in metrics

    def test_loci_emitted_equals_hit_count(self, scan_outdir: Path) -> None:
        qc_rows = read_tsv(scan_outdir / "qc.tsv")
        hits = read_tsv(scan_outdir / "hits.tsv")
        loci_qc = next(r for r in qc_rows if r["metric"] == "loci_emitted")
        assert int(loci_qc["value"]) == len(hits)

    def test_strictness_counts_present(self, scan_outdir: Path) -> None:
        rows = read_tsv(scan_outdir / "qc.tsv")
        metrics = {r["metric"] for r in rows}
        assert "strictness_strict_complete" in metrics
        assert "strictness_strict_target_missing" in metrics


# ---------------------------------------------------------------------------
# TestSampleSupportTsv
# ---------------------------------------------------------------------------


class TestSampleSupportTsv:
    def test_row_count(self, scan_outdir: Path) -> None:
        rows = read_tsv(scan_outdir / "sample_support.tsv")
        # 2 hits × 5 samples = 10 rows
        assert len(rows) == 10

    def test_t1_traverses_s2_target(self, scan_outdir: Path) -> None:
        rows = read_tsv(scan_outdir / "sample_support.tsv")
        # hit for s2_target is at chr1:8 → locus GPX00000001
        t1_s2 = next(
            r for r in rows
            if r["sample_id"] == "T1" and r["locus_id"] == "GPX00000001"
        )
        assert t1_s2["genotype"] == "traverses"
        assert t1_s2["allele_supported"] == "true"
        assert t1_s2["evidence_class"] == "support"

    def test_o1_absent_from_s2_target(self, scan_outdir: Path) -> None:
        rows = read_tsv(scan_outdir / "sample_support.tsv")
        o1_s2 = next(
            r for r in rows
            if r["sample_id"] == "O1" and r["locus_id"] == "GPX00000001"
        )
        assert o1_s2["genotype"] == "absent"
        assert o1_s2["allele_supported"] == "false"
        assert o1_s2["evidence_class"] == "absence"

    def test_t2_missing_from_s4_target(self, scan_outdir: Path) -> None:
        rows = read_tsv(scan_outdir / "sample_support.tsv")
        t2_s4 = next(
            r for r in rows
            if r["sample_id"] == "T2" and r["locus_id"] == "GPX00000002"
        )
        assert t2_s4["genotype"] == "missing"
        assert t2_s4["evidence_class"] == "uninformative"

    def test_depth_is_na(self, scan_outdir: Path) -> None:
        rows = read_tsv(scan_outdir / "sample_support.tsv")
        for r in rows:
            assert r["depth"] == "NA"
            assert r["allele_fraction"] == "NA"


# ---------------------------------------------------------------------------
# TestMinSegmentLength
# ---------------------------------------------------------------------------


class TestMinSegmentLength:
    def test_min_segment_length_7_excludes_s2_target(
        self, tmp_path: Path, cohort: CohortDefinition, cfg
    ) -> None:
        """s2_target is 10 bp → still passes min_len=7.
        s4_target is 7 bp → exactly at threshold, still passes."""
        out = tmp_path / "minlen7"
        out.mkdir()
        new_cfg = cfg.model_copy(
            update={"gfa": cfg.gfa.model_copy(update={"min_segment_length": 7})}
        )
        run_gfa_scan(gfa=GFA_PATH, cohort=cohort, cfg=new_cfg, outdir=out)
        rows = read_tsv(out / "hits.tsv")
        assert len(rows) == 2   # both still pass

    def test_min_segment_length_11_excludes_s4_target(
        self, tmp_path: Path, cohort: CohortDefinition, cfg
    ) -> None:
        """s4_target is 7 bp → excluded at min_len=11.  Only s2_target (10 bp) remains.
        Wait — s2_target is 10 bp, also excluded at 11.  So 0 hits."""
        out = tmp_path / "minlen11"
        out.mkdir()
        new_cfg = cfg.model_copy(
            update={"gfa": cfg.gfa.model_copy(update={"min_segment_length": 11})}
        )
        run_gfa_scan(gfa=GFA_PATH, cohort=cohort, cfg=new_cfg, outdir=out)
        rows = read_tsv(out / "hits.tsv")
        assert len(rows) == 0


# ---------------------------------------------------------------------------
# TestContigFilter
# ---------------------------------------------------------------------------


class TestContigFilter:
    def test_contig_filter_chr1_gives_all_hits(
        self, tmp_path: Path, cohort: CohortDefinition, cfg
    ) -> None:
        out = tmp_path / "ctg"
        out.mkdir()
        run_gfa_scan(gfa=GFA_PATH, cohort=cohort, cfg=cfg, outdir=out, contig="chr1")
        rows = read_tsv(out / "hits.tsv")
        assert len(rows) == 2

    def test_contig_filter_chrx_gives_no_hits(
        self, tmp_path: Path, cohort: CohortDefinition, cfg
    ) -> None:
        out = tmp_path / "noctg"
        out.mkdir()
        run_gfa_scan(gfa=GFA_PATH, cohort=cohort, cfg=cfg, outdir=out, contig="chrX")
        rows = read_tsv(out / "hits.tsv")
        assert len(rows) == 0

    def test_region_filter_bubble1_only(
        self, tmp_path: Path, cohort: CohortDefinition, cfg
    ) -> None:
        out = tmp_path / "reg"
        out.mkdir()
        run_gfa_scan(
            gfa=GFA_PATH, cohort=cohort, cfg=cfg, outdir=out,
            region="chr1:1-30",   # covers chr1:0-30 (0-based after conversion)
        )
        rows = read_tsv(out / "hits.tsv")
        assert len(rows) == 1
        assert rows[0]["start"] == "8"


# ---------------------------------------------------------------------------
# TestErrorHandling
# ---------------------------------------------------------------------------


class TestErrorHandling:
    def test_file_not_found(self, tmp_path: Path, cohort: CohortDefinition, cfg) -> None:
        with pytest.raises(FileNotFoundError):
            run_gfa_scan(
                gfa=tmp_path / "missing.gfa",
                cohort=cohort,
                cfg=cfg,
                outdir=tmp_path,
            )

    def test_no_target_samples_in_gfa(self, tmp_path: Path, cfg) -> None:
        wrong_cohort = CohortDefinition.from_lists(
            targets=["GHOST1", "GHOST2"],
            off_targets=["O1", "O2", "O3"],
        )
        with pytest.raises(ValueError, match="No target samples"):
            run_gfa_scan(
                gfa=GFA_PATH,
                cohort=wrong_cohort,
                cfg=cfg,
                outdir=tmp_path,
            )

    def test_no_offtarget_samples_in_gfa(self, tmp_path: Path, cfg) -> None:
        wrong_cohort = CohortDefinition.from_lists(
            targets=["T1", "T2"],
            off_targets=["GHOST3"],
        )
        with pytest.raises(ValueError, match="No off-target samples"):
            run_gfa_scan(
                gfa=GFA_PATH,
                cohort=wrong_cohort,
                cfg=cfg,
                outdir=tmp_path,
            )

    def test_unsupported_mode_raises(
        self, tmp_path: Path, cohort: CohortDefinition, cfg
    ) -> None:
        with pytest.raises(NotImplementedError):
            run_gfa_scan(
                gfa=GFA_PATH,
                cohort=cohort,
                cfg=cfg,
                outdir=tmp_path,
                mode="private_genotype",
            )

    def test_write_flags_suppress_files(
        self, tmp_path: Path, cohort: CohortDefinition, cfg
    ) -> None:
        out = tmp_path / "minimal"
        out.mkdir()
        run_gfa_scan(
            gfa=GFA_PATH, cohort=cohort, cfg=cfg, outdir=out,
            write_hits=False,
            write_regions=False,
            write_evidence=False,
            write_sample_support=False,
            write_qc=False,
            write_run_json=False,
        )
        assert not (out / "hits.tsv").exists()
        assert not (out / "run.json").exists()
