"""Integration tests for the BAM support layer.

Covers:
  - resolve_bam_sample_pairs()
  - annotate_loci_with_bam()
  - _classify_bam_evidence() (via annotate_loci_with_bam)
  - End-to-end: run_vcf_scan() with BAM inputs
"""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from privy.backends.bam_support import (
    BamAnnotationResult,
    HitLocusInfo,
    _classify_bam_evidence,
    annotate_loci_with_bam,
    resolve_bam_sample_pairs,
)
from privy.core.cohort import CohortDefinition
from privy.core.config import BamConfig, default_config
from privy.core.evidence import EvidenceClass


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_cfg(
    min_depth: int = 8,
    min_alt_count: int = 2,
    allele_fraction_min: float = 0.2,
    min_mapq: int = 20,
    min_baseq: int = 20,
) -> BamConfig:
    return BamConfig(
        min_depth=min_depth,
        min_alt_count=min_alt_count,
        allele_fraction_min=allele_fraction_min,
        min_mapq=min_mapq,
        min_baseq=min_baseq,
    )


def _snp_locus(
    locus_id: str = "PPX00000001",
    contig: str = "chr1",
    start: int = 99,
    end: int = 100,
) -> HitLocusInfo:
    return HitLocusInfo(
        locus_id=locus_id,
        contig=contig,
        start=start,
        end=end,
        variant_type="snp",
        ref_allele="A",
        alt_allele="T",
    )


def _indel_locus(locus_id: str = "PPX00000002") -> HitLocusInfo:
    return HitLocusInfo(
        locus_id=locus_id,
        contig="chr1",
        start=898,
        end=901,
        variant_type="deletion",
        ref_allele="AGG",
        alt_allele="A",
    )


@pytest.fixture
def small_cohort_t1_o1() -> CohortDefinition:
    return CohortDefinition.from_lists(targets=["T1"], off_targets=["O1"])


# ---------------------------------------------------------------------------
# resolve_bam_sample_pairs
# ---------------------------------------------------------------------------

class TestResolveBamSamplePairs:
    def test_resolves_from_bam_paths(
        self, bam_target_t1: Path, small_cohort_t1_o1: CohortDefinition
    ) -> None:
        pairs = resolve_bam_sample_pairs([bam_target_t1], None, small_cohort_t1_o1)
        assert len(pairs) == 1
        path, sample = pairs[0]
        assert path == bam_target_t1
        assert sample == "T1"

    def test_resolves_from_manifest(
        self,
        bam_target_t1: Path,
        bam_offtarget_o1: Path,
        tmp_path: Path,
        small_cohort_t1_o1: CohortDefinition,
    ) -> None:
        manifest = tmp_path / "manifest.tsv"
        manifest.write_text(
            f"bam_path\tsample_id\n{bam_target_t1}\tT1\n{bam_offtarget_o1}\tO1\n"
        )
        pairs = resolve_bam_sample_pairs(None, manifest, small_cohort_t1_o1)
        assert len(pairs) == 2
        samples = [s for _, s in pairs]
        assert "T1" in samples and "O1" in samples

    def test_manifest_takes_precedence_over_paths(
        self,
        bam_target_t1: Path,
        bam_offtarget_o1: Path,
        tmp_path: Path,
        small_cohort_t1_o1: CohortDefinition,
    ) -> None:
        manifest = tmp_path / "manifest.tsv"
        manifest.write_text(
            f"bam_path\tsample_id\n{bam_offtarget_o1}\tO1\n"
        )
        pairs = resolve_bam_sample_pairs(
            [bam_target_t1], manifest, small_cohort_t1_o1
        )
        assert len(pairs) == 1
        assert pairs[0][1] == "O1"

    def test_falls_back_to_filename_stem_when_no_sm(
        self, tmp_path: Path, small_cohort_t1_o1: CohortDefinition
    ) -> None:
        import pysam  # noqa: PLC0415

        bam_path = tmp_path / "MySample.bam"
        header = pysam.AlignmentHeader.from_dict({
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": "chr1", "LN": 1000}],
        })
        with pysam.AlignmentFile(str(bam_path), "wb", header=header):
            pass
        pysam.index(str(bam_path))
        pairs = resolve_bam_sample_pairs([bam_path], None, small_cohort_t1_o1)
        assert pairs[0][1] == "MySample"

    def test_empty_inputs_return_empty(
        self, small_cohort_t1_o1: CohortDefinition
    ) -> None:
        assert resolve_bam_sample_pairs(None, None, small_cohort_t1_o1) == []
        assert resolve_bam_sample_pairs([], None, small_cohort_t1_o1) == []


# ---------------------------------------------------------------------------
# _classify_bam_evidence
# ---------------------------------------------------------------------------

class TestClassifyBamEvidence:
    def _cfg(self, **kwargs) -> BamConfig:
        return _make_cfg(**kwargs)

    def test_low_depth_is_uninformative(self) -> None:
        cls, val = _classify_bam_evidence(
            ref_count=3, alt_count=2, depth=5,
            allele_fraction=0.4, cohort_role="target",
            cfg=_make_cfg(min_depth=8), is_snp=True,
        )
        assert cls == EvidenceClass.UNINFORMATIVE
        assert val == 0.0

    def test_indel_is_uninformative_regardless_of_depth(self) -> None:
        cls, val = _classify_bam_evidence(
            ref_count=12, alt_count=0, depth=12,
            allele_fraction=None, cohort_role="target",
            cfg=_make_cfg(min_depth=8), is_snp=False,
        )
        assert cls == EvidenceClass.UNINFORMATIVE
        assert val == 0.0

    def test_target_with_high_alt_is_support(self) -> None:
        cls, val = _classify_bam_evidence(
            ref_count=0, alt_count=12, depth=12,
            allele_fraction=1.0, cohort_role="target",
            cfg=_make_cfg(min_depth=8, min_alt_count=2, allele_fraction_min=0.2),
            is_snp=True,
        )
        assert cls == EvidenceClass.SUPPORT
        assert val == 1.0

    def test_target_with_no_alt_is_ambiguous(self) -> None:
        cls, val = _classify_bam_evidence(
            ref_count=12, alt_count=0, depth=12,
            allele_fraction=0.0, cohort_role="target",
            cfg=_make_cfg(min_depth=8), is_snp=True,
        )
        assert cls == EvidenceClass.AMBIGUOUS
        assert val == 0.3

    def test_offtarget_with_no_alt_is_absence(self) -> None:
        cls, val = _classify_bam_evidence(
            ref_count=12, alt_count=0, depth=12,
            allele_fraction=0.0, cohort_role="off_target",
            cfg=_make_cfg(min_depth=8), is_snp=True,
        )
        assert cls == EvidenceClass.ABSENCE
        assert val == 1.0

    def test_offtarget_with_alt_is_contradiction(self) -> None:
        cls, val = _classify_bam_evidence(
            ref_count=0, alt_count=12, depth=12,
            allele_fraction=1.0, cohort_role="off_target",
            cfg=_make_cfg(min_depth=8, min_alt_count=2, allele_fraction_min=0.2),
            is_snp=True,
        )
        assert cls == EvidenceClass.CONTRADICTION
        assert val == -1.0

    def test_alt_count_below_threshold_treated_as_absent(self) -> None:
        cls, _ = _classify_bam_evidence(
            ref_count=10, alt_count=1, depth=11,
            allele_fraction=0.09, cohort_role="off_target",
            cfg=_make_cfg(min_depth=8, min_alt_count=2, allele_fraction_min=0.2),
            is_snp=True,
        )
        assert cls == EvidenceClass.ABSENCE


# ---------------------------------------------------------------------------
# annotate_loci_with_bam — evidence classification
# ---------------------------------------------------------------------------

class TestAnnotateLociEvidenceClasses:
    def test_target_bam_produces_support(
        self,
        bam_target_t1: Path,
        small_cohort_t1_o1: CohortDefinition,
    ) -> None:
        loci = [_snp_locus()]
        result = annotate_loci_with_bam(
            loci, [(bam_target_t1, "T1")], small_cohort_t1_o1, _make_cfg()
        )
        classes = {er.evidence_class for er in result.evidence_records}
        assert EvidenceClass.SUPPORT in classes

    def test_offtarget_clean_produces_absence(
        self,
        bam_offtarget_o1: Path,
        small_cohort_t1_o1: CohortDefinition,
    ) -> None:
        loci = [_snp_locus()]
        result = annotate_loci_with_bam(
            loci, [(bam_offtarget_o1, "O1")], small_cohort_t1_o1, _make_cfg()
        )
        classes = {er.evidence_class for er in result.evidence_records}
        assert EvidenceClass.ABSENCE in classes

    def test_offtarget_with_alt_produces_contradiction(
        self,
        bam_offtarget_o1_with_alt: Path,
        small_cohort_t1_o1: CohortDefinition,
    ) -> None:
        loci = [_snp_locus()]
        result = annotate_loci_with_bam(
            loci,
            [(bam_offtarget_o1_with_alt, "O1")],
            small_cohort_t1_o1,
            _make_cfg(),
        )
        classes = {er.evidence_class for er in result.evidence_records}
        assert EvidenceClass.CONTRADICTION in classes

    def test_low_depth_produces_uninformative(
        self,
        bam_low_depth_t2: Path,
        small_cohort_t1_o1: CohortDefinition,
    ) -> None:
        cohort = CohortDefinition.from_lists(targets=["T2"], off_targets=["O1"])
        loci = [_snp_locus()]
        result = annotate_loci_with_bam(
            loci, [(bam_low_depth_t2, "T2")], cohort, _make_cfg(min_depth=8)
        )
        classes = {er.evidence_class for er in result.evidence_records}
        assert EvidenceClass.UNINFORMATIVE in classes

    def test_indel_locus_produces_uninformative(
        self,
        bam_target_t1: Path,
        small_cohort_t1_o1: CohortDefinition,
    ) -> None:
        loci = [_indel_locus()]
        result = annotate_loci_with_bam(
            loci, [(bam_target_t1, "T1")], small_cohort_t1_o1, _make_cfg()
        )
        assert all(
            er.evidence_class == EvidenceClass.UNINFORMATIVE
            for er in result.evidence_records
        )

    def test_missing_bam_index_is_skipped(
        self,
        tmp_path: Path,
        small_cohort_t1_o1: CohortDefinition,
    ) -> None:
        import pysam  # noqa: PLC0415

        bam_path = tmp_path / "noindex.bam"
        header = pysam.AlignmentHeader.from_dict({
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": "chr1", "LN": 1000}],
            "RG": [{"ID": "T1", "SM": "T1"}],
        })
        with pysam.AlignmentFile(str(bam_path), "wb", header=header):
            pass
        # No index created
        loci = [_snp_locus()]
        result = annotate_loci_with_bam(
            loci, [(bam_path, "T1")], small_cohort_t1_o1, _make_cfg()
        )
        assert result.evidence_records == []

    def test_ignored_sample_produces_no_records(
        self,
        bam_target_t1: Path,
        tmp_path: Path,
    ) -> None:
        cohort = CohortDefinition.from_lists(
            targets=["other"],
            off_targets=["O1"],
            ignored_samples=["T1"],
        )
        loci = [_snp_locus()]
        result = annotate_loci_with_bam(
            loci, [(bam_target_t1, "T1")], cohort, _make_cfg()
        )
        assert result.evidence_records == []


# ---------------------------------------------------------------------------
# annotate_loci_with_bam — support scores
# ---------------------------------------------------------------------------

class TestAnnotateLociSupportScores:
    def test_support_score_positive_when_target_confirmed(
        self,
        bam_target_t1: Path,
        small_cohort_t1_o1: CohortDefinition,
    ) -> None:
        loci = [_snp_locus()]
        result = annotate_loci_with_bam(
            loci, [(bam_target_t1, "T1")], small_cohort_t1_o1, _make_cfg()
        )
        score = result.support_score_by_locus["PPX00000001"]
        assert score > 0.0
        assert score <= 1.0

    def test_support_score_zero_when_uninformative_only(
        self,
        bam_low_depth_t2: Path,
    ) -> None:
        cohort = CohortDefinition.from_lists(targets=["T2"], off_targets=["O1"])
        loci = [_snp_locus()]
        result = annotate_loci_with_bam(
            loci, [(bam_low_depth_t2, "T2")], cohort, _make_cfg(min_depth=8)
        )
        score = result.support_score_by_locus["PPX00000001"]
        assert score == 0.0

    def test_support_score_negative_when_contradiction(
        self,
        bam_offtarget_o1_with_alt: Path,
        small_cohort_t1_o1: CohortDefinition,
    ) -> None:
        loci = [_snp_locus()]
        result = annotate_loci_with_bam(
            loci,
            [(bam_offtarget_o1_with_alt, "O1")],
            small_cohort_t1_o1,
            _make_cfg(),
        )
        score = result.support_score_by_locus["PPX00000001"]
        assert score < 0.0

    def test_all_loci_have_score_entry(
        self,
        bam_target_t1: Path,
        small_cohort_t1_o1: CohortDefinition,
    ) -> None:
        loci = [_snp_locus("L1"), _snp_locus("L2"), _indel_locus("L3")]
        result = annotate_loci_with_bam(
            loci, [(bam_target_t1, "T1")], small_cohort_t1_o1, _make_cfg()
        )
        for locus in loci:
            assert locus.locus_id in result.support_score_by_locus

    def test_no_bam_pairs_returns_zero_scores(
        self, small_cohort_t1_o1: CohortDefinition
    ) -> None:
        loci = [_snp_locus()]
        result = annotate_loci_with_bam(loci, [], small_cohort_t1_o1, _make_cfg())
        assert result.support_score_by_locus["PPX00000001"] == 0.0
        assert result.evidence_records == []


# ---------------------------------------------------------------------------
# annotate_loci_with_bam — bam_metrics
# ---------------------------------------------------------------------------

class TestAnnotateLociBamMetrics:
    def test_metrics_populated_for_covered_snp(
        self,
        bam_target_t1: Path,
        small_cohort_t1_o1: CohortDefinition,
    ) -> None:
        loci = [_snp_locus()]
        result = annotate_loci_with_bam(
            loci, [(bam_target_t1, "T1")], small_cohort_t1_o1, _make_cfg()
        )
        m = result.bam_metrics[("PPX00000001", "T1")]
        assert m["depth"] == "12"
        assert m["allele_fraction"] == "1.0"

    def test_metrics_depth_zero_for_uncovered_region(
        self,
        bam_target_t1: Path,
        small_cohort_t1_o1: CohortDefinition,
    ) -> None:
        locus = HitLocusInfo("L9", "chr1", 5000, 5001, "snp", "A", "T")
        result = annotate_loci_with_bam(
            [locus], [(bam_target_t1, "T1")], small_cohort_t1_o1, _make_cfg()
        )
        m = result.bam_metrics[("L9", "T1")]
        assert m["depth"] == "0"

    def test_indel_metrics_allele_fraction_is_na(
        self,
        bam_target_t1: Path,
        small_cohort_t1_o1: CohortDefinition,
    ) -> None:
        loci = [_indel_locus()]
        result = annotate_loci_with_bam(
            loci, [(bam_target_t1, "T1")], small_cohort_t1_o1, _make_cfg()
        )
        m = result.bam_metrics[("PPX00000002", "T1")]
        assert m["allele_fraction"] == "NA"

    def test_evidence_records_have_sample_id(
        self,
        bam_target_t1: Path,
        small_cohort_t1_o1: CohortDefinition,
    ) -> None:
        loci = [_snp_locus()]
        result = annotate_loci_with_bam(
            loci, [(bam_target_t1, "T1")], small_cohort_t1_o1, _make_cfg()
        )
        for er in result.evidence_records:
            assert er.sample_id == "T1"


# ---------------------------------------------------------------------------
# End-to-end: run_vcf_scan with BAM support
# ---------------------------------------------------------------------------

class TestVcfScanWithBam:
    def test_support_score_nonzero_with_target_bam(
        self,
        indexed_vcf: Path,
        small_cohort: CohortDefinition,
        bam_target_t1: Path,
        tmp_path: Path,
    ) -> None:
        from privy.backends.vcf_scan import run_vcf_scan  # noqa: PLC0415
        from privy.core.config import default_config  # noqa: PLC0415

        cfg = default_config()
        outdir = tmp_path / "out"
        outdir.mkdir()
        run_vcf_scan(
            vcf=indexed_vcf,
            cohort=small_cohort,
            cfg=cfg,
            outdir=outdir,
            bam=[bam_target_t1],
        )
        hits_path = outdir / "hits.tsv"
        assert hits_path.exists()
        with open(hits_path, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            rows = list(reader)

        # At least one hit at contig=chr1 pos 99 (VCF pos 100) should have
        # support_score > 0 because the target BAM has alt reads there
        scores = [float(r["support_score"]) for r in rows]
        assert any(s > 0.0 for s in scores)

    def test_evidence_tsv_has_bam_rows(
        self,
        indexed_vcf: Path,
        small_cohort: CohortDefinition,
        bam_target_t1: Path,
        tmp_path: Path,
    ) -> None:
        from privy.backends.vcf_scan import run_vcf_scan  # noqa: PLC0415

        outdir = tmp_path / "out"
        outdir.mkdir()
        run_vcf_scan(
            vcf=indexed_vcf,
            cohort=small_cohort,
            cfg=default_config(),
            outdir=outdir,
            bam=[bam_target_t1],
        )
        with open(outdir / "evidence.tsv", newline="") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        source_types = {r["source_type"] for r in rows}
        assert "bam" in source_types
        assert "vcf" in source_types

    def test_sample_support_tsv_has_depth_values(
        self,
        indexed_vcf: Path,
        small_cohort: CohortDefinition,
        bam_target_t1: Path,
        tmp_path: Path,
    ) -> None:
        from privy.backends.vcf_scan import run_vcf_scan  # noqa: PLC0415

        outdir = tmp_path / "out"
        outdir.mkdir()
        run_vcf_scan(
            vcf=indexed_vcf,
            cohort=small_cohort,
            cfg=default_config(),
            outdir=outdir,
            bam=[bam_target_t1],
        )
        with open(outdir / "sample_support.tsv", newline="") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        t1_rows = [r for r in rows if r["sample_id"] == "T1"]
        # T1 BAM covers position 99; at least one row should have non-NA depth
        depths = [r["depth"] for r in t1_rows]
        assert any(d != "NA" for d in depths)

    def test_scan_without_bam_unchanged(
        self,
        indexed_vcf: Path,
        small_cohort: CohortDefinition,
        tmp_path: Path,
    ) -> None:
        from privy.backends.vcf_scan import run_vcf_scan  # noqa: PLC0415

        outdir = tmp_path / "out"
        outdir.mkdir()
        run_vcf_scan(
            vcf=indexed_vcf,
            cohort=small_cohort,
            cfg=default_config(),
            outdir=outdir,
        )
        with open(outdir / "hits.tsv", newline="") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        # All support_scores should be 0 when no BAM provided
        scores = [float(r["support_score"]) for r in rows]
        assert all(s == 0.0 for s in scores)

    def test_contradiction_reduces_final_score(
        self,
        indexed_vcf: Path,
        small_cohort: CohortDefinition,
        bam_offtarget_o1_with_alt: Path,
        tmp_path: Path,
    ) -> None:
        from privy.backends.vcf_scan import run_vcf_scan  # noqa: PLC0415

        # Baseline without BAM
        outdir_base = tmp_path / "base"
        outdir_base.mkdir()
        run_vcf_scan(
            vcf=indexed_vcf,
            cohort=small_cohort,
            cfg=default_config(),
            outdir=outdir_base,
        )

        # With contradiction BAM
        outdir_bam = tmp_path / "bam"
        outdir_bam.mkdir()
        run_vcf_scan(
            vcf=indexed_vcf,
            cohort=small_cohort,
            cfg=default_config(),
            outdir=outdir_bam,
            bam=[bam_offtarget_o1_with_alt],
        )

        def _max_final(out: Path) -> float:
            with open(out / "hits.tsv", newline="") as fh:
                rows = list(csv.DictReader(fh, delimiter="\t"))
            return max(float(r["final_score"]) for r in rows) if rows else 0.0

        # At least the top score should be lower with the contradiction BAM
        assert _max_final(outdir_bam) <= _max_final(outdir_base)

    def test_scan_with_manifest(
        self,
        indexed_vcf: Path,
        small_cohort: CohortDefinition,
        bam_target_t1: Path,
        bam_offtarget_o1: Path,
        tmp_path: Path,
    ) -> None:
        from privy.backends.vcf_scan import run_vcf_scan  # noqa: PLC0415

        manifest = tmp_path / "manifest.tsv"
        manifest.write_text(
            f"bam_path\tsample_id\n{bam_target_t1}\tT1\n{bam_offtarget_o1}\tO1\n"
        )
        outdir = tmp_path / "out"
        outdir.mkdir()
        run_vcf_scan(
            vcf=indexed_vcf,
            cohort=small_cohort,
            cfg=default_config(),
            outdir=outdir,
            bam_manifest=manifest,
        )
        assert (outdir / "hits.tsv").exists()
        with open(outdir / "evidence.tsv", newline="") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        assert any(r["source_type"] == "bam" for r in rows)
