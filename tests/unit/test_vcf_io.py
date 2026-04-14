"""Unit tests for :mod:`privy.io.vcf` helper functions.

Pure-Python tests (``is_missing_genotype``, ``has_alt_allele``,
``classify_variant_type``, ``format_allele_key``) require no pysam
installation.  Tests that exercise VCF file access use the ``indexed_vcf``
fixture from ``conftest.py``.
"""

from __future__ import annotations

import pytest

from privy.io.vcf import (
    classify_variant_type,
    format_allele_key,
    has_alt_allele,
    is_missing_genotype,
)


# ---------------------------------------------------------------------------
# is_missing_genotype
# ---------------------------------------------------------------------------


class TestIsMissingGenotype:
    def test_none_is_missing(self) -> None:
        assert is_missing_genotype(None) is True

    def test_both_alleles_none(self) -> None:
        assert is_missing_genotype((None, None)) is True

    def test_first_allele_none(self) -> None:
        assert is_missing_genotype((None, 0)) is True

    def test_second_allele_none(self) -> None:
        assert is_missing_genotype((0, None)) is True

    def test_diploid_ref_not_missing(self) -> None:
        assert is_missing_genotype((0, 0)) is False

    def test_het_not_missing(self) -> None:
        assert is_missing_genotype((0, 1)) is False

    def test_hom_alt_not_missing(self) -> None:
        assert is_missing_genotype((1, 1)) is False

    def test_haploid_ref_not_missing(self) -> None:
        assert is_missing_genotype((0,)) is False

    def test_haploid_alt_not_missing(self) -> None:
        assert is_missing_genotype((1,)) is False


# ---------------------------------------------------------------------------
# has_alt_allele
# ---------------------------------------------------------------------------


class TestHasAltAllele:
    def test_het_first_alt_present(self) -> None:
        # GT 0/1 — alt_index=0 means allele value 1
        assert has_alt_allele((0, 1), 0) is True

    def test_hom_first_alt_present(self) -> None:
        assert has_alt_allele((1, 1), 0) is True

    def test_ref_hom_does_not_have_alt(self) -> None:
        assert has_alt_allele((0, 0), 0) is False

    def test_second_alt_present(self) -> None:
        # GT 0/2 — alt_index=1 means allele value 2
        assert has_alt_allele((0, 2), 1) is True

    def test_first_alt_present_does_not_satisfy_second(self) -> None:
        # GT 0/1 — alt_index=1 checks for allele value 2, which is absent
        assert has_alt_allele((0, 1), 1) is False

    def test_haploid_alt_present(self) -> None:
        assert has_alt_allele((1,), 0) is True

    def test_haploid_ref_absent(self) -> None:
        assert has_alt_allele((0,), 0) is False


# ---------------------------------------------------------------------------
# classify_variant_type
# ---------------------------------------------------------------------------


class TestClassifyVariantType:
    def test_snp(self) -> None:
        assert classify_variant_type("A", "T") == "snp"

    def test_snp_gc(self) -> None:
        assert classify_variant_type("G", "C") == "snp"

    def test_deletion(self) -> None:
        assert classify_variant_type("AGG", "A") == "indel"

    def test_insertion(self) -> None:
        assert classify_variant_type("A", "ATG") == "indel"

    def test_mnp_is_indel(self) -> None:
        # Both REF and ALT are len > 1 — classified as indel by our rule
        assert classify_variant_type("AC", "GT") == "indel"

    def test_symbolic_del(self) -> None:
        assert classify_variant_type("A", "<DEL>") == "sv"

    def test_symbolic_ins(self) -> None:
        assert classify_variant_type("A", "<INS>") == "sv"

    def test_symbolic_dup(self) -> None:
        assert classify_variant_type("N", "<DUP:TANDEM>") == "sv"


# ---------------------------------------------------------------------------
# format_allele_key
# ---------------------------------------------------------------------------


class TestFormatAlleleKey:
    def test_basic_snp(self) -> None:
        assert format_allele_key("chr1", 100, "A", "T") == "chr1:100:A:T"

    def test_indel(self) -> None:
        assert format_allele_key("chr1", 500, "AGG", "A") == "chr1:500:AGG:A"

    def test_pos_is_1based(self) -> None:
        # VCF POS is passed through unchanged (already 1-based)
        assert format_allele_key("scaffold1", 1, "G", "A") == "scaffold1:1:G:A"

    def test_long_ref_truncated(self) -> None:
        long_ref = "A" * 30
        key = format_allele_key("chr1", 1, long_ref, "T", max_allele_len=20)
        assert key.startswith("chr1:1:" + "A" * 20 + "...")

    def test_long_alt_truncated(self) -> None:
        long_alt = "T" * 25
        key = format_allele_key("chr2", 200, "A", long_alt, max_allele_len=20)
        assert "T" * 20 + "..." in key

    def test_exact_max_length_not_truncated(self) -> None:
        ref = "A" * 20
        key = format_allele_key("chr1", 1, ref, "T", max_allele_len=20)
        assert "..." not in key


# ---------------------------------------------------------------------------
# VCF file access (requires pysam + indexed_vcf fixture)
# ---------------------------------------------------------------------------


class TestGetVcfSamples:
    def test_returns_correct_ordered_samples(self, indexed_vcf) -> None:
        from privy.io.vcf import get_vcf_samples  # noqa: PLC0415

        samples = get_vcf_samples(indexed_vcf)
        assert samples == ["T1", "T2", "O1", "O2", "O3"]

    def test_returns_list(self, indexed_vcf) -> None:
        from privy.io.vcf import get_vcf_samples  # noqa: PLC0415

        assert isinstance(get_vcf_samples(indexed_vcf), list)


class TestGetVcfContigs:
    def test_chr1_present(self, indexed_vcf) -> None:
        from privy.io.vcf import get_vcf_contigs  # noqa: PLC0415

        contigs = get_vcf_contigs(indexed_vcf)
        assert "chr1" in contigs

    def test_returns_list(self, indexed_vcf) -> None:
        from privy.io.vcf import get_vcf_contigs  # noqa: PLC0415

        assert isinstance(get_vcf_contigs(indexed_vcf), list)


class TestValidateVcfIndex:
    def test_indexed_vcf_passes(self, indexed_vcf) -> None:
        from privy.io.vcf import validate_vcf_index  # noqa: PLC0415

        # Should not raise
        validate_vcf_index(indexed_vcf)

    def test_missing_index_raises(self, tmp_path) -> None:
        from privy.io.vcf import validate_vcf_index  # noqa: PLC0415

        fake_vcf = tmp_path / "no_index.vcf.gz"
        fake_vcf.touch()
        with pytest.raises(FileNotFoundError, match="VCF index not found"):
            validate_vcf_index(fake_vcf)


class TestExtractCohortCounts:
    """Verify per-record cohort counting against known VCF positions."""

    @staticmethod
    def _fetch_record(indexed_vcf, contig: str, pos_1based: int):
        """Return the pysam.VariantRecord at a specific 1-based VCF POS."""
        import pysam  # noqa: PLC0415

        with pysam.VariantFile(str(indexed_vcf)) as vf:
            for rec in vf.fetch(contig, pos_1based - 1, pos_1based):
                return rec
        return None

    def test_strict_complete_counts(self, indexed_vcf) -> None:
        from privy.io.vcf import extract_cohort_counts  # noqa: PLC0415

        rec = self._fetch_record(indexed_vcf, "chr1", 100)
        assert rec is not None
        ts, tt, os_, ot, tm, om = extract_cohort_counts(
            rec, ["T1", "T2"], ["O1", "O2", "O3"]
        )
        assert ts == 2  # both targets carry the alt
        assert tt == 2
        assert os_ == 0  # no off-targets carry the alt
        assert ot == 3
        assert tm == 0  # no missing targets
        assert om == 0  # no missing off-targets

    def test_target_missing_counts(self, indexed_vcf) -> None:
        from privy.io.vcf import extract_cohort_counts  # noqa: PLC0415

        # pos=200: T2=./. is missing
        rec = self._fetch_record(indexed_vcf, "chr1", 200)
        assert rec is not None
        ts, tt, os_, ot, tm, om = extract_cohort_counts(
            rec, ["T1", "T2"], ["O1", "O2", "O3"]
        )
        assert ts == 1   # only T1 has the alt
        assert tt == 2
        assert os_ == 0
        assert ot == 3
        assert tm == 1   # T2 is missing
        assert om == 0

    def test_offtarget_missing_counts(self, indexed_vcf) -> None:
        from privy.io.vcf import extract_cohort_counts  # noqa: PLC0415

        # pos=300: O3=./. is missing
        rec = self._fetch_record(indexed_vcf, "chr1", 300)
        assert rec is not None
        ts, tt, os_, ot, tm, om = extract_cohort_counts(
            rec, ["T1", "T2"], ["O1", "O2", "O3"]
        )
        assert ts == 2
        assert tt == 2
        assert os_ == 0
        assert ot == 3
        assert tm == 0
        assert om == 1   # O3 is missing

    def test_both_missing_counts(self, indexed_vcf) -> None:
        from privy.io.vcf import extract_cohort_counts  # noqa: PLC0415

        # pos=400: T2 and O3 are both missing
        rec = self._fetch_record(indexed_vcf, "chr1", 400)
        assert rec is not None
        ts, tt, os_, ot, tm, om = extract_cohort_counts(
            rec, ["T1", "T2"], ["O1", "O2", "O3"]
        )
        assert ts == 1
        assert tt == 2
        assert os_ == 0
        assert ot == 3
        assert tm == 1   # T2
        assert om == 1   # O3

    def test_contradicted_counts(self, indexed_vcf) -> None:
        from privy.io.vcf import extract_cohort_counts  # noqa: PLC0415

        # pos=500: O1=0/1 — off-target carries the allele
        rec = self._fetch_record(indexed_vcf, "chr1", 500)
        assert rec is not None
        ts, tt, os_, ot, tm, om = extract_cohort_counts(
            rec, ["T1", "T2"], ["O1", "O2", "O3"]
        )
        assert ts == 2
        assert os_ == 1   # O1 carries the alt → contradiction

    def test_absent_sample_counted_as_missing(self, indexed_vcf) -> None:
        from privy.io.vcf import extract_cohort_counts  # noqa: PLC0415

        # "GHOST" is not in the VCF header — should be counted as missing
        rec = self._fetch_record(indexed_vcf, "chr1", 100)
        assert rec is not None
        ts, tt, os_, ot, tm, om = extract_cohort_counts(
            rec, ["T1", "GHOST"], ["O1", "O2", "O3"]
        )
        assert tt == 2
        assert tm == 1   # GHOST absent from VCF → missing

    def test_second_alt_allele_not_counted_for_first(self, indexed_vcf) -> None:
        from privy.io.vcf import extract_cohort_counts  # noqa: PLC0415

        # pos=800: ALT=T,G; samples have GT 0/1 (carry T, not G)
        # For alt_index=1 (G), target support should be 0
        rec = self._fetch_record(indexed_vcf, "chr1", 800)
        assert rec is not None
        ts, tt, os_, ot, tm, om = extract_cohort_counts(
            rec, ["T1", "T2"], ["O1", "O2", "O3"], alt_index=1
        )
        assert ts == 0   # no sample has allele value 2 (G)
        assert os_ == 0
