"""Shared pytest fixtures for Panex Privus tests.

Provides:
    - ``indexed_vcf``: a bgzip-compressed, tabix-indexed synthetic VCF in a
      temporary directory.  Built with ``pysam.tabix_compress`` / ``tabix_index``
      so that bgzip/tabix do not need to be on PATH.
    - ``small_cohort``: a :class:`~privy.core.cohort.CohortDefinition` matching
      the samples in the synthetic VCF.
    - ``default_cfg``: a default :class:`~privy.core.config.PrivyConfig`.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from privy.core.cohort import CohortDefinition
from privy.core.config import default_config

# ---------------------------------------------------------------------------
# Synthetic VCF content
#
# Samples: T1, T2 (targets) and O1, O2, O3 (off-targets).
# Records are designed to exercise every strictness class, filter path, and
# variant type that the scan backend must handle.
# ---------------------------------------------------------------------------

_VCF_LINES: list[str] = [
    "##fileformat=VCFv4.2",
    '##FILTER=<ID=PASS,Description="All filters passed">',
    '##FILTER=<ID=FAIL,Description="Filter failed">',
    "##contig=<ID=chr1,length=10000>",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tT1\tT2\tO1\tO2\tO3",
    # pos=100  strict_complete: all targets carry alt, no off-targets, no missing
    "chr1\t100\t.\tA\tT\t50\tPASS\t.\tGT\t0/1\t0/1\t0/0\t0/0\t0/0",
    # pos=200  strict_target_missing: T2 is missing (./.)
    "chr1\t200\t.\tA\tT\t50\tPASS\t.\tGT\t0/1\t./.\t0/0\t0/0\t0/0",
    # pos=300  strict_offtarget_missing: O3 is missing
    "chr1\t300\t.\tA\tT\t50\tPASS\t.\tGT\t0/1\t0/1\t0/0\t0/0\t./.",
    # pos=400  strict_both_missing: T2 and O3 are both missing
    "chr1\t400\t.\tA\tT\t50\tPASS\t.\tGT\t0/1\t./.\t0/0\t0/0\t./.",
    # pos=500  contradicted: O1 also carries the allele → not emitted
    "chr1\t500\t.\tA\tT\t50\tPASS\t.\tGT\t0/1\t0/1\t0/1\t0/0\t0/0",
    # pos=600  FILTER=FAIL: skipped by pass_only (default)
    "chr1\t600\t.\tA\tT\t50\tFAIL\t.\tGT\t0/1\t0/1\t0/0\t0/0\t0/0",
    # pos=700  low QUAL (5): skipped only when min_qual > 5
    "chr1\t700\t.\tA\tT\t5\tPASS\t.\tGT\t0/1\t0/1\t0/0\t0/0\t0/0",
    # pos=800  multiallelic (ALT=T,G): skipped entirely when allow_multiallelic=False;
    #          ALT=G (alt_index=1) fails target-support threshold so only T is emitted
    "chr1\t800\t.\tA\tT,G\t50\tPASS\t.\tGT\t0/1\t0/1\t0/0\t0/0\t0/0",
    # pos=900  indel (AGG→A): strict_complete
    "chr1\t900\t.\tAGG\tA\t50\tPASS\t.\tGT\t0/1\t0/1\t0/0\t0/0\t0/0",
]

VCF_TEXT: str = "\n".join(_VCF_LINES) + "\n"


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="session")
def vcf_text() -> str:
    """Raw VCF text for the small synthetic cohort (session-scoped constant)."""
    return VCF_TEXT


@pytest.fixture
def indexed_vcf(tmp_path: Path, vcf_text: str) -> Path:
    """Write, bgzip-compress, and tabix-index a synthetic VCF in *tmp_path*.

    Uses ``pysam.tabix_compress`` and ``pysam.tabix_index`` so the test
    fixtures are self-contained and do not require bgzip/tabix on PATH.

    Returns:
        Path to the ``cohort.vcf.gz`` file.  The ``.tbi`` index is created
        alongside it automatically.
    """
    import pysam  # noqa: PLC0415

    plain = tmp_path / "cohort.vcf"
    plain.write_text(vcf_text)
    gz_path = str(tmp_path / "cohort.vcf.gz")
    pysam.tabix_compress(str(plain), gz_path, force=True)
    pysam.tabix_index(gz_path, preset="vcf", force=True)
    return Path(gz_path)


@pytest.fixture
def small_cohort() -> CohortDefinition:
    """CohortDefinition matching the T1/T2 targets and O1/O2/O3 off-targets in the VCF."""
    return CohortDefinition.from_lists(
        targets=["T1", "T2"],
        off_targets=["O1", "O2", "O3"],
    )


@pytest.fixture
def default_cfg():
    """Default PrivyConfig with no overrides applied."""
    return default_config()
