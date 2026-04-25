"""Shared pytest fixtures for Panex Privus tests.

Provides:
    - ``indexed_vcf``: a bgzip-compressed, tabix-indexed synthetic VCF in a
      temporary directory.  Built with ``pysam.tabix_compress`` / ``tabix_index``
      so that bgzip/tabix do not need to be on PATH.
    - ``small_cohort``: a :class:`~privy.core.cohort.CohortDefinition` matching
      the samples in the synthetic VCF.
    - ``default_cfg``: a default :class:`~privy.core.config.PrivyConfig`.
    - BAM fixtures: ``bam_target_t1``, ``bam_offtarget_o1``,
      ``bam_offtarget_o1_with_alt``, ``bam_low_depth_t2`` — sorted, indexed
      synthetic BAMs built with pysam for BAM support layer tests.
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


# ---------------------------------------------------------------------------
# BAM fixtures
#
# The synthetic VCF has SNP hits at (1-based) positions 100, 200, 300, 400,
# 700, 800 and an indel at 900.  In 0-based coords the primary test SNP is
# at position 99 (VCF pos 100, REF=A, ALT=T).
#
# Each BAM is created pre-sorted (reads emitted in coordinate order) so
# pysam.sort() is not required.  Only pysam.index() is called, which uses
# htslib directly.
# ---------------------------------------------------------------------------

_BAM_CONTIGS = [("chr1", 10000)]
# All test reads: length 30, starting at pos 85 (0-based).
# Offset 14 within the read lands on position 99 (0-based) = VCF pos 100.
_READ_START = 85
_READ_LEN = 30
_ALT_OFFSET = 14   # within the read → reference position 99


def _write_synthetic_bam(
    path: Path,
    sample: str,
    contigs: list[tuple[str, int]],
    reads: list[dict],
) -> Path:
    """Write a coordinate-sorted, indexed BAM with synthetic reads.

    ``reads`` must already be sorted by ``start`` or be at the same position.
    Each entry is a dict with keys: ``seq`` (str), ``start`` (int 0-based),
    and optionally ``mapq`` (int, default 60).

    Returns the path to the written BAM.
    """
    import pysam  # noqa: PLC0415

    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": name, "LN": length} for name, length in contigs],
        "RG": [{"ID": sample, "SM": sample}],
    })

    sorted_reads = sorted(reads, key=lambda r: r.get("start", 0))

    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for i, rd in enumerate(sorted_reads):
            seg = pysam.AlignedSegment(header)
            seg.query_name = f"read_{i:04d}"
            seg.query_sequence = rd["seq"]
            seg.flag = 0
            seg.reference_id = 0
            seg.reference_start = rd["start"]
            seg.mapping_quality = rd.get("mapq", 60)
            seg.cigar = [(0, len(rd["seq"]))]
            seg.query_qualities = pysam.qualitystring_to_array("I" * len(rd["seq"]))
            seg.set_tag("RG", sample)
            bam.write(seg)

    pysam.index(str(path))
    return path


@pytest.fixture
def bam_target_t1(tmp_path: Path) -> Path:
    """T1 BAM: 12 reads with ALT allele (T) at 0-based position 99 (VCF pos 100).

    depth=12, alt_count=12, allele_fraction=1.0 → should classify as SUPPORT.
    """
    seq_with_alt = "A" * _ALT_OFFSET + "T" + "A" * (_READ_LEN - _ALT_OFFSET - 1)
    reads = [{"seq": seq_with_alt, "start": _READ_START} for _ in range(12)]
    return _write_synthetic_bam(tmp_path / "T1.bam", "T1", _BAM_CONTIGS, reads)


@pytest.fixture
def bam_offtarget_o1(tmp_path: Path) -> Path:
    """O1 BAM: 12 reads with only REF allele (A) at position 99.

    depth=12, alt_count=0, allele_fraction=0.0 → should classify as ABSENCE.
    """
    seq_ref_only = "A" * _READ_LEN
    reads = [{"seq": seq_ref_only, "start": _READ_START} for _ in range(12)]
    return _write_synthetic_bam(tmp_path / "O1.bam", "O1", _BAM_CONTIGS, reads)


@pytest.fixture
def bam_offtarget_o1_with_alt(tmp_path: Path) -> Path:
    """O1 BAM: 12 reads carrying ALT allele (T) at position 99.

    Simulates an off-target sample that carries the private allele
    → should classify as CONTRADICTION.
    """
    seq_with_alt = "A" * _ALT_OFFSET + "T" + "A" * (_READ_LEN - _ALT_OFFSET - 1)
    reads = [{"seq": seq_with_alt, "start": _READ_START} for _ in range(12)]
    return _write_synthetic_bam(
        tmp_path / "O1_alt.bam", "O1", _BAM_CONTIGS, reads
    )


@pytest.fixture
def bam_low_depth_t2(tmp_path: Path) -> Path:
    """T2 BAM: only 3 reads at position 99 (below default min_depth=8).

    Should classify as UNINFORMATIVE regardless of allele content.
    """
    seq_with_alt = "A" * _ALT_OFFSET + "T" + "A" * (_READ_LEN - _ALT_OFFSET - 1)
    reads = [{"seq": seq_with_alt, "start": _READ_START} for _ in range(3)]
    return _write_synthetic_bam(tmp_path / "T2.bam", "T2", _BAM_CONTIGS, reads)
