"""VCF/BCF streaming reader for Panex Privus.

Wraps ``pysam.VariantFile`` to provide indexed, chunked iteration over
multisample VCF records.  All coordinate output uses 0-based half-open
convention.  VCF 1-based POS values are converted on read.

Design notes:
    - ``pysam`` is imported lazily inside each function so the rest of the
      package remains importable in environments without pysam (e.g., for
      unit-testing pure-Python logic).
    - ``stream_vcf_records`` uses a ``try/finally`` pattern rather than a
      ``with`` statement to avoid the well-known
      ``generator + context manager = premature close`` antipattern.
"""

from __future__ import annotations

import logging
from collections.abc import Iterable, Iterator
from pathlib import Path
from typing import Protocol, cast

log = logging.getLogger("privy.io.vcf")

Genotype = tuple[int | None, ...]


class VariantSampleCall(Protocol):
    """Subset of the pysam sample-call interface used by this project."""

    def __getitem__(self, key: str) -> Genotype | None: ...


class VariantSamples(Protocol):
    """Mapping-like sample container exposed by pysam records."""

    def __getitem__(self, sample: str) -> VariantSampleCall: ...


class VariantRecordLike(Protocol):
    """Structural type for the pysam record fields consumed by the scan path."""

    @property
    def filter(self) -> Iterable[str]: ...

    @property
    def qual(self) -> float | None: ...

    @property
    def alts(self) -> tuple[str, ...] | None: ...

    @property
    def ref(self) -> str: ...

    @property
    def pos(self) -> int: ...

    @property
    def chrom(self) -> str: ...

    @property
    def samples(self) -> VariantSamples: ...


# ---------------------------------------------------------------------------
# Index validation
# ---------------------------------------------------------------------------

def validate_vcf_index(vcf_path: Path) -> None:
    """Raise ``FileNotFoundError`` if no tabix or CSI index is found.

    Checks for both ``<file>.tbi`` and ``<file>.csi`` alongside the VCF.

    Args:
        vcf_path: Path to the ``.vcf.gz`` (or ``.bcf``) file.

    Raises:
        FileNotFoundError: If neither index variant is present.
    """
    tbi = Path(str(vcf_path) + ".tbi")
    csi = Path(str(vcf_path) + ".csi")
    if not tbi.exists() and not csi.exists():
        raise FileNotFoundError(
            f"VCF index not found for: {vcf_path}\n"
            "Create one with:\n"
            f"  bgzip {vcf_path.with_suffix('')}\n"
            f"  tabix -p vcf {vcf_path}"
        )


# ---------------------------------------------------------------------------
# Header queries
# ---------------------------------------------------------------------------

def get_vcf_samples(vcf_path: Path) -> list[str]:
    """Return the sample names declared in the VCF header.

    Args:
        vcf_path: Path to the VCF/BCF file.

    Returns:
        Ordered list of sample names.
    """
    import pysam  # noqa: PLC0415

    with pysam.VariantFile(str(vcf_path)) as vf:
        return list(vf.header.samples)


def get_vcf_contigs(vcf_path: Path) -> list[str]:
    """Return the contig names declared in the VCF header.

    Args:
        vcf_path: Path to the VCF/BCF file.

    Returns:
        List of contig names in header order.
    """
    import pysam  # noqa: PLC0415

    with pysam.VariantFile(str(vcf_path)) as vf:
        return [str(contig) for contig in vf.header.contigs]


# ---------------------------------------------------------------------------
# Record streaming
# ---------------------------------------------------------------------------

def stream_vcf_records(
    vcf_path: Path,
    contig: str | None = None,
    start: int | None = None,
    end: int | None = None,
) -> Iterator[VariantRecordLike]:
    """Stream VCF records from an (optionally indexed) VCF/BCF file.

    Uses ``pysam.VariantFile.fetch()`` for indexed access when a contig is
    specified; otherwise iterates all records sequentially.

    Coordinates follow the pysam/BED convention: ``start`` and ``end``
    are **0-based, half-open**.  pysam's ``fetch(contig, start, end)`` also
    uses 0-based half-open, so no conversion is needed.

    Args:
        vcf_path: Path to the VCF/BCF file.
        contig: Optional contig name for region-restricted fetch.
        start: Optional 0-based start position (requires *contig*).
        end: Optional 0-based end position (requires *contig*).

    Yields:
        ``pysam.VariantRecord`` objects.

    Note:
        The function manages the file handle with ``try/finally`` rather than
        a ``with`` block to avoid the generator-context-manager interaction
        that would close the handle before the caller has consumed all records.
    """
    import pysam  # noqa: PLC0415

    vf = pysam.VariantFile(str(vcf_path))
    try:
        if contig is not None:
            if start is not None and end is not None:
                iterator = vf.fetch(contig, start, end)
            elif start is not None:
                iterator = vf.fetch(contig, start)
            else:
                iterator = vf.fetch(contig)
        else:
            iterator = vf.fetch()

        for record in iterator:
            yield cast(VariantRecordLike, record)
    finally:
        vf.close()


# ---------------------------------------------------------------------------
# Genotype helpers
# ---------------------------------------------------------------------------

def is_missing_genotype(gt: Genotype | None) -> bool:
    """Return True if all alleles in *gt* are missing.

    Args:
        gt: Genotype tuple from ``pysam`` (values are ``int | None``).
            ``None`` values represent the ``.`` allele in ``./. `` calls.

    Returns:
        ``True`` if *gt* is ``None`` or any allele is ``None``.
    """
    return gt is None or any(a is None for a in gt)


def has_alt_allele(gt: Genotype, alt_index: int) -> bool:
    """Return True if *gt* contains at least one copy of the specified alt allele.

    In pysam GT tuples, allele values are 0-based integer indices:
    ``0`` = REF, ``1`` = first ALT, ``2`` = second ALT, ``None`` = missing.

    Args:
        gt: Genotype tuple (e.g., ``(0, 1)`` for heterozygous first ALT).
        alt_index: 0-based alternate allele index (0 = first ALT in VCF).

    Returns:
        ``True`` if the GT value ``alt_index + 1`` appears in *gt*.
    """
    return (alt_index + 1) in gt


def classify_variant_type(ref: str, alt: str) -> str:
    """Return a variant type string for a REF/ALT pair.

    Returns:
        ``"sv"`` for symbolic alleles (``<DEL>``, ``<INS>``, etc.),
        ``"snp"`` for single-nucleotide changes,
        ``"indel"`` otherwise.
    """
    if alt.startswith("<"):
        return "sv"
    if len(ref) == 1 and len(alt) == 1:
        return "snp"
    return "indel"


def format_allele_key(
    contig: str,
    pos: int,
    ref: str,
    alt: str,
    max_allele_len: int = 20,
) -> str:
    """Format a canonical allele key string.

    Format: ``contig:pos:ref:alt`` where *pos* is the **1-based VCF POS**
    (not converted to 0-based) so the key is human-readable and matches
    standard VCF notation.

    Long alleles are truncated with ``...`` for readability.

    Args:
        contig: Contig name.
        pos: 1-based VCF POS.
        ref: Reference allele string.
        alt: Alternate allele string.
        max_allele_len: Maximum allele length before truncation.
    """
    ref_str = ref[:max_allele_len] + ("..." if len(ref) > max_allele_len else "")
    alt_str = alt[:max_allele_len] + ("..." if len(alt) > max_allele_len else "")
    return f"{contig}:{pos}:{ref_str}:{alt_str}"


# ---------------------------------------------------------------------------
# Cohort counting — the core per-record function
# ---------------------------------------------------------------------------

def extract_cohort_counts(
    record: VariantRecordLike,
    target_samples: list[str],
    offtarget_samples: list[str],
    alt_index: int = 0,
) -> tuple[int, int, int, int, int, int]:
    """Count target and off-target support, absence, and missingness for one alt allele.

    This is the function that feeds directly into
    :func:`~privy.core.patterns.classify_strictness`.

    Args:
        record: A ``pysam.VariantRecord`` from :func:`stream_vcf_records`.
        target_samples: List of target sample names.
        offtarget_samples: List of off-target sample names.
        alt_index: 0-based index of the alternate allele to evaluate
            (0 = first ALT in the VCF record).

    Returns:
        Tuple of::

            (target_support_n, target_total_n,
             offtarget_support_n, offtarget_total_n,
             target_missing_n, offtarget_missing_n)

        where ``target_total_n`` = ``len(target_samples)`` and similarly for
        off-targets.  Samples absent from the VCF are counted as missing.
    """
    target_support = 0
    target_missing = 0
    offtarget_support = 0
    offtarget_missing = 0

    samples = record.samples

    for sample in target_samples:
        try:
            gt = samples[sample]["GT"]
        except (KeyError, TypeError):
            # pysam raises KeyError for sample names absent from the VCF header
            target_missing += 1
            continue
        if is_missing_genotype(gt):
            target_missing += 1
        elif gt is not None and has_alt_allele(gt, alt_index):
            target_support += 1
        # else: absent (0/0 or different alt) — counted implicitly

    for sample in offtarget_samples:
        try:
            gt = samples[sample]["GT"]
        except (KeyError, TypeError):
            offtarget_missing += 1
            continue
        if is_missing_genotype(gt):
            offtarget_missing += 1
        elif gt is not None and has_alt_allele(gt, alt_index):
            offtarget_support += 1

    return (
        target_support,
        len(target_samples),
        offtarget_support,
        len(offtarget_samples),
        target_missing,
        offtarget_missing,
    )
