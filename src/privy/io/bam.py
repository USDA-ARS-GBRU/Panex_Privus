"""BAM depth and allele-level queries for Panex Privus.

BAM is a support layer — it provides read-level evidence at known candidate
loci, not a de-novo variant caller.  Every function performs a focused,
position-specific query rather than scanning the entire file.

Key design choices:
  - pysam is imported lazily inside each function to avoid import-time cost
    when BAM support is disabled.
  - :func:`query_allele_counts_at_locus` returns ``(total_depth, 0, 0)``
    for indels to signal depth-only mode without raising an exception.
  - Unmapped, secondary, and supplementary reads are always excluded.
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path
from typing import Optional

log = logging.getLogger("privy.io.bam")


def validate_bam_index(bam_path: Path) -> None:
    """Raise :exc:`FileNotFoundError` if no BAM index (.bai or .csi) is found.

    Args:
        bam_path: Path to the BAM file.

    Raises:
        FileNotFoundError: If neither ``<bam_path>.bai`` nor
            ``<bam_path>.csi`` exists alongside the BAM file.
    """
    bai = Path(str(bam_path) + ".bai")
    csi = Path(str(bam_path) + ".csi")
    if not bai.exists() and not csi.exists():
        raise FileNotFoundError(
            f"BAM index not found for {bam_path}. "
            "Run 'samtools index' to create a .bai or .csi index."
        )


def get_bam_sample_name(bam_path: Path) -> Optional[str]:
    """Return the SM tag from the first @RG read-group header, or *None*.

    Args:
        bam_path: Path to the BAM file (does not need to be indexed).

    Returns:
        Sample name string, or ``None`` if no @RG header or SM tag is present.
    """
    import pysam  # noqa: PLC0415

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        header_dict = bam.header.to_dict()
        rg_list = header_dict.get("RG", [])
        if rg_list:
            return rg_list[0].get("SM")
    return None


def load_bam_manifest(manifest_path: Path) -> list[dict[str, str]]:
    """Load a TSV BAM manifest mapping file paths to sample identifiers.

    The manifest must have a header row.  Lines beginning with ``#`` are
    treated as comments and skipped.  Required columns: ``bam_path``,
    ``sample_id``.  An optional ``group`` column is preserved but not
    interpreted by this function.

    Args:
        manifest_path: Path to the TSV manifest file.

    Returns:
        List of row dicts, one per non-comment data row.  Each dict contains
        at least ``"bam_path"`` and ``"sample_id"`` keys.

    Raises:
        FileNotFoundError: If *manifest_path* does not exist.
        ValueError: If the required columns are absent.
    """
    if not manifest_path.exists():
        raise FileNotFoundError(f"BAM manifest not found: {manifest_path}")

    rows: list[dict[str, str]] = []
    with open(manifest_path, newline="", encoding="utf-8") as fh:
        non_comment = (line for line in fh if not line.startswith("#"))
        reader = csv.DictReader(non_comment, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("BAM manifest is missing a header row.")
        fieldnames = set(reader.fieldnames)
        if "bam_path" not in fieldnames:
            raise ValueError("BAM manifest must contain a 'bam_path' column.")
        if "sample_id" not in fieldnames:
            raise ValueError("BAM manifest must contain a 'sample_id' column.")
        for row in reader:
            rows.append(dict(row))

    log.debug("Loaded %d entries from BAM manifest %s", len(rows), manifest_path)
    return rows


def query_position_depth(
    bam_path: Path,
    contig: str,
    start: int,
    end: int,
    min_mapq: int = 20,
) -> list[int]:
    """Return per-position read depth over [start, end) in a BAM file.

    Reads with mapping quality below *min_mapq* are excluded.
    Unmapped, secondary, and supplementary reads are always excluded.

    Args:
        bam_path: Path to an indexed BAM file.
        contig: Contig/chromosome name.
        start: 0-based start position (inclusive).
        end: 0-based end position (exclusive).
        min_mapq: Minimum mapping quality for read inclusion.

    Returns:
        List of integer depth values of length ``end - start``.  Returns all
        zeros if the contig is absent from the BAM or the region has no reads.
    """
    import pysam  # noqa: PLC0415

    width = max(end - start, 0)
    if width == 0:
        return []

    def _mapq_filter(read: "pysam.AlignedSegment") -> bool:  # type: ignore[name-defined]
        return (
            not read.is_unmapped
            and not read.is_secondary
            and not read.is_supplementary
            and read.mapping_quality >= min_mapq
        )

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        try:
            a, c, g, t = bam.count_coverage(
                contig, start, end,
                quality_threshold=0,
                read_callback=_mapq_filter,
            )
        except (ValueError, KeyError):
            return [0] * width
        return [a[i] + c[i] + g[i] + t[i] for i in range(len(a))]


def query_allele_counts_at_locus(
    bam_path: Path,
    contig: str,
    pos: int,
    ref_allele: str,
    alt_allele: str,
    min_mapq: int = 20,
    min_baseq: int = 20,
) -> tuple[int, int, int]:
    """Return ``(ref_count, alt_count, other_count)`` at a SNP position.

    For non-SNP variants (ref or alt allele length > 1) this function
    returns ``(total_depth, 0, 0)`` to signal depth-only mode, because
    indel/SV alleles cannot be distinguished from depth alone.

    Args:
        bam_path: Path to an indexed BAM file.
        contig: Contig/chromosome name.
        pos: 0-based reference position of the SNP.
        ref_allele: Reference allele string.
        alt_allele: Alternate allele string.
        min_mapq: Minimum mapping quality for read inclusion.
        min_baseq: Minimum base quality at the query position.

    Returns:
        ``(ref_count, alt_count, other_count)`` for SNPs, or
        ``(total_depth, 0, 0)`` for indels.
    """
    import pysam  # noqa: PLC0415

    is_snp = len(ref_allele) == 1 and len(alt_allele) == 1
    if not is_snp:
        depths = query_position_depth(bam_path, contig, pos, pos + 1, min_mapq)
        return (depths[0] if depths else 0, 0, 0)

    ref_base = ref_allele.upper()
    alt_base = alt_allele.upper()
    ref_count = 0
    alt_count = 0
    other_count = 0

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        try:
            for pileup_col in bam.pileup(
                contig, pos, pos + 1,
                min_mapping_quality=min_mapq,
                min_base_quality=min_baseq,
                truncate=True,
            ):
                if pileup_col.reference_pos != pos:
                    continue
                for pileup_read in pileup_col.pileups:
                    if pileup_read.is_del or pileup_read.is_refskip:
                        other_count += 1
                        continue
                    qpos = pileup_read.query_position
                    if qpos is None:
                        continue
                    seq = pileup_read.alignment.query_sequence
                    if seq is None:
                        continue
                    base = seq[qpos].upper()
                    if base == ref_base:
                        ref_count += 1
                    elif base == alt_base:
                        alt_count += 1
                    else:
                        other_count += 1
        except (ValueError, KeyError):
            pass

    return (ref_count, alt_count, other_count)
