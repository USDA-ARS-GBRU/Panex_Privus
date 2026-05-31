"""Per-sample allele dosage from haplotype paths (polyploid genotypes).

Groups a sample's PanSN haplotype paths (``sample#hap#contig``) into one genotype
and counts allele dosage at a microhaplotype locus — yielding 0..ploidy alt-allele
dosage, observed heterozygosity, and the dosage encoding consumed by genomic-
prediction tools (rrBLUP / BGLR / VanRaden GRM).

Citations / definitions: scratch/notes/70_popgen_estimators_math.md.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

from privy.microhap.model import Microhaplotype
from privy.synteny.model import split_pansn


def group_paths_by_sample(
    paths: Sequence[str],
    delimiter: str = "#",
) -> dict[str, list[str]]:
    """Group PanSN path ids by sample, preserving first-seen order.

    ``sample0#0#chr1`` and ``sample0#1#chr1`` group under ``sample0``.
    """
    groups: dict[str, list[str]] = {}
    for path in paths:
        sample = split_pansn(path, delimiter)[0]
        groups.setdefault(sample, []).append(path)
    return groups


def observed_ploidy(sample_paths: Sequence[str]) -> int:
    """Number of haplotype paths a sample contributes (its observed ploidy)."""
    return len(sample_paths)


def sample_allele_dosage(
    mh: Microhaplotype,
    sample_paths: Sequence[str],
) -> dict[str, int]:
    """Count of each allele id across a sample's haplotype paths at locus *mh*.

    Paths with no call at this locus are skipped (so the dosages sum to the number
    of *called* haplotypes, which may be < ploidy under missing data).
    """
    counts: dict[str, int] = {}
    for path in sample_paths:
        allele = mh.alleles.get(path)
        if allele is None:
            continue
        counts[allele] = counts.get(allele, 0) + 1
    return counts


def alt_dosage(
    mh: Microhaplotype,
    sample_paths: Sequence[str],
    ref_allele: str | None = None,
) -> int | None:
    """Alternate-allele dosage (count of non-reference alleles) for a sample.

    Uses *ref_allele* if given, else the locus's reference allele.  Returns the
    {0..ploidy} dosage GP tools expect, or None when the sample has no call here.
    """
    ref = ref_allele if ref_allele is not None else mh.ref_allele
    counts = sample_allele_dosage(mh, sample_paths)
    if not counts:
        return None
    total = sum(counts.values())
    if ref is None:
        return total   # no reference allele defined -> all copies count as alternate
    return total - counts.get(ref, 0)


def is_heterozygous(mh: Microhaplotype, sample_paths: Sequence[str]) -> bool:
    """True when a sample carries more than one distinct allele at *mh*."""
    alleles = {mh.alleles[p] for p in sample_paths if p in mh.alleles}
    return len(alleles) > 1


def observed_heterozygosity(
    mh: Microhaplotype,
    samples_to_paths: Mapping[str, Sequence[str]],
) -> float:
    """Fraction of samples (with ≥1 called haplotype) that are heterozygous at *mh*."""
    total = 0
    het = 0
    for sample_paths in samples_to_paths.values():
        called = [p for p in sample_paths if p in mh.alleles]
        if not called:
            continue
        total += 1
        if len({mh.alleles[p] for p in called}) > 1:
            het += 1
    return het / total if total else 0.0
