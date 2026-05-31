"""Allelic diversity estimators for multi-allelic microhaplotype loci.

Pure-Python implementations of the standard multi-allelic diversity statistics
(Nei gene diversity / expected heterozygosity, effective number of alleles,
combined alternative allele frequency).  Operates over the per-genome allele
calls carried by :class:`~privy.microhap.model.Microhaplotype`.

Citations: Nei (1973) PNAS 70:3321 for gene diversity.  See
scratch/notes/70_popgen_estimators_math.md.  Validate against
hierfstat / adegenet on shared inputs.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

from privy.microhap.model import Microhaplotype
from privy.polyploid.dosage import observed_heterozygosity


def allele_frequencies(
    alleles: Mapping[str, str],
    genomes: Sequence[str] | None = None,
) -> dict[str, float]:
    """Allele-id → frequency over the (optionally cohort-restricted) genomes with data.

    Args:
        alleles: genome/path id → allele id.
        genomes: restrict to these genomes (default: all in *alleles*).
    """
    if genomes is None:
        calls = list(alleles.values())
    else:
        keep = set(genomes)
        calls = [a for g, a in alleles.items() if g in keep]
    n = len(calls)
    if n == 0:
        return {}
    counts: dict[str, int] = {}
    for a in calls:
        counts[a] = counts.get(a, 0) + 1
    return {a: c / n for a, c in counts.items()}


def nei_gene_diversity(freqs: Mapping[str, float]) -> float:
    """Nei gene diversity / expected heterozygosity ``He = 1 - Σ p_u²``."""
    return 1.0 - sum(p * p for p in freqs.values())


def effective_n_alleles(freqs: Mapping[str, float]) -> float:
    """Effective number of alleles ``Ne = 1 / Σ p_u²`` (= 1/(1-He))."""
    homozygosity = sum(p * p for p in freqs.values())
    if homozygosity <= 0:
        return 0.0
    return 1.0 / homozygosity


@dataclass(frozen=True)
class LocusDiversity:
    """Diversity summary for one microhaplotype locus over a genome set."""

    locus_id: str
    n_genomes: int
    n_alleles: int
    gene_diversity: float     # Nei He
    effective_alleles: float
    aaf: float                # combined alternative allele frequency


def locus_diversity(
    mh: Microhaplotype,
    genomes: Sequence[str] | None = None,
) -> LocusDiversity:
    """Compute diversity metrics for *mh*, optionally restricted to *genomes*."""
    freqs = allele_frequencies(mh.alleles, genomes)
    keep = None if genomes is None else set(genomes)
    n_genomes = sum(1 for g in mh.alleles if keep is None or g in keep)
    he = nei_gene_diversity(freqs)
    if mh.ref_allele is not None and mh.ref_allele in freqs:
        aaf = 1.0 - freqs[mh.ref_allele]
    elif freqs:
        aaf = 1.0 - max(freqs.values())
    else:
        aaf = 0.0
    return LocusDiversity(
        locus_id=mh.locus_id,
        n_genomes=n_genomes,
        n_alleles=len(freqs),
        gene_diversity=he,
        effective_alleles=effective_n_alleles(freqs),
        aaf=aaf,
    )


def inbreeding_fis(
    mh: Microhaplotype,
    samples_to_paths: Mapping[str, Sequence[str]],
) -> float | None:
    """Inbreeding coefficient F_IS = 1 − H_o / H_e at a locus.

    H_o is the observed heterozygosity across per-sample genotypes (a sample is
    heterozygous when its haplotype paths carry >1 distinct allele); H_e is the
    expected heterozygosity (Nei gene diversity) from allele frequencies over all
    gene copies.  Returns None when H_e is 0 (monomorphic — F_IS undefined).
    """
    he = nei_gene_diversity(allele_frequencies(mh.alleles))
    if he <= 0:
        return None
    ho = observed_heterozygosity(mh, samples_to_paths)
    return 1.0 - ho / he
