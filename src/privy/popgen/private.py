"""Private-allele metrics: counts and rarefied private allelic richness.

Ties Privy's core target-private question to standard population-genetics: how
many alleles are private to each cohort, and (rarefied for fair comparison across
unequal sample sizes) the expected private allelic richness.

Citations: Kalinowski (2004) Conserv Genet 5:539; Szpiech et al. (2008) ADZE,
Bioinformatics 24:2498.  See scratch/notes/70_popgen_estimators_math.md.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from math import comb

from privy.microhap.model import Microhaplotype


def _cohort_allele_counts(mh: Microhaplotype, paths: Sequence[str]) -> dict[str, int]:
    """Count each allele's gene copies (one per called haplotype path) in a cohort."""
    counts: dict[str, int] = {}
    for path in paths:
        allele = mh.alleles.get(path)
        if allele is None:
            continue
        counts[allele] = counts.get(allele, 0) + 1
    return counts


def private_allele_counts(
    loci: Sequence[Microhaplotype],
    cohorts: Mapping[str, Sequence[str]],
) -> dict[str, int]:
    """Count alleles private to each cohort (present in it, absent from all others)."""
    result = {cohort: 0 for cohort in cohorts}
    for mh in loci:
        counts = {c: _cohort_allele_counts(mh, paths) for c, paths in cohorts.items()}
        for cohort in cohorts:
            others = [c for c in cohorts if c != cohort]
            for allele in counts[cohort]:
                if all(allele not in counts[o] for o in others):
                    result[cohort] += 1
    return result


def rarefied_private_allelic_richness(
    loci: Sequence[Microhaplotype],
    cohorts: Mapping[str, Sequence[str]],
    g: int | None = None,
) -> dict[str, float]:
    """Expected private allelic richness per cohort, rarefied to *g* gene copies.

    For each locus and allele, sums the probability the allele *is* drawn in the
    focal cohort's size-*g* subsample times the probability it is *absent* from
    every other cohort's size-*g* subsample (Kalinowski 2004 / ADZE).  Loci where a
    cohort has fewer than *g* called gene copies are skipped for that cohort.

    Args:
        g: rarefaction subsample size in gene copies; defaults to the smallest
            cohort gene-copy count observed at each locus (computed per locus).
    """
    result = {cohort: 0.0 for cohort in cohorts}
    for mh in loci:
        counts = {c: _cohort_allele_counts(mh, paths) for c, paths in cohorts.items()}
        totals = {c: sum(v.values()) for c, v in counts.items()}
        locus_g = g if g is not None else min((n for n in totals.values() if n > 0), default=0)
        if locus_g <= 0:
            continue
        alleles = {a for c in counts.values() for a in c}
        for focal in cohorts:
            n_focal = totals[focal]
            if n_focal < locus_g:
                continue
            others = [c for c in cohorts if c != focal and totals[c] >= locus_g]
            for allele in alleles:
                p_present = _prob_present(n_focal, counts[focal].get(allele, 0), locus_g)
                if p_present <= 0:
                    continue
                p_absent_elsewhere = 1.0
                for other in others:
                    p_absent_elsewhere *= _prob_absent(
                        totals[other], counts[other].get(allele, 0), locus_g
                    )
                result[focal] += p_present * p_absent_elsewhere
    return result


def _prob_absent(n: int, n_allele: int, g: int) -> float:
    """Probability allele (count *n_allele* of *n*) is absent from a size-*g* draw."""
    if g > n:
        return 0.0
    denom = comb(n, g)
    if denom == 0:
        return 1.0
    return comb(n - n_allele, g) / denom


def _prob_present(n: int, n_allele: int, g: int) -> float:
    """Probability allele is observed at least once in a size-*g* draw."""
    return 1.0 - _prob_absent(n, n_allele, g)
