"""Cohort differentiation for multi-allelic microhaplotype loci.

Quantifies how strongly each locus separates the target cohort from the
off-target cohort — the breeder-actionable signal: a locus with G_ST/F_ST ≈ 1 (no
shared allele) is a fully diagnostic / private marker for the target group.

Multi-allelic-native estimators (Nei's G_ST and Jost's D); genome-wide values use
the ratio-of-averages convention.  Citations: Nei (1973); Jost (2008) Mol Ecol
17:4015; Meirmans & Hedrick (2011).  See scratch/notes/70_popgen_estimators_math.md.
Validate against mmod / hierfstat.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

from privy.microhap.model import Microhaplotype
from privy.popgen.diversity import allele_frequencies, nei_gene_diversity


@dataclass(frozen=True)
class LocusDifferentiation:
    """Differentiation of one locus between the target and off-target cohorts."""

    locus_id: str
    target_n: int
    offtarget_n: int
    h_s: float            # mean within-cohort gene diversity
    h_t: float            # total gene diversity (mean allele freqs)
    gst: float            # Nei G_ST = (H_T - H_S)/H_T
    jost_d: float         # Jost's D
    is_diagnostic: bool   # cohorts share no allele (fully private split)


def _present_freqs(mh: Microhaplotype, genomes: Sequence[str]) -> dict[str, float]:
    return allele_frequencies(mh.alleles, genomes)


def locus_differentiation(
    mh: Microhaplotype,
    targets: Sequence[str],
    off_targets: Sequence[str],
) -> LocusDifferentiation | None:
    """Differentiation of *mh* between target and off-target cohorts.

    Returns None when either cohort has no called genome at this locus.
    """
    f_t = _present_freqs(mh, targets)
    f_o = _present_freqs(mh, off_targets)
    n_t = sum(1 for g in targets if g in mh.alleles)
    n_o = sum(1 for g in off_targets if g in mh.alleles)
    if not f_t or not f_o:
        return None

    he_t = nei_gene_diversity(f_t)
    he_o = nei_gene_diversity(f_o)
    h_s = (he_t + he_o) / 2.0

    alleles = set(f_t) | set(f_o)
    mean_freq = {a: (f_t.get(a, 0.0) + f_o.get(a, 0.0)) / 2.0 for a in alleles}
    h_t = nei_gene_diversity(mean_freq)

    gst = (h_t - h_s) / h_t if h_t > 0 else 0.0
    denom = 1.0 - h_s
    jost_d = ((h_t - h_s) / denom) * 2.0 if denom > 0 else 0.0  # k=2
    shared = set(f_t) & set(f_o)
    return LocusDifferentiation(
        locus_id=mh.locus_id,
        target_n=n_t,
        offtarget_n=n_o,
        h_s=h_s,
        h_t=h_t,
        gst=max(0.0, min(1.0, gst)),
        jost_d=max(0.0, min(1.0, jost_d)),
        is_diagnostic=not shared,
    )


def genome_wide_fst(
    loci: Sequence[Microhaplotype],
    targets: Sequence[str],
    off_targets: Sequence[str],
) -> float:
    """Genome-wide G_ST across *loci* by ratio of averages: ΣΔ / ΣH_T.

    Returns 0.0 when no informative locus is present.
    """
    num = 0.0
    den = 0.0
    for mh in loci:
        d = locus_differentiation(mh, targets, off_targets)
        if d is None:
            continue
        num += d.h_t - d.h_s
        den += d.h_t
    return num / den if den > 0 else 0.0
