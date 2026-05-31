"""Breeder-facing population genetics on microhaplotype alleles.

Turns Privy's multi-allelic loci into actionable summaries for plant breeders:
allelic diversity, cohort differentiation (which loci separate the target group
from the rest), and private-allele signal — all multi-allelic-aware and
target/off-target-aware.  Pure-Python (numpy/sklearn are optional extras used only
by the heavier structure/relationship methods).

Estimator definitions and citations: scratch/notes/70_popgen_estimators_math.md.
"""

from __future__ import annotations

from privy.popgen.differentiation import (
    LocusDifferentiation,
    genome_wide_fst,
    locus_differentiation,
)
from privy.popgen.diversity import (
    LocusDiversity,
    allele_frequencies,
    effective_n_alleles,
    inbreeding_fis,
    locus_diversity,
    nei_gene_diversity,
)
from privy.popgen.private import (
    private_allele_counts,
    rarefied_private_allelic_richness,
)
from privy.popgen.relationship import (
    DosageMatrix,
    build_dosage_matrix,
    vanraden_grm,
)
from privy.popgen.structure import (
    DapcResult,
    PcaResult,
    dapc,
    labels_from_cohorts,
    pca,
)

__all__ = [
    "DapcResult",
    "DosageMatrix",
    "LocusDifferentiation",
    "LocusDiversity",
    "PcaResult",
    "allele_frequencies",
    "build_dosage_matrix",
    "dapc",
    "effective_n_alleles",
    "genome_wide_fst",
    "inbreeding_fis",
    "labels_from_cohorts",
    "locus_differentiation",
    "locus_diversity",
    "nei_gene_diversity",
    "pca",
    "private_allele_counts",
    "rarefied_private_allelic_richness",
    "vanraden_grm",
]
