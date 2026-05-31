"""Polyploid / homeolog handling for Panex Privus.

A sample in a polyploid (or simply phased) pangenome contributes several
haplotype paths.  This package groups those paths back into per-sample genotypes
and computes allele dosage (0..ploidy) at microhaplotype loci — the substrate for
dosage-aware genomic prediction, observed heterozygosity, and homeolog-aware
analysis in crops.
"""

from __future__ import annotations

from privy.polyploid.dosage import (
    alt_dosage,
    group_paths_by_sample,
    is_heterozygous,
    observed_heterozygosity,
    observed_ploidy,
    sample_allele_dosage,
)

__all__ = [
    "alt_dosage",
    "group_paths_by_sample",
    "is_heterozygous",
    "observed_heterozygosity",
    "observed_ploidy",
    "sample_allele_dosage",
]
