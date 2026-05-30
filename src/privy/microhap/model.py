"""Data model for microhaplotypes (multi-allelic local loci)."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field


@dataclass(frozen=True)
class Microhaplotype:
    """A multi-allelic locus: per-genome local alleles between shared flanks.

    Alleles are content-hashed (MD5 of the orientation-aware local sequence), the
    PHG/hVCF convention — identical sequence ⇒ identical allele id across genomes
    and runs.

    Attributes:
        locus_id: Stable id, e.g. ``chr1_0000000010``.
        contig / start / end: Reference span of the locus (0-based half-open).
        alleles: genome/path id → allele id (MD5 hex).  Genomes missing a flank
            are simply absent from this mapping.
        ref_allele: Allele id carried by the reference genome (None if absent).
        allele_tokens: allele id → human-readable token (segment list / sequence)
            for inspection.
    """

    locus_id: str
    contig: str
    start: int
    end: int
    alleles: Mapping[str, str]
    ref_allele: str | None = None
    allele_tokens: Mapping[str, str] = field(default_factory=dict)

    @property
    def n_genomes(self) -> int:
        """Number of genomes with an allele call here (non-missing)."""
        return len(self.alleles)

    @property
    def n_alleles(self) -> int:
        """Number of distinct alleles observed."""
        return len(set(self.alleles.values()))

    @property
    def is_multiallelic(self) -> bool:
        """True when more than one distinct allele is present."""
        return self.n_alleles > 1

    def allele_counts(self) -> dict[str, int]:
        """Count of genomes carrying each allele id."""
        counts: dict[str, int] = {}
        for allele in self.alleles.values():
            counts[allele] = counts.get(allele, 0) + 1
        return counts

    def allele_frequencies(self) -> dict[str, float]:
        """Allele id → frequency over genomes with data."""
        n = self.n_genomes
        if n == 0:
            return {}
        return {a: c / n for a, c in self.allele_counts().items()}

    def aaf(self) -> float:
        """Combined alternative allele frequency (1 − frequency of the reference allele).

        When the reference allele is unknown, the most common allele is treated as
        the reference proxy.  Returns 0.0 with no data.
        """
        freqs = self.allele_frequencies()
        if not freqs:
            return 0.0
        if self.ref_allele is not None and self.ref_allele in freqs:
            return 1.0 - freqs[self.ref_allele]
        return 1.0 - max(freqs.values())

    def private_alleles(
        self,
        targets: Sequence[str],
        off_targets: Sequence[str],
    ) -> list[str]:
        """Allele ids present in ≥1 target genome and 0 off-target genomes."""
        target_set = set(targets)
        offtarget_set = set(off_targets)
        in_target: set[str] = set()
        in_offtarget: set[str] = set()
        for genome, allele in self.alleles.items():
            if genome in target_set:
                in_target.add(allele)
            if genome in offtarget_set:
                in_offtarget.add(allele)
        return sorted(in_target - in_offtarget)

    def is_target_private(
        self,
        targets: Sequence[str],
        off_targets: Sequence[str],
    ) -> bool:
        """True when at least one allele is private to the target cohort."""
        return bool(self.private_alleles(targets, off_targets))
