"""Genomic relationship matrix (GRM) and dosage matrix for genomic prediction.

Builds the per-sample alt-allele dosage matrix (samples × microhaplotype loci,
0..ploidy) and the VanRaden (2008) genomic relationship matrix — the standard GP
inputs consumed by rrBLUP / BGLR / sommer.  Privy's job ends at producing inputs;
model fitting stays in those tools.

Polyploid-aware via the ploidy generalization (centre by ``ploidy·p``, scale by
``Σ ploidy·p(1−p)``).  Citation: VanRaden 2008 J Dairy Sci 91:4414; polyploid form
per Ashraf/Endelman.  See scratch/notes/70_popgen_estimators_math.md.  numpy-backed
(numpy is a core dependency).
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

from privy.microhap.model import Microhaplotype
from privy.polyploid.dosage import alt_dosage, observed_ploidy


@dataclass(frozen=True)
class DosageMatrix:
    """Samples × loci alt-allele dosage (0..ploidy); None marks a missing call."""

    samples: tuple[str, ...]
    locus_ids: tuple[str, ...]
    ploidy: int
    matrix: list[list[int | None]]   # rows = samples, cols = loci


def build_dosage_matrix(
    loci: Sequence[Microhaplotype],
    samples_to_paths: Mapping[str, Sequence[str]],
) -> DosageMatrix:
    """Build a samples × loci alt-allele dosage matrix from grouped haplotype paths."""
    samples = tuple(samples_to_paths)
    locus_ids = tuple(mh.locus_id for mh in loci)
    ploidy = max((observed_ploidy(p) for p in samples_to_paths.values()), default=0)
    matrix: list[list[int | None]] = []
    for sample in samples:
        paths = samples_to_paths[sample]
        matrix.append([alt_dosage(mh, paths) for mh in loci])
    return DosageMatrix(samples=samples, locus_ids=locus_ids, ploidy=ploidy, matrix=matrix)


def vanraden_grm(dm: DosageMatrix) -> tuple[tuple[str, ...], list[list[float]]]:
    """Compute the VanRaden genomic relationship matrix from a dosage matrix.

    Missing dosages are mean-imputed per locus (to ``ploidy·p̂``); monomorphic loci
    contribute nothing.  Returns ``(sample_ids, G)`` with G a labelled square matrix.

    Raises:
        ValueError: If the matrix has no samples or no loci.
    """
    import numpy as np  # noqa: PLC0415

    if not dm.samples or not dm.locus_ids:
        raise ValueError("dosage matrix must have at least one sample and one locus")
    ploidy = dm.ploidy or 1

    m = np.array(
        [[np.nan if v is None else float(v) for v in row] for row in dm.matrix],
        dtype=float,
    )
    # Per-locus alt allele frequency from non-missing calls.
    col_mean = np.nanmean(m, axis=0)
    p = np.divide(col_mean, ploidy, out=np.zeros_like(col_mean), where=ploidy != 0)
    # Mean-impute missing dosages to ploidy*p (the column mean).
    inds = np.where(np.isnan(m))
    m[inds] = np.take(col_mean, inds[1])

    centred = m - (ploidy * p)
    denom = float(np.sum(ploidy * p * (1.0 - p)))
    if denom <= 0:
        # No polymorphic loci: relationships are undefined; return zeros.
        n = len(dm.samples)
        return dm.samples, [[0.0] * n for _ in range(n)]
    grm = (centred @ centred.T) / denom
    return dm.samples, [[round(float(v), 6) for v in row] for row in grm]
