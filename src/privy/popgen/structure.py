"""Population structure: PCA (Patterson scaling) and DAPC.

PCA on the per-sample dosage matrix uses the Patterson/EIGENSTRAT scaling
(centre by the locus mean, scale by ``sqrt(p(1-p))``) and is numpy-only (numpy is
a core dependency).  DAPC (Jombart 2010) is an optional add-on requiring
scikit-learn (Tier-1 extra) and degrades gracefully when it is absent.

Citations: Patterson, Price & Reich (2006) PLoS Genet 2:e190; Jombart, Devillard
& Balloux (2010) BMC Genet 11:94.  See scratch/notes/70_popgen_estimators_math.md.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

from privy.popgen.relationship import DosageMatrix
from privy.utils.optional import require


@dataclass(frozen=True)
class PcaResult:
    """PCA coordinates per sample plus explained-variance ratios."""

    samples: tuple[str, ...]
    coords: list[list[float]]              # n_samples × n_components
    explained_variance_ratio: list[float]
    n_components: int


def _imputed_array(dm: DosageMatrix) -> tuple[Any, Any, int]:
    """Return (numpy array with missing mean-imputed, column means, ploidy)."""
    import numpy as np  # noqa: PLC0415

    m = np.array(
        [[np.nan if v is None else float(v) for v in row] for row in dm.matrix],
        dtype=float,
    )
    if m.size == 0:
        return m, np.zeros(0), max(dm.ploidy, 1)
    col_mean = np.nanmean(m, axis=0)
    inds = np.where(np.isnan(m))
    m[inds] = np.take(col_mean, inds[1])
    return m, col_mean, max(dm.ploidy, 1)


def pca(dm: DosageMatrix, n_components: int = 2) -> PcaResult:
    """Patterson-scaled PCA of samples over loci.

    Monomorphic loci (zero variance) are dropped.  ``n_components`` is capped at
    the rank of the scaled matrix.  Returns empty coords when no locus is informative.
    """
    import numpy as np  # noqa: PLC0415

    m, col_mean, ploidy = _imputed_array(dm)
    if m.size == 0:
        return PcaResult(dm.samples, [], [], 0)

    p = col_mean / ploidy
    scale = np.sqrt(p * (1.0 - p))
    keep = scale > 0
    if not keep.any():
        return PcaResult(dm.samples, [[] for _ in dm.samples], [], 0)
    z = (m[:, keep] - col_mean[keep]) / scale[keep]

    # SVD; principal coordinates = U * S.
    u, s, _vt = np.linalg.svd(z, full_matrices=False)
    k = min(n_components, s.shape[0], (s > 1e-12).sum())
    if k == 0:
        return PcaResult(dm.samples, [[] for _ in dm.samples], [], 0)
    coords = (u[:, :k] * s[:k])
    total = float(np.sum(s**2))
    evr = [float(v) for v in (s[:k] ** 2 / total)] if total > 0 else [0.0] * k
    return PcaResult(
        samples=dm.samples,
        coords=[[round(float(v), 6) for v in row] for row in coords],
        explained_variance_ratio=[round(v, 6) for v in evr],
        n_components=k,
    )


@dataclass(frozen=True)
class DapcResult:
    """DAPC discriminant coordinates + per-sample assignment."""

    samples: tuple[str, ...]
    coords: list[list[float]]
    assigned: list[str]


def dapc(
    dm: DosageMatrix,
    labels: Mapping[str, str],
    *,
    n_pcs: int | None = None,
) -> DapcResult:
    """Discriminant Analysis of Principal Components (PCA → LDA on group labels).

    Requires scikit-learn (optional Tier-1 extra).

    Args:
        labels: sample id → group label (samples missing a label are dropped).
        n_pcs: PCs retained before LDA (default: min(n_samples-1, n_loci, 10)).

    Raises:
        MissingDependencyError: If scikit-learn is not installed.
        ValueError: If fewer than two labelled groups are present.
    """
    require("sklearn", feature="DAPC population structure")
    import numpy as np  # noqa: PLC0415
    from sklearn.decomposition import PCA  # noqa: PLC0415
    from sklearn.discriminant_analysis import LinearDiscriminantAnalysis  # noqa: PLC0415

    m, _col_mean, _ploidy = _imputed_array(dm)
    idx = [i for i, s in enumerate(dm.samples) if s in labels]
    if len(idx) < 2:
        raise ValueError("DAPC needs at least two labelled samples")
    y = [labels[dm.samples[i]] for i in idx]
    if len(set(y)) < 2:
        raise ValueError("DAPC needs at least two distinct groups")
    x = m[idx, :]

    n_keep = n_pcs if n_pcs is not None else min(len(idx) - 1, x.shape[1], 10)
    n_keep = max(1, n_keep)
    pcs = PCA(n_components=n_keep).fit_transform(x)
    lda = LinearDiscriminantAnalysis()
    coords = lda.fit_transform(pcs, y)
    assigned = list(lda.predict(pcs))
    return DapcResult(
        samples=tuple(dm.samples[i] for i in idx),
        coords=[[round(float(v), 6) for v in row] for row in np.atleast_2d(coords)],
        assigned=[str(a) for a in assigned],
    )


def labels_from_cohorts(
    targets: Sequence[str],
    off_targets: Sequence[str],
    delimiter: str = "#",
) -> dict[str, str]:
    """Build a sample→label map ("target"/"offtarget") from cohort path ids."""
    from privy.synteny.model import split_pansn  # noqa: PLC0415

    labels: dict[str, str] = {}
    for path in targets:
        labels[split_pansn(path, delimiter)[0]] = "target"
    for path in off_targets:
        labels.setdefault(split_pansn(path, delimiter)[0], "offtarget")
    return labels
