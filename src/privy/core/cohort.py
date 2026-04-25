"""CohortDefinition — the biological grouping model for Panex Privus.

A :class:`CohortDefinition` is the authoritative record of which samples are
targets, which are off-targets, and which to ignore.  It is constructed once
per run and passed to all backends.

Key invariants:
    - A sample cannot appear in both ``targets`` and ``off_targets``.
    - Both ``targets`` and ``off_targets`` must be non-empty.
    - Samples in ``ignored_samples`` are excluded from all analyses.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class CohortDefinition:
    """Defines the biological grouping of samples for a Panex Privus run.

    Attributes:
        targets: Sample names that form the focal (target) cohort.
        off_targets: Sample names that form the background/comparison cohort.
        ignored_samples: Samples present in an input file but excluded from
            all analyses.
        metadata: Optional key-value annotations (e.g., species, project).

    Raises:
        ValueError: If ``targets`` and ``off_targets`` overlap, or either is
            empty.
    """

    targets: tuple[str, ...]
    off_targets: tuple[str, ...]
    ignored_samples: tuple[str, ...] = field(default_factory=tuple)
    metadata: dict[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        overlap = set(self.targets) & set(self.off_targets)
        if overlap:
            raise ValueError(
                f"Samples appear in both targets and off_targets: {sorted(overlap)}"
            )
        if not self.targets:
            raise ValueError("targets must contain at least one sample name.")
        if not self.off_targets:
            raise ValueError("off_targets must contain at least one sample name.")

    # ----------------------------------------------------------------- factory

    @classmethod
    def from_lists(
        cls,
        targets: list[str],
        off_targets: list[str],
        ignored_samples: list[str] | None = None,
        metadata: dict[str, str] | None = None,
    ) -> CohortDefinition:
        """Construct from plain Python lists (convenience wrapper)."""
        return cls(
            targets=tuple(targets),
            off_targets=tuple(off_targets),
            ignored_samples=tuple(ignored_samples or []),
            metadata=metadata or {},
        )

    # ----------------------------------------------------------------- lookups

    def is_target(self, sample: str) -> bool:
        """Return True if *sample* is in the target cohort."""
        return sample in self.targets

    def is_off_target(self, sample: str) -> bool:
        """Return True if *sample* is in the off-target cohort."""
        return sample in self.off_targets

    def is_ignored(self, sample: str) -> bool:
        """Return True if *sample* is explicitly ignored."""
        return sample in self.ignored_samples

    def is_known(self, sample: str) -> bool:
        """Return True if *sample* is assigned to any role."""
        return (
            sample in self.targets
            or sample in self.off_targets
            or sample in self.ignored_samples
        )

    # ---------------------------------------------------------------- computed

    @property
    def n_targets(self) -> int:
        """Number of target samples."""
        return len(self.targets)

    @property
    def n_off_targets(self) -> int:
        """Number of off-target samples."""
        return len(self.off_targets)

    @property
    def all_active_samples(self) -> frozenset[str]:
        """All samples that participate in analysis (not ignored)."""
        return frozenset(self.targets) | frozenset(self.off_targets)

    # ------------------------------------------------------------------- repr

    def __repr__(self) -> str:
        return (
            f"CohortDefinition("
            f"n_targets={self.n_targets}, "
            f"n_off_targets={self.n_off_targets}, "
            f"n_ignored={len(self.ignored_samples)})"
        )
