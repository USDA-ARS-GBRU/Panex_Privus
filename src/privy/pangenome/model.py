"""Shared data model for pangenome analyses.

The central abstraction is a feature-by-sample presence matrix.  A feature can
come from a GFA segment today or a VCF allele in a later adapter; once it is in
this shape, the same cohort-level analyses can be reused.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class FeatureRecord:
    """One analyzable pangenome feature."""

    feature_id: str
    source_type: str
    feature_type: str
    length: int = 1
    contig: str | None = None
    start: int | None = None
    end: int | None = None


@dataclass(frozen=True)
class FeatureMatrix:
    """Feature-by-sample presence matrix in sparse set form."""

    source_type: str
    features: tuple[FeatureRecord, ...]
    samples: tuple[str, ...]
    presence: dict[str, frozenset[str]]

    def samples_for_feature(self, feature_id: str) -> frozenset[str]:
        """Return samples where *feature_id* is present."""
        return self.presence.get(feature_id, frozenset())


@dataclass(frozen=True)
class PangenomeGroups:
    """Built-in pangenome groups used by Privy."""

    full: tuple[str, ...]
    target: tuple[str, ...]
    off_target: tuple[str, ...]
    ignored: tuple[str, ...] = ()

    def as_dict(self) -> dict[str, tuple[str, ...]]:
        """Return groups in output order."""
        return {
            "full": self.full,
            "target": self.target,
            "off_target": self.off_target,
        }


def resolve_pangenome_groups(
    all_samples: list[str] | tuple[str, ...] | set[str],
    targets: list[str] | tuple[str, ...],
    off_targets: list[str] | tuple[str, ...] | None = None,
    ignored_samples: list[str] | tuple[str, ...] | None = None,
) -> PangenomeGroups:
    """Resolve target/off-target pangenome groups from input samples.

    If targets are provided but off-targets are omitted, all non-target,
    non-ignored samples in the input become off-targets.
    """
    ignored = tuple(dict.fromkeys(ignored_samples or ()))
    ignored_set = set(ignored)
    all_ordered = tuple(s for s in dict.fromkeys(all_samples) if s not in ignored_set)
    target = tuple(s for s in dict.fromkeys(targets) if s not in ignored_set)

    if not target:
        raise ValueError("At least one target sample is required for pangenome analysis.")

    unknown_targets = sorted(set(target) - set(all_ordered))
    if unknown_targets:
        raise ValueError(
            "Target samples were not found in the input: "
            + ", ".join(unknown_targets)
        )

    if off_targets:
        off_target = tuple(s for s in dict.fromkeys(off_targets) if s not in ignored_set)
        unknown_offtargets = sorted(set(off_target) - set(all_ordered))
        if unknown_offtargets:
            raise ValueError(
                "Off-target samples were not found in the input: "
                + ", ".join(unknown_offtargets)
            )
    else:
        target_set = set(target)
        off_target = tuple(s for s in all_ordered if s not in target_set)

    overlap = set(target) & set(off_target)
    if overlap:
        raise ValueError(
            "Samples appear in both target and off-target groups: "
            + ", ".join(sorted(overlap))
        )
    if not off_target:
        raise ValueError(
            "At least one off-target sample is required. Provide --off-targets or "
            "include non-target samples in the input so Privy can infer them."
        )

    return PangenomeGroups(
        full=tuple(s for s in all_ordered if s in set(target) | set(off_target)),
        target=target,
        off_target=off_target,
        ignored=ignored,
    )
