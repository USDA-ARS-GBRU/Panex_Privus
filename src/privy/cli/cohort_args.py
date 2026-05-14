"""Shared CLI helpers for cohort sample options and cohort files."""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

GROUPED_COHORT_FLAGS = {
    "--targets": "targets",
    "--target": "targets",
    "--off-targets": "off_targets",
    "--off-target": "off_targets",
    "--ignore-samples": "ignore_samples",
    "--ignore-sample": "ignore_samples",
}

GROUPED_COHORT_DISPLAY = {
    "targets": "--targets",
    "off_targets": "--off-targets",
    "ignore_samples": "--ignore-samples",
}


@dataclass(frozen=True)
class CohortSampleSpec:
    """Loose cohort samples parsed from CLI-facing inputs.

    Unlike :class:`privy.core.cohort.CohortDefinition`, this object allows
    missing off-targets so context modules can infer them from the input.
    Discovery commands still validate that both active groups are complete.
    """

    targets: tuple[str, ...] = field(default_factory=tuple)
    off_targets: tuple[str, ...] = field(default_factory=tuple)
    ignored_samples: tuple[str, ...] = field(default_factory=tuple)

    @classmethod
    def from_values(
        cls,
        *,
        targets: list[str] | tuple[str, ...] | None = None,
        off_targets: list[str] | tuple[str, ...] | None = None,
        ignored_samples: list[str] | tuple[str, ...] | None = None,
    ) -> CohortSampleSpec:
        """Build a cohort sample spec while preserving order and removing duplicates."""
        return cls(
            targets=tuple(_dedupe_samples(targets or ())),
            off_targets=tuple(_dedupe_samples(off_targets or ())),
            ignored_samples=tuple(_dedupe_samples(ignored_samples or ())),
        )


def parse_grouped_cohort_args(args: list[str]) -> dict[str, list[str] | None]:
    """Parse grouped cohort sample options left over after Click parsing.

    Typer's native ``list[str]`` option handling requires users to repeat the
    flag for every value. Privy's public CLI accepts the friendlier grouped
    form, for example ``--targets T1 T2 --off-targets O1 O2``.
    """
    values: dict[str, list[str]] = {
        "targets": [],
        "off_targets": [],
        "ignore_samples": [],
    }
    provided: set[str] = set()
    active_group: str | None = None

    for token in args:
        if token.startswith("--"):
            flag, separator, inline_value = token.partition("=")
            group = GROUPED_COHORT_FLAGS.get(flag)
            if group is None:
                raise ValueError(f"No such option: {flag}")

            active_group = group
            provided.add(group)
            if separator:
                if not inline_value:
                    raise ValueError(f"{flag} requires at least one sample name.")
                values[group].append(inline_value)
            continue

        if active_group is None:
            raise ValueError(
                f"Unexpected argument {token!r}. Sample names must follow "
                "--targets, --off-targets, or --ignore-samples."
            )
        values[active_group].append(token)

    for group in provided:
        if not values[group]:
            raise ValueError(
                f"{GROUPED_COHORT_DISPLAY[group]} requires at least one sample name."
            )

    return {
        group: values[group] if group in provided else None
        for group in values
    }


def collect_sample_values(
    values: list[str] | None,
    path: Path | None = None,
) -> list[str] | None:
    """Merge grouped CLI values and an optional one-sample-per-line file.

    ``None`` means the role was not supplied at all. An empty list means the
    role was supplied explicitly, but no samples were read after stripping blank
    and comment lines.
    """
    if values is None and path is None:
        return None

    merged: list[str] = []
    for value in values or []:
        sample = value.strip()
        if sample:
            merged.append(sample)

    if path is not None:
        merged.extend(read_sample_list(path))

    return _dedupe_samples(merged)


def load_cohort_file(path: Path) -> CohortSampleSpec:
    """Load a cohort sample specification from YAML or TSV."""
    if not path.exists():
        raise FileNotFoundError(f"Cohort file not found: {path}")

    suffix = path.suffix.lower()
    if suffix in {".yaml", ".yml"}:
        return _load_cohort_yaml(path)
    if suffix == ".tsv":
        return _load_cohort_tsv(path)
    raise ValueError(
        f"Unsupported cohort file format: {path.suffix!r}. Use .yaml, .yml, or .tsv."
    )


def read_sample_list(path: Path) -> list[str]:
    """Read one sample name per line, ignoring blank lines and comments."""
    if not path.exists():
        raise FileNotFoundError(f"Sample list file not found: {path}")

    samples: list[str] = []
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        samples.append(line.split()[0])
    return samples


def _load_cohort_yaml(path: Path) -> CohortSampleSpec:
    """Load a cohort sample specification from YAML."""
    with open(path, encoding="utf-8") as fh:
        raw = yaml.safe_load(fh) or {}

    if not isinstance(raw, dict):
        raise ValueError("Cohort YAML must contain a mapping.")

    cohort_raw: dict[str, Any] = raw
    nested = cohort_raw.get("cohorts")
    if isinstance(nested, dict):
        cohort_raw = nested

    return CohortSampleSpec.from_values(
        targets=_as_sample_list(
            _first_present(cohort_raw, ("targets", "target")),
            field="targets",
        ),
        off_targets=_as_sample_list(
            _first_present(
                cohort_raw,
                ("off_targets", "off-targets", "offtargets", "background"),
            ),
            field="off_targets",
        ),
        ignored_samples=_as_sample_list(
            _first_present(
                cohort_raw,
                ("ignored_samples", "ignore_samples", "ignored", "ignore"),
            ),
            field="ignored_samples",
        ),
    )


def _load_cohort_tsv(path: Path) -> CohortSampleSpec:
    """Load a cohort sample specification from a TSV with sample/role columns."""
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Cohort TSV is missing a header row.")

        fieldnames = set(reader.fieldnames)
        sample_field = _first_available_field(
            fieldnames,
            ("sample_id", "sample", "sample_name"),
        )
        role_field = _first_available_field(
            fieldnames,
            ("cohort_role", "role", "group"),
        )

        if sample_field is None or role_field is None:
            raise ValueError(
                "Cohort TSV must contain sample_id/sample and cohort_role/role columns."
            )

        targets: list[str] = []
        off_targets: list[str] = []
        ignored: list[str] = []

        for row in reader:
            sample = (row.get(sample_field) or "").strip()
            role = _normalize_role(row.get(role_field) or "")
            if not sample:
                continue
            if role in {"target", "targets"}:
                targets.append(sample)
            elif role in {
                "off_target",
                "off_targets",
                "offtarget",
                "offtargets",
                "background",
                "backgrounds",
                "comparison",
                "comparisons",
            }:
                off_targets.append(sample)
            elif role in {"ignored", "ignore", "excluded", "exclude"}:
                ignored.append(sample)
            else:
                raise ValueError(
                    f"Unsupported cohort role {role!r} for sample {sample!r} in {path}."
                )

    return CohortSampleSpec.from_values(
        targets=targets,
        off_targets=off_targets,
        ignored_samples=ignored,
    )


def _first_present(raw: dict[str, Any], keys: tuple[str, ...]) -> Any:
    """Return the first present mapping value for one of several accepted keys."""
    for key in keys:
        if key in raw:
            return raw[key]
    return []


def _as_sample_list(value: Any, *, field: str) -> list[str]:
    """Normalize a YAML cohort value into a sample-name list."""
    if value is None:
        return []
    if isinstance(value, str):
        sample = value.strip()
        return [sample] if sample else []
    if not isinstance(value, (list, tuple)):
        raise ValueError(f"Cohort YAML field {field!r} must be a string or list.")

    samples: list[str] = []
    for item in value:
        sample = str(item).strip()
        if sample:
            samples.append(sample)
    return samples


def _first_available_field(
    fieldnames: set[str],
    candidates: tuple[str, ...],
) -> str | None:
    """Return the first available field name from a set of TSV aliases."""
    for candidate in candidates:
        if candidate in fieldnames:
            return candidate
    return None


def _normalize_role(role: str) -> str:
    """Normalize cohort role spelling across CLI cohort TSVs."""
    return role.strip().lower().replace("-", "_").replace(" ", "_")


def _dedupe_samples(values: list[str] | tuple[str, ...]) -> list[str]:
    """Return non-empty sample names in first-seen order."""
    return list(dict.fromkeys(sample for sample in values if sample))
