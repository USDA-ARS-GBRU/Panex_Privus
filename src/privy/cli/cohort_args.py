"""Shared CLI parsing helpers for grouped cohort sample options."""

from __future__ import annotations

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
