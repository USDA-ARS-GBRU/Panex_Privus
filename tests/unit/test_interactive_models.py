"""Unit tests for interactive dashboard models."""

from __future__ import annotations

import pytest

from privy.interactive.models import parse_focus_region


def test_parse_focus_region_accepts_commas() -> None:
    region = parse_focus_region("Gm15:1-4,000,000")

    assert region.contig == "Gm15"
    assert region.start == 1
    assert region.end == 4_000_000
    assert region.length == 4_000_000
    assert region.slug == "focus_Gm15_1_4000000"


@pytest.mark.parametrize(
    "value",
    [
        "Gm15",
        "Gm15:0-10",
        "Gm15:200-100",
        ":1-100",
    ],
)
def test_parse_focus_region_rejects_invalid_values(value: str) -> None:
    with pytest.raises(ValueError):
        parse_focus_region(value)
