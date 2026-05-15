"""Project metadata regression tests."""

from __future__ import annotations

from pathlib import Path


def test_cli_dependency_keeps_typer_and_click_compatible() -> None:
    """Typer releases that use ``click.Choice[...]`` need modern Click."""
    pyproject = Path(__file__).resolve().parents[2] / "pyproject.toml"

    assert '"click>=8.2"' in pyproject.read_text(encoding="utf-8")
