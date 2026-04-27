"""Pangenome analysis plots."""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from statistics import mean
from typing import Any


def plot_pangenome_growth(
    growth_rows: list[dict[str, Any]],
    outdir: Path,
    width: float = 9.0,
    height: float = 5.5,
    dpi: int = 150,
    output_format: str = "png",
) -> Path:
    """Plot mean feature growth curves with permutation intervals."""
    import matplotlib  # noqa: PLC0415
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # noqa: PLC0415
    import numpy as np  # noqa: PLC0415

    from privy.plot.themes import (  # noqa: PLC0415
        PANGENOME_GROUP_COLOURS,
        PANGENOME_GROUP_ORDER,
        apply_privy_theme,
    )

    apply_privy_theme()
    outdir.mkdir(parents=True, exist_ok=True)
    grouped: dict[str, dict[int, list[float]]] = defaultdict(lambda: defaultdict(list))
    for row in growth_rows:
        grouped[str(row["group"])][int(row["n"])].append(float(row["features"]))

    fig, ax = plt.subplots(figsize=(width, height))
    if not grouped:
        ax.text(0.5, 0.5, "No pangenome growth rows", ha="center", va="center",
                transform=ax.transAxes, fontsize=12, color="#888888")
    else:
        for group in PANGENOME_GROUP_ORDER:
            if group not in grouped:
                continue
            x = sorted(grouped[group])
            values = [grouped[group][n] for n in x]
            y = [mean(v) for v in values]
            lo = [float(np.quantile(v, 0.025)) for v in values]
            hi = [float(np.quantile(v, 0.975)) for v in values]
            color = PANGENOME_GROUP_COLOURS[group]
            ax.plot(x, y, color=color, linewidth=2.2, marker="o", markersize=3.5, label=group)
            ax.fill_between(x, lo, hi, color=color, alpha=0.16, linewidth=0)

    ax.set_xlabel("Number of samples")
    ax.set_ylabel("Observed features")
    ax.set_title("Pangenome growth")
    ax.legend(loc="upper left")
    ax.grid(axis="y", color="#d9d9d9", linestyle="--", linewidth=0.8, alpha=0.55)

    outpath = outdir / f"pangenome_growth.{output_format}"
    fig.savefig(outpath, dpi=dpi)
    plt.close(fig)
    return outpath


def plot_pangenome_coverage(
    coverage_rows: list[dict[str, Any]],
    outdir: Path,
    width: float = 9.0,
    height: float = 5.0,
    dpi: int = 150,
    output_format: str = "png",
) -> Path:
    """Plot feature coverage histograms for full, target, and off-target groups."""
    import matplotlib  # noqa: PLC0415
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # noqa: PLC0415

    from privy.plot.themes import (  # noqa: PLC0415
        PANGENOME_GROUP_COLOURS,
        PANGENOME_GROUP_ORDER,
        apply_privy_theme,
    )

    apply_privy_theme()
    outdir.mkdir(parents=True, exist_ok=True)
    grouped: dict[str, dict[int, int]] = defaultdict(dict)
    for row in coverage_rows:
        grouped[str(row["group"])][int(row["coverage"])] = int(row["n_features"])

    fig, ax = plt.subplots(figsize=(width, height))
    if not grouped:
        ax.text(0.5, 0.5, "No coverage rows", ha="center", va="center",
                transform=ax.transAxes, fontsize=12, color="#888888")
    else:
        for group in PANGENOME_GROUP_ORDER:
            if group not in grouped:
                continue
            x = sorted(grouped[group])
            y = [grouped[group][n] for n in x]
            ax.plot(
                x, y,
                color=PANGENOME_GROUP_COLOURS[group],
                linewidth=2.0,
                marker="o",
                markersize=3.5,
                label=group,
            )

    ax.set_xlabel("Samples containing feature")
    ax.set_ylabel("Number of features")
    ax.set_title("Pangenome feature coverage")
    ax.legend(loc="upper right")
    ax.grid(axis="y", color="#d9d9d9", linestyle="--", linewidth=0.8, alpha=0.55)

    outpath = outdir / f"pangenome_coverage.{output_format}"
    fig.savefig(outpath, dpi=dpi)
    plt.close(fig)
    return outpath


def plot_pangenome_composition(
    composition_rows: list[dict[str, Any]],
    outdir: Path,
    width: float = 8.0,
    height: float = 5.0,
    dpi: int = 150,
    output_format: str = "png",
) -> Path:
    """Plot core/accessory/private/absent composition by group."""
    import matplotlib  # noqa: PLC0415
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # noqa: PLC0415

    from privy.plot.themes import (  # noqa: PLC0415
        PANGENOME_CATEGORY_COLOURS,
        PANGENOME_CATEGORY_ORDER,
        PANGENOME_GROUP_ORDER,
        apply_privy_theme,
    )

    apply_privy_theme()
    outdir.mkdir(parents=True, exist_ok=True)
    counts: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for row in composition_rows:
        counts[str(row["group"])][str(row["category"])] = int(row["n_features"])

    fig, ax = plt.subplots(figsize=(width, height))
    groups = [g for g in PANGENOME_GROUP_ORDER if g in counts]
    if not groups:
        ax.text(0.5, 0.5, "No composition rows", ha="center", va="center",
                transform=ax.transAxes, fontsize=12, color="#888888")
    else:
        x = list(range(len(groups)))
        bottom = [0] * len(groups)
        for category in PANGENOME_CATEGORY_ORDER:
            values = [counts[group].get(category, 0) for group in groups]
            if not any(values):
                continue
            ax.bar(
                x,
                values,
                bottom=bottom,
                color=PANGENOME_CATEGORY_COLOURS[category],
                edgecolor="none",
                alpha=0.9,
                label=category,
            )
            bottom = [b + v for b, v in zip(bottom, values, strict=False)]
        ax.set_xticks(x)
        ax.set_xticklabels(groups)

    ax.set_ylabel("Number of features")
    ax.set_title("Pangenome composition")
    ax.legend(loc="upper right")
    ax.grid(axis="y", color="#d9d9d9", linestyle="--", linewidth=0.8, alpha=0.55)

    outpath = outdir / f"pangenome_composition.{output_format}"
    fig.savefig(outpath, dpi=dpi)
    plt.close(fig)
    return outpath


def plot_all_pangenome(
    coverage_rows: list[dict[str, Any]],
    composition_rows: list[dict[str, Any]],
    growth_rows: list[dict[str, Any]],
    outdir: Path,
    width: float = 9.0,
    height: float = 5.5,
    dpi: int = 150,
    output_format: str = "png",
) -> list[Path]:
    """Generate all first-pass pangenome plots."""
    return [
        plot_pangenome_growth(
            growth_rows, outdir,
            width=width, height=height, dpi=dpi, output_format=output_format,
        ),
        plot_pangenome_coverage(
            coverage_rows, outdir,
            width=width, height=min(height, 5.0), dpi=dpi, output_format=output_format,
        ),
        plot_pangenome_composition(
            composition_rows, outdir,
            width=min(width, 8.0), height=min(height, 5.0),
            dpi=dpi, output_format=output_format,
        ),
    ]
