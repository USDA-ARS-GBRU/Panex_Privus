"""Summary diagnostic plot generators for Panex Privus.

Each function accepts pre-parsed TSV rows (list of dicts), writes one figure
to *outdir*, and returns the output path.  All functions use the Agg backend
so they are safe in headless/CI environments.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from pathlib import Path


def plot_strictness_bar(
    hits_rows: list[dict[str, str]],
    outdir: Path,
    width: float = 8.0,
    height: float = 4.0,
    dpi: int = 150,
    output_format: str = "png",
) -> Path:
    """Horizontal bar chart of strictness class distribution.

    Args:
        hits_rows: Rows from hits.tsv.
        outdir: Output directory.
        width: Figure width in inches.
        height: Figure height in inches.
        dpi: Figure DPI for raster formats.
        output_format: ``"png"``, ``"svg"``, or ``"pdf"``.

    Returns:
        Path to the written figure file.
    """
    import matplotlib  # noqa: PLC0415
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # noqa: PLC0415

    from privy.plot.themes import (  # noqa: PLC0415
        STRICTNESS_COLOURS,
        STRICTNESS_ORDER,
        apply_privy_theme,
    )

    apply_privy_theme()

    counts: Counter[str] = Counter(
        r.get("strictness_class", "unknown") for r in hits_rows
    )
    labels = [s for s in STRICTNESS_ORDER if counts.get(s, 0) > 0]
    if not labels:
        labels = sorted(counts.keys())
    values = [counts[s] for s in labels]
    colors = [STRICTNESS_COLOURS.get(s, "#aaaaaa") for s in labels]

    fig, ax = plt.subplots(figsize=(width, height))
    y_pos = list(range(len(labels)))
    bars = ax.barh(y_pos, values, color=colors, edgecolor="none", height=0.6)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels)
    ax.set_xlabel("Number of loci")
    ax.set_title(f"Strictness class distribution  (n={len(hits_rows)})")
    ax.invert_yaxis()

    for bar, val in zip(bars, values, strict=False):
        ax.text(
            bar.get_width() + max(values) * 0.01,
            bar.get_y() + bar.get_height() / 2,
            str(val), va="center", ha="left", fontsize=9,
        )

    outpath = outdir / f"strictness_bar.{output_format}"
    fig.savefig(outpath, dpi=dpi)
    plt.close(fig)
    return outpath


def plot_score_distribution(
    hits_rows: list[dict[str, str]],
    outdir: Path,
    width: float = 8.0,
    height: float = 5.0,
    dpi: int = 150,
    output_format: str = "png",
) -> Path:
    """Stacked histogram of ``final_score`` values coloured by strictness class.

    Args:
        hits_rows: Rows from hits.tsv.
        outdir: Output directory.
        width: Figure width in inches.
        height: Figure height in inches.
        dpi: Figure DPI.
        output_format: ``"png"``, ``"svg"``, or ``"pdf"``.

    Returns:
        Path to the written figure file.
    """
    import matplotlib  # noqa: PLC0415
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # noqa: PLC0415
    import numpy as np  # noqa: PLC0415

    from privy.plot.themes import (  # noqa: PLC0415
        STRICTNESS_COLOURS,
        STRICTNESS_ORDER,
        apply_privy_theme,
    )

    apply_privy_theme()

    scores_by_class: dict[str, list[float]] = defaultdict(list)
    for row in hits_rows:
        sc = row.get("strictness_class", "unknown")
        try:
            scores_by_class[sc].append(float(row.get("final_score", 0.0)))
        except ValueError:
            pass

    all_scores = [s for vals in scores_by_class.values() for s in vals]
    if not all_scores:
        fig, ax = plt.subplots(figsize=(width, height))
        ax.text(0.5, 0.5, "No scores to display", ha="center", va="center",
                transform=ax.transAxes, fontsize=12, color="#888888")
        outpath = outdir / f"score_distribution.{output_format}"
        fig.savefig(outpath, dpi=dpi)
        plt.close(fig)
        return outpath

    bin_max = max(all_scores) + 0.05
    bins = np.linspace(0.0, bin_max, min(30, len(all_scores) + 2))

    fig, ax = plt.subplots(figsize=(width, height))
    bottom = np.zeros(len(bins) - 1)

    for sc in STRICTNESS_ORDER:
        vals = scores_by_class.get(sc, [])
        if not vals:
            continue
        hist, _ = np.histogram(vals, bins=bins)
        color = STRICTNESS_COLOURS.get(sc, "#aaaaaa")
        ax.bar(
            bins[:-1], hist, width=np.diff(bins), bottom=bottom,
            color=color, edgecolor="none", alpha=0.88, label=sc, align="edge",
        )
        bottom = bottom + hist

    ax.set_xlabel("final_score")
    ax.set_ylabel("Number of loci")
    ax.set_title(f"Score distribution by strictness class  (n={len(all_scores)})")
    ax.legend(loc="upper left")

    outpath = outdir / f"score_distribution.{output_format}"
    fig.savefig(outpath, dpi=dpi)
    plt.close(fig)
    return outpath


def plot_support_bar(
    evidence_rows: list[dict[str, str]],
    outdir: Path,
    width: float = 8.0,
    height: float = 4.0,
    dpi: int = 150,
    output_format: str = "png",
) -> Path:
    """Stacked bar chart of evidence class distribution by source type.

    Args:
        evidence_rows: Rows from evidence.tsv.
        outdir: Output directory.
        width: Figure width in inches.
        height: Figure height in inches.
        dpi: Figure DPI.
        output_format: ``"png"``, ``"svg"``, or ``"pdf"``.

    Returns:
        Path to the written figure file.
    """
    import matplotlib  # noqa: PLC0415
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # noqa: PLC0415

    from privy.plot.themes import EVIDENCE_COLOURS, apply_privy_theme  # noqa: PLC0415

    apply_privy_theme()

    evidence_order = ["support", "absence", "ambiguous", "contradiction", "uninformative"]
    counts: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for row in evidence_rows:
        src = row.get("source_type", "unknown")
        ec = row.get("evidence_class", "unknown")
        counts[src][ec] += 1

    sources = sorted(counts.keys())
    if not sources:
        fig, ax = plt.subplots(figsize=(width, height))
        ax.text(0.5, 0.5, "No evidence records", ha="center", va="center",
                transform=ax.transAxes, fontsize=12, color="#888888")
        outpath = outdir / f"support_bar.{output_format}"
        fig.savefig(outpath, dpi=dpi)
        plt.close(fig)
        return outpath

    fig, ax = plt.subplots(figsize=(width, height))
    x = list(range(len(sources)))
    bottom = [0] * len(sources)

    for ec in evidence_order:
        vals = [counts[src].get(ec, 0) for src in sources]
        if not any(vals):
            continue
        color = EVIDENCE_COLOURS.get(ec, "#aaaaaa")
        ax.bar(
            x, vals, bottom=bottom, color=color, edgecolor="none",
            alpha=0.88, label=ec,
        )
        bottom = [b + v for b, v in zip(bottom, vals, strict=False)]

    ax.set_xticks(x)
    ax.set_xticklabels(sources)
    ax.set_ylabel("Evidence records")
    ax.set_title("Evidence class distribution by source")
    ax.legend(loc="upper right")

    outpath = outdir / f"support_bar.{output_format}"
    fig.savefig(outpath, dpi=dpi)
    plt.close(fig)
    return outpath


def plot_compare_summary(
    compare_rows: list[dict[str, str]],
    outdir: Path,
    width: float = 8.0,
    height: float = 4.0,
    dpi: int = 150,
    output_format: str = "png",
) -> Path:
    """Horizontal bar chart of comparison match class distribution.

    Args:
        compare_rows: Rows from compare.tsv.
        outdir: Output directory.
        width: Figure width in inches.
        height: Figure height in inches.
        dpi: Figure DPI.
        output_format: ``"png"``, ``"svg"``, or ``"pdf"``.

    Returns:
        Path to the written figure file.
    """
    import matplotlib  # noqa: PLC0415
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # noqa: PLC0415

    from privy.plot.themes import MATCH_COLOURS, MATCH_ORDER, apply_privy_theme  # noqa: PLC0415

    apply_privy_theme()

    counts: Counter[str] = Counter(
        r.get("match_class", "unknown") for r in compare_rows
    )
    labels = [m for m in MATCH_ORDER if counts.get(m, 0) > 0]
    if not labels:
        labels = sorted(counts.keys())
    values = [counts[m] for m in labels]
    colors = [MATCH_COLOURS.get(m, "#aaaaaa") for m in labels]

    fig, ax = plt.subplots(figsize=(width, height))
    y_pos = list(range(len(labels)))
    bars = ax.barh(y_pos, values, color=colors, edgecolor="none", height=0.6)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels)
    ax.set_xlabel("Number of loci")
    ax.set_title(f"Comparison match class distribution  (n={len(compare_rows)})")
    ax.invert_yaxis()

    for bar, val in zip(bars, values, strict=False):
        ax.text(
            bar.get_width() + max(values) * 0.01,
            bar.get_y() + bar.get_height() / 2,
            str(val), va="center", ha="left", fontsize=9,
        )

    outpath = outdir / f"compare_summary.{output_format}"
    fig.savefig(outpath, dpi=dpi)
    plt.close(fig)
    return outpath
