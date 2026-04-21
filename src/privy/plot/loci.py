"""Per-locus ranked panel and main plot dispatcher for Panex Privus.

:func:`run_plot` is the single entry point called by ``privy plot``.
It dispatches to the appropriate plot function(s) based on *plot_type*
and the optional input files that are present.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

from privy.core.config import PrivyConfig

log = logging.getLogger("privy.plot.loci")


def plot_locus_panel(
    hits_rows: list[dict[str, str]],
    outdir: Path,
    top_n: int = 30,
    width: float = 10.0,
    height: float = 5.0,
    dpi: int = 150,
    output_format: str = "png",
    show_labels: bool = True,
) -> Path:
    """Ranked lollipop panel of top-N hits by ``final_score``.

    Each hit is drawn as a vertical stem terminating in a coloured dot.
    Dot colour maps to strictness class using the Panex Privus palette.

    Args:
        hits_rows: Rows from hits.tsv.
        outdir: Output directory.
        top_n: Number of top-ranked loci to display.
        width: Figure width in inches.
        height: Figure height in inches.
        dpi: Figure DPI.
        output_format: ``"png"``, ``"svg"``, or ``"pdf"``.
        show_labels: Annotate x-axis with locus IDs (only when top_n ≤ 20).

    Returns:
        Path to the written figure file.
    """
    import matplotlib  # noqa: PLC0415
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # noqa: PLC0415
    import matplotlib.lines as mlines  # noqa: PLC0415

    from privy.plot.themes import (  # noqa: PLC0415
        STRICTNESS_COLOURS,
        STRICTNESS_ORDER,
        apply_privy_theme,
    )

    apply_privy_theme()

    rows = sorted(
        hits_rows,
        key=lambda r: float(r.get("final_score", 0.0)),
        reverse=True,
    )[:top_n]

    fig, ax = plt.subplots(figsize=(width, height))

    if not rows:
        ax.text(0.5, 0.5, "No hits to display", ha="center", va="center",
                transform=ax.transAxes, fontsize=12, color="#888888")
        outpath = outdir / f"locus_panel.{output_format}"
        fig.savefig(outpath, dpi=dpi)
        plt.close(fig)
        return outpath

    x = list(range(1, len(rows) + 1))
    scores = [float(r.get("final_score", 0.0)) for r in rows]
    colors = [STRICTNESS_COLOURS.get(r.get("strictness_class", ""), "#aaaaaa") for r in rows]

    for xi, yi in zip(x, scores):
        ax.vlines(xi, 0, yi, color="#cccccc", linewidth=0.9, zorder=1)
    for xi, yi, color in zip(x, scores, colors):
        ax.scatter(xi, yi, color=color, s=45, zorder=2, edgecolors="none")

    ax.set_xlabel(f"Rank (top {len(rows)} by final_score)")
    ax.set_ylabel("final_score")
    ax.set_title("Top-scored private loci")
    ax.set_xlim(0.5, len(rows) + 0.5)
    ax.set_ylim(bottom=0.0)

    if show_labels and len(rows) <= 20:
        ax.set_xticks(x)
        ax.set_xticklabels(
            [r.get("locus_id", "") for r in rows],
            rotation=45, ha="right", fontsize=7,
        )
    else:
        ax.set_xticks([1, len(rows) // 2, len(rows)])

    present_classes = {r.get("strictness_class", "") for r in rows}
    handles = [
        mlines.Line2D(
            [], [], marker="o", color="w",
            markerfacecolor=STRICTNESS_COLOURS.get(sc, "#aaaaaa"),
            markersize=8, label=sc,
        )
        for sc in STRICTNESS_ORDER if sc in present_classes
    ]
    if handles:
        ax.legend(handles=handles, loc="upper right")

    outpath = outdir / f"locus_panel.{output_format}"
    fig.savefig(outpath, dpi=dpi)
    plt.close(fig)
    return outpath


def run_plot(
    hits: Path,
    regions: Optional[Path],
    evidence: Optional[Path],
    vcf: Optional[Path],
    bam: Optional[list[Path]],
    bam_manifest: Optional[Path],
    gfa: Optional[Path],
    xmfa: Optional[Path],  # accepted for CLI compat; not used
    compare: Optional[Path],
    cfg: PrivyConfig,
    locus_id: Optional[str],
    region_id: Optional[str],
    top_n: Optional[int],
    contig: Optional[str],
    region: Optional[str],
    plot_type: str,
    width: float,
    height: float,
    dpi: int,
    output_format: str,
    show_labels: bool,
    outdir: Path,
) -> list[Path]:
    """Generate focused plots and return the list of written file paths.

    The *plot_type* parameter controls which plots are generated:

    - ``"all"`` (default): every applicable plot given the available inputs.
    - ``"locus_panel"``: ranked lollipop of top-N hits.
    - ``"strictness_bar"``: strictness class distribution.
    - ``"score_distribution"``: final_score histogram by strictness class.
    - ``"support_bar"``: evidence class breakdown by source (requires --evidence).
    - ``"compare_summary"``: match class distribution (requires --compare).

    Args:
        hits: Path to hits.tsv (required).
        regions: Optional path to regions.tsv (not yet used by any plot type).
        evidence: Optional path to evidence.tsv (needed for support_bar).
        vcf: Accepted for CLI compat; not used in this version.
        bam: Accepted for CLI compat; not used in this version.
        bam_manifest: Accepted for CLI compat; not used in this version.
        gfa: Accepted for CLI compat; not used in this version.
        xmfa: Accepted for CLI compat; not used.
        compare: Optional path to compare.tsv (needed for compare_summary).
        cfg: Resolved PrivyConfig (not currently used but kept for future options).
        locus_id: Restrict to a single locus ID (not yet used).
        region_id: Restrict to a single region ID (not yet used).
        top_n: Number of top loci to show in locus_panel.
        contig: Restrict to one contig (not yet used).
        region: Restrict to a coordinate range (not yet used).
        plot_type: Which plot(s) to generate.
        width: Figure width in inches.
        height: Figure height in inches.
        dpi: Figure DPI.
        output_format: ``"png"``, ``"svg"``, or ``"pdf"``.
        show_labels: Show locus-id labels on the x-axis of locus_panel.
        outdir: Output directory.

    Returns:
        List of Path objects for each written figure file.
    """
    outdir.mkdir(parents=True, exist_ok=True)

    from privy.io.tsv import read_tsv  # noqa: PLC0415
    from privy.plot.summaries import (  # noqa: PLC0415
        plot_compare_summary,
        plot_score_distribution,
        plot_strictness_bar,
        plot_support_bar,
    )

    hits_rows = read_tsv(hits)
    effective_top_n = top_n if top_n is not None else 30
    do_all = (plot_type == "all")
    generated: list[Path] = []

    if do_all or plot_type == "locus_panel":
        generated.append(
            plot_locus_panel(
                hits_rows, outdir,
                top_n=effective_top_n,
                width=width, height=height,
                dpi=dpi, output_format=output_format,
                show_labels=show_labels,
            )
        )

    if do_all or plot_type == "strictness_bar":
        generated.append(
            plot_strictness_bar(
                hits_rows, outdir,
                width=width, height=min(height, 5.0),
                dpi=dpi, output_format=output_format,
            )
        )

    if do_all or plot_type == "score_distribution":
        generated.append(
            plot_score_distribution(
                hits_rows, outdir,
                width=width, height=height,
                dpi=dpi, output_format=output_format,
            )
        )

    if evidence is not None and evidence.exists():
        if do_all or plot_type == "support_bar":
            evidence_rows = read_tsv(evidence)
            generated.append(
                plot_support_bar(
                    evidence_rows, outdir,
                    width=width, height=min(height, 5.0),
                    dpi=dpi, output_format=output_format,
                )
            )
    elif plot_type == "support_bar":
        raise ValueError("--evidence is required for plot_type=support_bar")

    if compare is not None and compare.exists():
        if do_all or plot_type == "compare_summary":
            compare_rows = read_tsv(compare)
            generated.append(
                plot_compare_summary(
                    compare_rows, outdir,
                    width=width, height=min(height, 5.0),
                    dpi=dpi, output_format=output_format,
                )
            )
    elif plot_type == "compare_summary":
        raise ValueError("--compare is required for plot_type=compare_summary")

    log.info("Generated %d plot(s) in %s", len(generated), outdir)
    return generated
