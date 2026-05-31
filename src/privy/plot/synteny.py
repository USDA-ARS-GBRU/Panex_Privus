"""Static comparative-synteny figures: riparian braids and dotplots.

Renders the synteny artifacts written by ``privy synteny`` (``synteny_blocks.tsv``)
into publication-quality static figures with matplotlib — the static counterpart
to the interactive dashboard (P5).  Blocks are coloured by type
(collinear / inversion / translocation / duplication), the single biggest reason
riparian plots read well, and the layout separates DATA (the block rows) from
RENDER so figures can be re-skinned without recomputation.

All functions follow the Privy plot convention: row dicts in, a saved figure
:class:`~pathlib.Path` out.  Coordinates are 0-based half-open bp.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Any


def _block_rows_to_floats(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Coerce the numeric columns of synteny_blocks rows to ints."""
    out: list[dict[str, Any]] = []
    for row in rows:
        out.append({
            **row,
            "query_genome": str(row["query_genome"]),
            "query_start": int(row["query_start"]),
            "query_end": int(row["query_end"]),
            "ref_genome": str(row["ref_genome"]),
            "ref_contig": str(row["ref_contig"]),
            "ref_start": int(row["ref_start"]),
            "ref_end": int(row["ref_end"]),
            "block_type": str(row["block_type"]),
        })
    return out


def plot_riparian(
    block_rows: list[dict[str, Any]],
    outdir: Path,
    *,
    private_region_ids: set[str] | None = None,
    width: float = 11.0,
    height: float = 6.0,
    dpi: int = 150,
    output_format: str = "png",
) -> Path:
    """Render a static riparian braid plot from synteny block rows.

    The reference genome is the common ``ref_genome`` track at the bottom; each
    query genome is stacked above it.  Every block is drawn as a filled ribbon
    from its query-track interval to its reference-track interval, coloured by
    block type (slanted/crossing ribbons reveal rearrangements).
    """
    import matplotlib  # noqa: PLC0415
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # noqa: PLC0415
    from matplotlib.patches import Patch, Rectangle  # noqa: PLC0415

    from privy.plot.themes import (  # noqa: PLC0415
        SYNTENY_BLOCK_COLOURS,
        SYNTENY_BLOCK_ORDER,
        apply_privy_theme,
    )

    apply_privy_theme()
    outdir.mkdir(parents=True, exist_ok=True)
    rows = _block_rows_to_floats(block_rows)

    fig, ax = plt.subplots(figsize=(width, height))
    if not rows:
        ax.text(0.5, 0.5, "No synteny blocks", ha="center", va="center",
                transform=ax.transAxes, fontsize=12, color="#888888")
        return _save(fig, plt, outdir, "riparian", output_format, dpi)

    ref_genome = rows[0]["ref_genome"]
    query_genomes = sorted({r["query_genome"] for r in rows})
    # reference track y=0; queries stacked above
    track_y = {ref_genome: 0.0}
    for i, g in enumerate(query_genomes, start=1):
        track_y[g] = float(i)
    bar_h = 0.18

    seen_types: set[str] = set()
    for r in rows:
        qy = track_y[r["query_genome"]]
        ry = track_y[ref_genome]
        colour = SYNTENY_BLOCK_COLOURS.get(r["block_type"], "#999999")
        seen_types.add(r["block_type"])
        # braid: quadrilateral query-top -> reference-top
        poly = [
            (r["query_start"], qy),
            (r["query_end"], qy),
            (r["ref_end"], ry + bar_h),
            (r["ref_start"], ry + bar_h),
        ]
        ax.fill(*zip(*poly, strict=True), color=colour, alpha=0.45, linewidth=0)

    # chromosome bars per track (span of that genome's blocks)
    for genome, y in track_y.items():
        if genome == ref_genome:
            starts = [r["ref_start"] for r in rows]
            ends = [r["ref_end"] for r in rows]
        else:
            gr = [r for r in rows if r["query_genome"] == genome]
            starts = [r["query_start"] for r in gr]
            ends = [r["query_end"] for r in gr]
        if not starts:
            continue
        lo, hi = min(starts), max(ends)
        ax.add_patch(Rectangle((lo, y), hi - lo, bar_h, color="#333333", zorder=3))
        ax.text(lo, y + bar_h + 0.06, genome, fontsize=8, color="#333333", va="bottom")

    ax.set_ylim(-0.4, len(query_genomes) + 0.6)
    ax.set_xlabel(f"Position (reference {ref_genome})")
    ax.set_yticks([])
    ax.set_title("Riparian synteny")
    legend = [
        Patch(facecolor=SYNTENY_BLOCK_COLOURS[t], alpha=0.6, label=t)
        for t in SYNTENY_BLOCK_ORDER if t in seen_types
    ]
    if legend:
        ax.legend(handles=legend, loc="upper right", title="block type")
    return _save(fig, plt, outdir, "riparian", output_format, dpi)


def plot_dotplot(
    block_rows: list[dict[str, Any]],
    outdir: Path,
    *,
    width: float = 7.0,
    height: float = 7.0,
    dpi: int = 150,
    output_format: str = "png",
) -> Path:
    """Render a query-vs-reference dotplot: one line segment per block, typed by colour."""
    import matplotlib  # noqa: PLC0415
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # noqa: PLC0415
    from matplotlib.patches import Patch  # noqa: PLC0415

    from privy.plot.themes import (  # noqa: PLC0415
        SYNTENY_BLOCK_COLOURS,
        SYNTENY_BLOCK_ORDER,
        apply_privy_theme,
    )

    apply_privy_theme()
    outdir.mkdir(parents=True, exist_ok=True)
    rows = _block_rows_to_floats(block_rows)

    fig, ax = plt.subplots(figsize=(width, height))
    if not rows:
        ax.text(0.5, 0.5, "No synteny blocks", ha="center", va="center",
                transform=ax.transAxes, fontsize=12, color="#888888")
        return _save(fig, plt, outdir, "dotplot", output_format, dpi)

    seen_types: set[str] = set()
    for r in rows:
        colour = SYNTENY_BLOCK_COLOURS.get(r["block_type"], "#999999")
        seen_types.add(r["block_type"])
        # forward block -> ascending diagonal; inversion -> descending
        if r["block_type"] == "inversion":
            xs = [r["query_start"], r["query_end"]]
            ys = [r["ref_end"], r["ref_start"]]
        else:
            xs = [r["query_start"], r["query_end"]]
            ys = [r["ref_start"], r["ref_end"]]
        ax.plot(xs, ys, color=colour, linewidth=2.2, alpha=0.85, solid_capstyle="round")

    ax.set_xlabel(f"Query ({rows[0]['query_genome']} …)")
    ax.set_ylabel(f"Reference ({rows[0]['ref_genome']})")
    ax.set_title("Synteny dotplot")
    legend = [
        Patch(facecolor=SYNTENY_BLOCK_COLOURS[t], label=t)
        for t in SYNTENY_BLOCK_ORDER if t in seen_types
    ]
    if legend:
        ax.legend(handles=legend, loc="best", title="block type")
    return _save(fig, plt, outdir, "dotplot", output_format, dpi)


def plot_block_density(
    block_rows: list[dict[str, Any]],
    outdir: Path,
    *,
    window: int | None = None,
    step: int | None = None,
    width: float = 11.0,
    height: float = 3.6,
    dpi: int = 150,
    output_format: str = "png",
) -> Path:
    """Stacked-area density of synteny block types along the reference.

    Windows the primary reference contig and renders the per-window proportion of
    each block type (collinear/inversion/translocation/duplication) + uncovered
    ("missing"), via hierarchical base assignment — showing where rearrangements
    and private structure concentrate.
    """
    import matplotlib  # noqa: PLC0415
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # noqa: PLC0415

    from privy.core.intervals import (  # noqa: PLC0415
        class_proportions,
        hierarchical_base_assignment,
    )
    from privy.plot.themes import (  # noqa: PLC0415
        SYNTENY_BLOCK_COLOURS,
        SYNTENY_BLOCK_ORDER,
        apply_privy_theme,
    )

    apply_privy_theme()
    outdir.mkdir(parents=True, exist_ok=True)
    rows = _block_rows_to_floats(block_rows)

    fig, ax = plt.subplots(figsize=(width, height))
    if not rows:
        ax.text(0.5, 0.5, "No synteny blocks", ha="center", va="center",
                transform=ax.transAxes, fontsize=12, color="#888888")
        return _save(fig, plt, outdir, "block_density", output_format, dpi)

    # Primary reference contig (the one with the most blocks).
    contig = Counter(r["ref_contig"] for r in rows).most_common(1)[0][0]
    crows = [r for r in rows if r["ref_contig"] == contig]
    lo = min(r["ref_start"] for r in crows)
    hi = max(r["ref_end"] for r in crows)
    span = max(1, hi - lo)
    w = window or max(1, span // 40)
    st = step or w

    classes = [*SYNTENY_BLOCK_ORDER, "missing"]
    feats = [(r["ref_start"], r["ref_end"], r["block_type"]) for r in crows]
    xs: list[float] = []
    series: dict[str, list[float]] = {c: [] for c in classes}
    pos = lo
    while pos < hi:
        seg_end = min(pos + w, hi)
        props = class_proportions(
            hierarchical_base_assignment(pos, seg_end, feats, SYNTENY_BLOCK_ORDER)
        )
        xs.append((pos + seg_end) / 2)
        for c in classes:
            series[c].append(props.get(c, 0.0))
        pos += st

    colours = [
        "#eeeeee" if c == "missing" else SYNTENY_BLOCK_COLOURS.get(c, "#cccccc")
        for c in classes
    ]
    ax.stackplot(xs, [series[c] for c in classes], labels=classes, colors=colours, alpha=0.85)
    ax.set_xlabel(f"Reference {contig} position")
    ax.set_ylabel("block-type proportion")
    ax.set_ylim(0, 1)
    ax.set_title("Synteny block-type density")
    ax.legend(loc="upper right", fontsize=7, ncol=2)
    return _save(fig, plt, outdir, "block_density", output_format, dpi)


def _save(fig: Any, plt: Any, outdir: Path, name: str, output_format: str, dpi: int) -> Path:
    outpath = outdir / f"{name}.{output_format}"
    fig.savefig(outpath, dpi=dpi)
    plt.close(fig)
    return outpath
