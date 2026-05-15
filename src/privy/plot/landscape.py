"""Landscape analysis plots."""

from __future__ import annotations

import logging
import re
from collections import defaultdict
from pathlib import Path
from typing import Any

log = logging.getLogger("privy.plot.landscape")


def plot_landscape_heatmap(
    sample_rows: list[dict[str, Any]],
    outdir: Path,
    value_column: str,
    filename_stem: str,
    title: str,
    cmap: str = "viridis",
    width: float = 10.0,
    height: float = 5.5,
    dpi: int = 150,
    output_format: str = "png",
) -> Path:
    """Plot a numeric sample-by-window landscape heatmap."""
    import matplotlib  # noqa: PLC0415
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # noqa: PLC0415
    import numpy as np  # noqa: PLC0415

    from privy.plot.themes import apply_privy_theme  # noqa: PLC0415

    apply_privy_theme()
    outdir.mkdir(parents=True, exist_ok=True)
    samples = _sample_order(sample_rows)
    windows = _window_order_from_sample_rows(sample_rows)

    fig_width = max(width, min(18.0, 2.8 + len(windows) * 0.28))
    fig_height = max(height, min(16.0, 2.2 + len(samples) * 0.28))
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    if not samples or not windows:
        ax.text(0.5, 0.5, "No landscape rows", ha="center", va="center",
                transform=ax.transAxes, fontsize=12, color="#888888")
    else:
        sample_index = {sample: i for i, sample in enumerate(samples)}
        window_index = {window: i for i, window in enumerate(windows)}
        matrix = np.full((len(samples), len(windows)), np.nan)
        for row in sample_rows:
            value = _to_optional_float(row.get(value_column))
            if value is None:
                continue
            matrix[sample_index[str(row["sample"])], window_index[str(row["window_id"])]] = value
        image = ax.imshow(matrix, aspect="auto", interpolation="nearest", cmap=cmap)
        fig.colorbar(image, ax=ax, fraction=0.025, pad=0.02, label=value_column)
        _style_window_axis(ax, windows, sample_rows)
        ax.set_yticks(range(len(samples)))
        ax.set_yticklabels(samples)

    ax.set_title(title)
    ax.set_xlabel("Window")
    ax.set_ylabel("Sample")
    outpath = outdir / f"{filename_stem}.{output_format}"
    fig.savefig(outpath, dpi=dpi)
    plt.close(fig)
    return outpath


def plot_local_background_map(
    sample_rows: list[dict[str, Any]],
    outdir: Path,
    filename_stem: str = "local_background_map",
    title: str = "Local background map",
    width: float = 10.0,
    height: float = 5.5,
    dpi: int = 150,
    output_format: str = "png",
) -> Path:
    """Plot sample-by-window nearest-background assignments."""
    import matplotlib  # noqa: PLC0415
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # noqa: PLC0415
    import numpy as np  # noqa: PLC0415
    from matplotlib.colors import ListedColormap  # noqa: PLC0415
    from matplotlib.lines import Line2D  # noqa: PLC0415

    from privy.plot.themes import apply_privy_theme  # noqa: PLC0415

    apply_privy_theme()
    outdir.mkdir(parents=True, exist_ok=True)
    samples = _sample_order(sample_rows)
    windows = _window_order_from_sample_rows(sample_rows)
    backgrounds = ["unassigned", *samples]
    background_to_code = {background: i for i, background in enumerate(backgrounds)}

    colors = ["#d9d9d9", *_palette(len(samples))]
    cmap = ListedColormap(colors)
    fig_width = max(width, min(18.0, 2.8 + len(windows) * 0.28))
    fig_height = max(height, min(16.0, 2.2 + len(samples) * 0.28))
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    if not samples or not windows:
        ax.text(0.5, 0.5, "No background assignments", ha="center", va="center",
                transform=ax.transAxes, fontsize=12, color="#888888")
    else:
        sample_index = {sample: i for i, sample in enumerate(samples)}
        window_index = {window: i for i, window in enumerate(windows)}
        matrix = np.zeros((len(samples), len(windows)), dtype=int)
        for row in sample_rows:
            nearest = str(row.get("nearest_background", "NA"))
            if nearest == "NA" or nearest not in background_to_code:
                nearest = "unassigned"
            matrix[sample_index[str(row["sample"])], window_index[str(row["window_id"])]] = (
                background_to_code[nearest]
            )
        ax.imshow(matrix, aspect="auto", interpolation="nearest", cmap=cmap,
                  vmin=0, vmax=len(backgrounds) - 1)
        _style_window_axis(ax, windows, sample_rows)
        ax.set_yticks(range(len(samples)))
        ax.set_yticklabels(samples)
        if len(backgrounds) <= 21:
            handles = [
                Line2D([0], [0], marker="s", linestyle="", markersize=7,
                       markerfacecolor=colors[i], markeredgewidth=0, label=background)
                for i, background in enumerate(backgrounds)
            ]
            ax.legend(handles=handles, title="Nearest background", loc="upper left",
                      bbox_to_anchor=(1.01, 1.0), borderaxespad=0)

    ax.set_title(title)
    ax.set_xlabel("Window")
    ax.set_ylabel("Sample")
    outpath = outdir / f"{filename_stem}.{output_format}"
    fig.savefig(outpath, dpi=dpi)
    plt.close(fig)
    return outpath


def plot_similarity_cluster_map(
    similarity_rows: list[dict[str, Any]],
    outdir: Path,
    filename_stem: str = "similarity_cluster_map",
    title: str = "Local similarity cluster map",
    width: float = 7.0,
    height: float = 6.0,
    dpi: int = 150,
    output_format: str = "png",
) -> Path:
    """Plot an average local genotype-similarity heatmap with clustered sample order."""
    import matplotlib  # noqa: PLC0415
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # noqa: PLC0415
    import numpy as np  # noqa: PLC0415

    from privy.plot.themes import apply_privy_theme  # noqa: PLC0415

    apply_privy_theme()
    outdir.mkdir(parents=True, exist_ok=True)
    samples = _sample_order_from_similarity(similarity_rows)
    matrix = _mean_similarity_matrix(samples, similarity_rows)
    order = _cluster_order(matrix) if len(samples) > 2 else list(range(len(samples)))
    ordered_samples = [samples[i] for i in order]
    ordered_matrix = matrix[np.ix_(order, order)] if len(samples) else matrix

    fig_size = max(width, min(12.0, 2.5 + len(samples) * 0.45))
    fig, ax = plt.subplots(figsize=(fig_size, max(height, fig_size * 0.85)))
    if not samples:
        ax.text(0.5, 0.5, "No similarity rows", ha="center", va="center",
                transform=ax.transAxes, fontsize=12, color="#888888")
    else:
        image = ax.imshow(ordered_matrix, aspect="equal", interpolation="nearest",
                          cmap="magma", vmin=0, vmax=1)
        fig.colorbar(image, ax=ax, fraction=0.035, pad=0.02, label="mean similarity")
        ax.set_xticks(range(len(ordered_samples)))
        ax.set_xticklabels(ordered_samples, rotation=45, ha="right")
        ax.set_yticks(range(len(ordered_samples)))
        ax.set_yticklabels(ordered_samples)

    ax.set_title(title)
    outpath = outdir / f"{filename_stem}.{output_format}"
    fig.savefig(outpath, dpi=dpi)
    plt.close(fig)
    return outpath


def plot_all_landscape(
    sample_rows: list[dict[str, Any]],
    window_rows: list[dict[str, Any]],
    similarity_rows: list[dict[str, Any]],
    outdir: Path,
    output_format: str = "png",
) -> list[Path]:
    """Generate all first-pass landscape plots."""
    paths, _index_rows = plot_landscape_set(
        sample_rows=sample_rows,
        window_rows=window_rows,
        similarity_rows=similarity_rows,
        outdir=outdir,
        output_format=output_format,
        plot_scope="genome",
        contigs=None,
    )
    return paths


def plot_landscape_set(
    sample_rows: list[dict[str, Any]],
    window_rows: list[dict[str, Any]],
    similarity_rows: list[dict[str, Any]],
    outdir: Path,
    output_format: str = "png",
    plot_scope: str = "chromosome",
    contigs: list[str] | None = None,
) -> tuple[list[Path], list[dict[str, Any]]]:
    """Generate landscape plots and an index describing each figure."""
    if plot_scope not in {"chromosome", "genome", "both"}:
        raise ValueError("--plot-scope must be one of: chromosome, genome, both.")
    if plot_scope == "genome" and contigs:
        raise ValueError("--contig and --contigs require chromosome or both plot scope.")

    paths: list[Path] = []
    index_rows: list[dict[str, Any]] = []
    selected_contigs = _selected_contigs(sample_rows, contigs)
    window_stats = _window_stats(window_rows)

    def append_plot(
        path: Path,
        plot_type: str,
        scope: str,
        contig: str,
    ) -> None:
        paths.append(path)
        stats = (
            window_stats.get(contig, {"n_windows": 0, "start": 0, "end": 0})
            if scope == "chromosome"
            else _genome_stats(window_stats)
        )
        index_rows.append({
            "plot_type": plot_type,
            "plot_scope": scope,
            "contig": contig if scope == "chromosome" else "genome",
            "path": path.name,
            "n_windows": stats["n_windows"],
            "start": stats["start"],
            "end": stats["end"],
            "output_format": output_format,
        })

    if plot_scope in {"chromosome", "both"}:
        for contig in selected_contigs:
            safe_contig = _safe_filename_component(contig)
            contig_sample_rows = [
                row for row in sample_rows if str(row.get("contig")) == contig
            ]
            contig_similarity_rows = [
                row for row in similarity_rows if str(row.get("contig")) == contig
            ]
            if not contig_sample_rows:
                continue

            log.info("Rendering missingness heatmap | contig=%s", contig)
            append_plot(
                plot_landscape_heatmap(
                    contig_sample_rows,
                    outdir,
                    value_column="missing_rate",
                    filename_stem=f"missingness_heatmap.{safe_contig}",
                    title=f"Windowed missingness: {contig}",
                    cmap="mako" if _has_colormap("mako") else "viridis",
                    output_format=output_format,
                ),
                "missingness_heatmap",
                "chromosome",
                contig,
            )
            log.info("Rendering private ALT burden heatmap | contig=%s", contig)
            append_plot(
                plot_landscape_heatmap(
                    contig_sample_rows,
                    outdir,
                    value_column="private_alt_rate",
                    filename_stem=f"private_burden_heatmap.{safe_contig}",
                    title=f"Windowed private ALT burden: {contig}",
                    cmap="rocket" if _has_colormap("rocket") else "plasma",
                    output_format=output_format,
                ),
                "private_burden_heatmap",
                "chromosome",
                contig,
            )
            log.info("Rendering local background map | contig=%s", contig)
            append_plot(
                plot_local_background_map(
                    contig_sample_rows,
                    outdir,
                    filename_stem=f"local_background_map.{safe_contig}",
                    title=f"Local background map: {contig}",
                    output_format=output_format,
                ),
                "local_background_map",
                "chromosome",
                contig,
            )
            if contig_similarity_rows:
                log.info("Rendering similarity cluster map | contig=%s", contig)
                append_plot(
                    plot_similarity_cluster_map(
                        contig_similarity_rows,
                        outdir,
                        filename_stem=f"similarity_cluster_map.{safe_contig}",
                        title=f"Local similarity cluster map: {contig}",
                        output_format=output_format,
                    ),
                    "similarity_cluster_map",
                    "chromosome",
                    contig,
                )

    if plot_scope in {"genome", "both"}:
        log.info("Rendering missingness heatmap | scope=genome")
        append_plot(
            plot_landscape_heatmap(
                sample_rows,
                outdir,
                value_column="missing_rate",
                filename_stem="missingness_heatmap",
                title="Windowed missingness",
                cmap="mako" if _has_colormap("mako") else "viridis",
                output_format=output_format,
            ),
            "missingness_heatmap",
            "genome",
            "genome",
        )
        log.info("Rendering private ALT burden heatmap | scope=genome")
        append_plot(
            plot_landscape_heatmap(
                sample_rows,
                outdir,
                value_column="private_alt_rate",
                filename_stem="private_burden_heatmap",
                title="Windowed private ALT burden",
                cmap="rocket" if _has_colormap("rocket") else "plasma",
                output_format=output_format,
            ),
            "private_burden_heatmap",
            "genome",
            "genome",
        )
        log.info("Rendering local background map | scope=genome")
        append_plot(
            plot_local_background_map(sample_rows, outdir, output_format=output_format),
            "local_background_map",
            "genome",
            "genome",
        )
        log.info("Rendering similarity cluster map | scope=genome")
        append_plot(
            plot_similarity_cluster_map(similarity_rows, outdir, output_format=output_format),
            "similarity_cluster_map",
            "genome",
            "genome",
        )
    return paths, index_rows


def _selected_contigs(
    sample_rows: list[dict[str, Any]],
    requested: list[str] | None,
) -> list[str]:
    available = _ordered_contigs(sample_rows)
    if requested is None or not requested:
        return available
    available_set = set(available)
    missing = [contig for contig in requested if contig not in available_set]
    if missing:
        raise ValueError(
            "Requested contig(s) not found in landscape output: "
            + ", ".join(missing)
        )
    requested_set = set(requested)
    return [contig for contig in available if contig in requested_set]


def _ordered_contigs(rows: list[dict[str, Any]]) -> list[str]:
    contigs: list[str] = []
    seen: set[str] = set()
    for row in rows:
        contig = str(row["contig"])
        if contig not in seen:
            seen.add(contig)
            contigs.append(contig)
    return contigs


def _window_stats(rows: list[dict[str, Any]]) -> dict[str, dict[str, int]]:
    stats: dict[str, dict[str, int]] = {}
    seen_windows: set[tuple[str, str]] = set()
    for row in rows:
        contig = str(row["contig"])
        window = str(row["window_id"])
        if (contig, window) in seen_windows:
            continue
        seen_windows.add((contig, window))
        start = int(row["start"])
        end = int(row["end"])
        current = stats.setdefault(
            contig,
            {"n_windows": 0, "start": start, "end": end},
        )
        current["n_windows"] += 1
        current["start"] = min(current["start"], start)
        current["end"] = max(current["end"], end)
    return stats


def _genome_stats(stats_by_contig: dict[str, dict[str, int]]) -> dict[str, int]:
    if not stats_by_contig:
        return {"n_windows": 0, "start": 0, "end": 0}
    return {
        "n_windows": sum(stats["n_windows"] for stats in stats_by_contig.values()),
        "start": min(stats["start"] for stats in stats_by_contig.values()),
        "end": max(stats["end"] for stats in stats_by_contig.values()),
    }


def _safe_filename_component(value: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("._")
    return safe or "contig"


def _sample_order(rows: list[dict[str, Any]]) -> list[str]:
    samples: list[str] = []
    seen: set[str] = set()
    for row in rows:
        sample = str(row["sample"])
        if sample not in seen:
            seen.add(sample)
            samples.append(sample)
    return samples


def _window_order_from_sample_rows(rows: list[dict[str, Any]]) -> list[str]:
    windows: list[str] = []
    seen: set[str] = set()
    for row in rows:
        window = str(row["window_id"])
        if window not in seen:
            seen.add(window)
            windows.append(window)
    return windows


def _sample_order_from_similarity(rows: list[dict[str, Any]]) -> list[str]:
    samples: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for key in ("sample_a", "sample_b"):
            sample = str(row[key])
            if sample not in seen:
                seen.add(sample)
                samples.append(sample)
    return samples


def _style_window_axis(ax: Any, windows: list[str], sample_rows: list[dict[str, Any]]) -> None:
    if not windows:
        return
    contigs: dict[str, str] = {}
    starts: dict[str, int] = {}
    for row in sample_rows:
        window = str(row["window_id"])
        contigs.setdefault(window, str(row["contig"]))
        starts.setdefault(window, int(row["start"]))
        if len(contigs) == len(windows):
            break
    tick_positions: list[int] = []
    tick_labels: list[str] = []
    last_contig: str | None = None
    for i, window in enumerate(windows):
        contig = contigs.get(window, "")
        if contig != last_contig:
            tick_positions.append(i)
            tick_labels.append(contig)
            if i > 0:
                ax.axvline(i - 0.5, color="#ffffff", linewidth=1.2)
            last_contig = contig
    if len(tick_positions) == 1 and len(windows) > 1:
        step = max(1, len(windows) // 6)
        tick_positions = list(range(0, len(windows), step))
        tick_labels = [f"{starts.get(windows[i], 0) // 1000} kb" for i in tick_positions]
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels, rotation=45, ha="right")


def _mean_similarity_matrix(samples: list[str], rows: list[dict[str, Any]]) -> Any:
    import numpy as np  # noqa: PLC0415

    sample_to_index = {sample: i for i, sample in enumerate(samples)}
    matrix = np.eye(len(samples))
    totals: dict[tuple[int, int], float] = defaultdict(float)
    counts: dict[tuple[int, int], int] = defaultdict(int)
    for row in rows:
        value = _to_optional_float(row.get("similarity"))
        if value is None:
            continue
        left = sample_to_index[str(row["sample_a"])]
        right = sample_to_index[str(row["sample_b"])]
        pair = (min(left, right), max(left, right))
        totals[pair] += value
        counts[pair] += 1
    for (left, right), total in totals.items():
        value = total / counts[(left, right)]
        matrix[left, right] = value
        matrix[right, left] = value
    return matrix


def _cluster_order(matrix: Any) -> list[int]:
    import numpy as np  # noqa: PLC0415

    n = matrix.shape[0]
    clusters: list[list[int]] = [[i] for i in range(n)]
    while len(clusters) > 1:
        best_pair: tuple[int, int] | None = None
        best_distance = float("inf")
        for i, left in enumerate(clusters):
            for j in range(i + 1, len(clusters)):
                right = clusters[j]
                values = [1 - float(matrix[a, b]) for a in left for b in right]
                distance = float(np.mean(values)) if values else float("inf")
                if distance < best_distance:
                    best_distance = distance
                    best_pair = (i, j)
        if best_pair is None:
            break
        left_i, right_i = best_pair
        merged = [*clusters[left_i], *clusters[right_i]]
        clusters = [
            cluster
            for index, cluster in enumerate(clusters)
            if index not in {left_i, right_i}
        ]
        clusters.append(merged)
    return clusters[0] if clusters else []


def _palette(n: int) -> list[str]:
    base = [
        "#0868ac",
        "#b35806",
        "#542788",
        "#1b9e77",
        "#d95f02",
        "#7570b3",
        "#e7298a",
        "#66a61e",
        "#e6ab02",
        "#a6761d",
        "#1f78b4",
        "#33a02c",
        "#fb9a99",
        "#e31a1c",
        "#fdbf6f",
        "#ff7f00",
        "#cab2d6",
        "#6a3d9a",
        "#b2df8a",
        "#a6cee3",
    ]
    if n <= len(base):
        return base[:n]
    return [base[i % len(base)] for i in range(n)]


def _has_colormap(name: str) -> bool:
    try:
        import matplotlib.pyplot as plt  # noqa: PLC0415

        plt.get_cmap(name)
    except ValueError:
        return False
    return True


def _to_optional_float(value: object) -> float | None:
    if value == "NA" or value is None:
        return None
    try:
        return float(str(value))
    except (TypeError, ValueError):
        return None
