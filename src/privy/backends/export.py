"""Export privy scan results to downstream genome-tool formats.

The export backend reads existing TSV outputs from ``privy scan`` and writes
format-specific files without re-opening the original VCF or GFA.  The first
supported formats are BED and GFF3, which are useful for IGV, bedtools,
annotation joins, and quick interval inspection.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Literal
from urllib.parse import quote

from privy.io.tsv import read_tsv

ExportKind = Literal["hits", "regions", "both"]
ExportFormat = Literal["bed", "gff3"]


def run_export(
    hits_path: Path | None,
    regions_path: Path | None,
    outdir: Path,
    export_format: ExportFormat = "bed",
    export_kind: ExportKind = "both",
    track_name: str = "Panex Privus",
    include_header: bool = True,
) -> list[Path]:
    """Export scan outputs and return paths written.

    Args:
        hits_path: Optional path to ``hits.tsv``.
        regions_path: Optional path to ``regions.tsv``.
        outdir: Output directory for exported files.
        export_format: Output format: ``"bed"`` or ``"gff3"``.
        export_kind: Which inputs to export: ``"hits"``, ``"regions"``, or ``"both"``.
        track_name: BED track name used when ``include_header`` is true.
        include_header: Write a UCSC-style BED track line.

    Raises:
        ValueError: If inputs are missing for the requested export kind or
            if an unsupported format/kind is requested.
    """
    if export_format not in {"bed", "gff3"}:
        raise ValueError(
            f"Unsupported export format: {export_format!r}. Use 'bed' or 'gff3'."
        )
    if export_kind not in {"hits", "regions", "both"}:
        raise ValueError(
            f"Unsupported export kind: {export_kind!r}. Use 'hits', 'regions', or 'both'."
        )

    outdir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []

    if export_kind in {"hits", "both"}:
        if hits_path is None:
            raise ValueError("--hits is required when exporting hits.")
        hit_rows = read_tsv(hits_path)
        path = outdir / f"hits.{export_format}"
        if export_format == "bed":
            write_hits_bed(
                hit_rows, path, track_name=track_name, include_header=include_header
            )
        else:
            write_hits_gff3(hit_rows, path, include_header=include_header)
        written.append(path)

    if export_kind in {"regions", "both"}:
        if regions_path is None:
            raise ValueError("--regions is required when exporting regions.")
        region_rows = read_tsv(regions_path)
        path = outdir / f"regions.{export_format}"
        if export_format == "bed":
            write_regions_bed(
                region_rows,
                path,
                track_name=f"{track_name} regions",
                include_header=include_header,
            )
        else:
            write_regions_gff3(region_rows, path, include_header=include_header)
        written.append(path)

    meta_path = outdir / "export.json"
    meta_path.write_text(
        json.dumps(
            {
                "tool": "privy export",
                "format": export_format,
                "kind": export_kind,
                "hits_path": str(hits_path) if hits_path is not None else None,
                "regions_path": str(regions_path) if regions_path is not None else None,
                "outputs": [str(p) for p in written],
                "timestamp": datetime.now(timezone.utc).isoformat(),
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    written.append(meta_path)
    return written


def write_hits_bed(
    rows: list[dict[str, str]],
    path: Path,
    track_name: str = "Panex Privus hits",
    include_header: bool = True,
) -> None:
    """Write hit rows to BED6 + score/detail columns."""
    with open(path, "w", encoding="utf-8") as fh:
        if include_header:
            fh.write(f'track name="{track_name}" description="Panex Privus private loci"\n')
        for row in rows:
            contig, start, end = _parse_interval(row, id_field="locus_id")
            name = row.get("locus_id", "hit")
            score = _bed_score(row.get("final_score", "0"))
            fields = [
                contig,
                str(start),
                str(end),
                name,
                str(score),
                ".",
                row.get("strictness_class", "NA"),
                row.get("variant_type", "NA"),
                row.get("allele_key", "NA"),
                row.get("final_score", "NA"),
            ]
            fh.write("\t".join(fields) + "\n")


def write_regions_bed(
    rows: list[dict[str, str]],
    path: Path,
    track_name: str = "Panex Privus regions",
    include_header: bool = True,
) -> None:
    """Write region rows to BED6 + summary columns."""
    with open(path, "w", encoding="utf-8") as fh:
        if include_header:
            fh.write(f'track name="{track_name}" description="Panex Privus private regions"\n')
        for row in rows:
            contig, start, end = _parse_interval(row, id_field="region_id")
            name = row.get("region_id", "region")
            score = _bed_score(row.get("final_score", "0"))
            fields = [
                contig,
                str(start),
                str(end),
                name,
                str(score),
                ".",
                row.get("dominant_strictness_class", "NA"),
                row.get("variant_types", "NA"),
                row.get("n_loci", "NA"),
                row.get("final_score", "NA"),
            ]
            fh.write("\t".join(fields) + "\n")


def write_hits_gff3(
    rows: list[dict[str, str]],
    path: Path,
    include_header: bool = True,
) -> None:
    """Write hit rows to GFF3 feature records."""
    with open(path, "w", encoding="utf-8") as fh:
        if include_header:
            fh.write("##gff-version 3\n")
        for row in rows:
            contig, start, end = _parse_interval(row, id_field="locus_id")
            gff_start, gff_end = _gff3_interval(start, end)
            locus_id = row.get("locus_id", "hit")
            attrs = _gff3_attrs({
                "ID": locus_id,
                "Name": locus_id,
                "strictness_class": row.get("strictness_class", "NA"),
                "variant_type": row.get("variant_type", "NA"),
                "allele_key": row.get("allele_key", "NA"),
                "final_score": row.get("final_score", "NA"),
            })
            fields = [
                contig,
                "privy",
                "sequence_variant",
                str(gff_start),
                str(gff_end),
                _gff3_score(row.get("final_score", "")),
                ".",
                ".",
                attrs,
            ]
            fh.write("\t".join(fields) + "\n")


def write_regions_gff3(
    rows: list[dict[str, str]],
    path: Path,
    include_header: bool = True,
) -> None:
    """Write region rows to GFF3 feature records."""
    with open(path, "w", encoding="utf-8") as fh:
        if include_header:
            fh.write("##gff-version 3\n")
        for row in rows:
            contig, start, end = _parse_interval(row, id_field="region_id")
            gff_start, gff_end = _gff3_interval(start, end)
            region_id = row.get("region_id", "region")
            attrs = _gff3_attrs({
                "ID": region_id,
                "Name": region_id,
                "dominant_strictness_class": row.get("dominant_strictness_class", "NA"),
                "variant_types": row.get("variant_types", "NA"),
                "n_loci": row.get("n_loci", "NA"),
                "final_score": row.get("final_score", "NA"),
            })
            fields = [
                contig,
                "privy",
                "region",
                str(gff_start),
                str(gff_end),
                _gff3_score(row.get("final_score", "")),
                ".",
                ".",
                attrs,
            ]
            fh.write("\t".join(fields) + "\n")


def _parse_interval(row: dict[str, str], id_field: str) -> tuple[str, int, int]:
    """Return ``(contig, start, end)`` from a TSV row."""
    try:
        contig = row["contig"]
        start = int(row["start"])
        end = int(row["end"])
    except (KeyError, ValueError) as exc:
        row_id = row.get(id_field, "<unknown>")
        raise ValueError(f"Malformed interval row for {row_id}: {row!r}") from exc

    if start < 0:
        row_id = row.get(id_field, "<unknown>")
        raise ValueError(f"BED start must be non-negative for {row_id}: {start}")
    if end < start:
        row_id = row.get(id_field, "<unknown>")
        raise ValueError(f"BED end is before start for {row_id}: {start}-{end}")
    return contig, start, end


def _bed_score(raw_score: str) -> int:
    """Convert a floating score to the BED integer score range [0, 1000]."""
    try:
        value = float(raw_score)
    except ValueError:
        return 0
    scaled = round(value * 1000)
    return max(0, min(1000, scaled))


def _gff3_interval(start: int, end: int) -> tuple[int, int]:
    """Convert a 0-based half-open interval to 1-based closed GFF3 coordinates."""
    return start + 1, max(end, start + 1)


def _gff3_score(raw_score: str) -> str:
    """Return a valid GFF3 score field."""
    try:
        value = float(raw_score)
    except ValueError:
        return "."
    return f"{value:.6g}"


def _gff3_attrs(values: dict[str, str]) -> str:
    """Format a GFF3 attribute field with URL-escaped values."""
    attrs = []
    for key, value in values.items():
        if value == "":
            continue
        attrs.append(f"{key}={quote(value, safe='._:-')}")
    return ";".join(attrs)
