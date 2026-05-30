"""Build the interactive synteny dashboard (Model C: data-injected single-file HTML).

Reads the artifacts written by ``privy synteny`` (``synteny_blocks.tsv``,
``synteny_regions.tsv``, ``synteny.json``), assembles a compact JSON data model,
and injects it into the prebuilt single-file dashboard bundle shipped as a package
asset.  The result is one self-contained ``.html`` — no server, no Node, fully
offline — that renders linked riparian + dotplot views with target-private
highlighting (the front end is built from ``web/`` by developers; users never need
Node).  See scratch/notes/51_frontend_stack_modelC.md.
"""

from __future__ import annotations

import csv
import json
import logging
from pathlib import Path
from typing import Any

log = logging.getLogger("privy.interactive.synteny_dashboard")

_ASSET = Path(__file__).parent / "assets" / "synteny_dashboard.html"
_PLACEHOLDER = '{"__privy_placeholder__": true}'


def build_synteny_dashboard(
    input_dir: Path,
    *,
    outdir: Path | None = None,
    output_name: str = "synteny_dashboard.html",
) -> Path:
    """Build a self-contained interactive synteny dashboard from a synteny directory.

    Args:
        input_dir: A ``privy synteny`` output directory (must contain
            ``synteny_blocks.tsv``).
        outdir: Where to write the dashboard (default: *input_dir*).
        output_name: Output file name.

    Returns:
        Path to the written ``.html``.

    Raises:
        FileNotFoundError: If ``synteny_blocks.tsv`` or the bundle asset is missing.
        RuntimeError: If the asset's data placeholder cannot be found.
    """
    input_dir = Path(input_dir)
    blocks_tsv = input_dir / "synteny_blocks.tsv"
    if not blocks_tsv.exists():
        raise FileNotFoundError(f"synteny_blocks.tsv not found in {input_dir}")
    if not _ASSET.exists():
        raise FileNotFoundError(
            f"dashboard bundle asset missing: {_ASSET} (rebuild from web/ with `npm run build`)"
        )

    data = _build_data_model(input_dir)
    template = _ASSET.read_text(encoding="utf-8")
    if _PLACEHOLDER not in template:
        raise RuntimeError(
            "dashboard asset has no data placeholder; rebuild web/ bundle from current source."
        )
    # Escape '</' so the JSON cannot break out of the <script> element.
    payload = json.dumps(data, separators=(",", ":")).replace("</", "<\\/")
    rendered = template.replace(_PLACEHOLDER, payload, 1)

    dest_dir = Path(outdir) if outdir is not None else input_dir
    dest_dir.mkdir(parents=True, exist_ok=True)
    out_path = dest_dir / output_name
    out_path.write_text(rendered, encoding="utf-8")
    log.info(
        "Synteny dashboard written | blocks=%d regions=%d -> %s",
        len(data["blocks"]), len(data["regions"]), out_path,
    )
    return out_path


def _build_data_model(input_dir: Path) -> dict[str, Any]:
    """Assemble the JSON data model the dashboard front end consumes."""
    blocks = _read_blocks(input_dir / "synteny_blocks.tsv")
    regions = _read_regions(input_dir / "synteny_regions.tsv")
    meta_path = input_dir / "synteny.json"
    file_meta = json.loads(meta_path.read_text(encoding="utf-8")) if meta_path.exists() else {}

    reference = blocks[0]["ref_genome"] if blocks else file_meta.get("reference", "")
    genomes: list[str] = []
    if reference:
        genomes.append(reference)
    for b in blocks:
        if b["query_genome"] not in genomes:
            genomes.append(b["query_genome"])

    n_private = sum(1 for r in regions if r["target_private"] is True)
    meta = {
        "n_blocks": len(blocks),
        "n_regions": len(regions),
        "n_target_private": n_private,
        "block_type_counts": file_meta.get("block_type_counts", {}),
    }
    return {"reference": reference, "genomes": genomes, "blocks": blocks,
            "regions": regions, "meta": meta}


def _read_blocks(path: Path) -> list[dict[str, Any]]:
    with open(path, encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return [
        {
            "block_id": r["block_id"],
            "query_genome": r["query_genome"],
            "query_start": int(r["query_start"]),
            "query_end": int(r["query_end"]),
            "ref_genome": r["ref_genome"],
            "ref_contig": r["ref_contig"],
            "ref_start": int(r["ref_start"]),
            "ref_end": int(r["ref_end"]),
            "strand": r["strand"],
            "block_type": r["block_type"],
        }
        for r in rows
    ]


def _read_regions(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    with open(path, encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return [
        {
            "region_id": r["region_id"],
            "ref_contig": r["ref_contig"],
            "ref_start": int(r["ref_start"]),
            "ref_end": int(r["ref_end"]),
            "target_private": str(r.get("target_private", "")).lower() == "true",
        }
        for r in rows
    ]
