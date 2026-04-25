"""Annotation classification engine for privy annotate (v0.7).

Reads hits.tsv produced by ``privy scan``, intersects each locus with a GFF3
annotation, and writes annotated_hits.tsv + annotation_summary.tsv.

Annotation class hierarchy (most → least specific):
    CDS → UTR → exonic → intronic → intergenic
"""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

from privy.io.gff import (
    UTR_FEATURES,
    AnnotationIndex,
    build_annotation_index,
    load_contig_alias,
    query_genes,
    query_sub_feature,
)
from privy.io.tsv import (
    ANNOTATED_HITS_COLUMNS,
    ANNOTATION_SUMMARY_COLUMNS,
    TsvWriter,
    read_tsv,
)

# Canonical ordering for annotation classes in summary output
ANNOTATION_ORDER: list[str] = ["CDS", "UTR", "exonic", "intronic", "intergenic"]


# ---------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------

def classify_locus(
    idx: AnnotationIndex,
    contig: str,
    start: int,
    end: int,
) -> tuple[str, str, str, int, int]:
    """Classify a single locus against the annotation index.

    Args:
        idx: Populated :class:`~privy.io.gff.AnnotationIndex`.
        contig: Contig name (must already be normalised to GFF3 namespace).
        start: 0-based inclusive start.
        end: 0-based exclusive end.

    Returns:
        Tuple of ``(annotation_class, gene_id, gene_strand, gene_start, gene_end)``.
        When intergenic, gene fields are empty strings and coords are -1.
    """
    overlapping_genes = query_genes(idx, contig, start, end)

    if not overlapping_genes:
        return ("intergenic", "", "", -1, -1)

    # Use the first (leftmost) overlapping gene for annotation
    g_start, g_end, gene_id, strand = overlapping_genes[0]

    # Walk the hierarchy: CDS → UTR → exonic → intronic
    if query_sub_feature(idx, contig, "CDS", start, end):
        return ("CDS", gene_id, strand, g_start, g_end)

    for utr_type in UTR_FEATURES:
        if query_sub_feature(idx, contig, utr_type, start, end):
            return ("UTR", gene_id, strand, g_start, g_end)

    if query_sub_feature(idx, contig, "exon", start, end):
        return ("exonic", gene_id, strand, g_start, g_end)

    # Within gene body but no exon → intronic
    return ("intronic", gene_id, strand, g_start, g_end)


# ---------------------------------------------------------------------------
# Main runner
# ---------------------------------------------------------------------------

def run_annotate(
    hits_path: Path,
    gff_path: Path,
    outdir: Path,
    contig_alias_path: Path | None = None,
    hits_contig_to_gff: bool = False,
) -> None:
    """Annotate private loci from *hits_path* with GFF3 gene features.

    Reads hits.tsv, builds a GFF3 index, classifies each locus, writes:
      - ``annotated_hits.tsv`` — all hits columns + annotation columns
      - ``annotation_summary.tsv`` — annotation class counts and percentages
      - ``annotate.json`` — run metadata

    Args:
        hits_path: Path to a hits.tsv file produced by ``privy scan``.
        gff_path: Path to a GFF3 annotation file (plain or .gz).
        outdir: Directory where output files are written.
        contig_alias_path: Optional two-column TSV mapping hits contig names
            to GFF3 contig names.  If None, names are used as-is.
        hits_contig_to_gff: If True, apply alias in the hits→GFF direction
            (hits name → GFF name).  If False, alias maps GFF names → hits
            names (GFF name → hits name, inverted at load time).
    """
    outdir.mkdir(parents=True, exist_ok=True)

    # Load alias map (hits contig → gff contig)
    alias: dict[str, str] = {}
    if contig_alias_path is not None:
        raw_alias = load_contig_alias(contig_alias_path)
        if hits_contig_to_gff:
            alias = raw_alias
        else:
            # alias file maps gff→hits; invert so we can look up hits→gff
            alias = {v: k for k, v in raw_alias.items()}

    # Build GFF3 index
    idx = build_annotation_index(gff_path)

    # Read hits
    hits_rows = read_tsv(hits_path)

    annotated: list[dict[str, object]] = []
    for row in hits_rows:
        raw_contig = row["contig"]
        gff_contig = alias.get(raw_contig, raw_contig)

        try:
            start = int(row["start"])
            end = int(row["end"])
        except (KeyError, ValueError):
            start = 0
            end = 1

        ann_class, gene_id, strand, g_start, g_end = classify_locus(
            idx, gff_contig, start, end
        )
        annotated.append({
            **row,
            "annotation_class": ann_class,
            "gene_id": gene_id,
            "gene_strand": strand,
            "gene_start": g_start if g_start >= 0 else "",
            "gene_end": g_end if g_end >= 0 else "",
        })

    # Write annotated_hits.tsv
    annotated_hits_path = outdir / "annotated_hits.tsv"
    with TsvWriter(annotated_hits_path, ANNOTATED_HITS_COLUMNS) as w:
        w.write_rows(annotated)

    # Compute annotation_summary.tsv
    counts: Counter[str] = Counter(str(r["annotation_class"]) for r in annotated)
    total = len(annotated)
    summary_rows = []
    for cls in ANNOTATION_ORDER:
        n = counts.get(cls, 0)
        pct = round(100.0 * n / total, 2) if total > 0 else 0.0
        summary_rows.append({"annotation_class": cls, "n_loci": n, "pct_total": pct})

    annotation_summary_path = outdir / "annotation_summary.tsv"
    with TsvWriter(annotation_summary_path, ANNOTATION_SUMMARY_COLUMNS) as w:
        w.write_rows(summary_rows)

    # Write annotate.json metadata
    meta = {
        "hits_path": str(hits_path),
        "gff_path": str(gff_path),
        "n_hits": total,
        "annotation_counts": dict(counts),
        "contig_alias_path": str(contig_alias_path) if contig_alias_path else None,
    }
    with open(outdir / "annotate.json", "w", encoding="utf-8") as fh:
        json.dump(meta, fh, indent=2)
