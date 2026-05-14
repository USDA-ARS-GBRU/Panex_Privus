"""TSV output writers for Panex Privus.

All output files use tab-separated values with explicit column headers.
Writers are context managers so handles are closed reliably even on error.

Column schemas here are the canonical definitions for all TSV outputs.
They must match the column specs in the architecture documentation and README.
"""

from __future__ import annotations

import csv
from collections.abc import Iterable
from pathlib import Path
from typing import Any

# ---------------------------------------------------------------------------
# Column schemas — canonical, matches docs/architecture.md and README.md
# ---------------------------------------------------------------------------

HITS_COLUMNS: list[str] = [
    "locus_id",
    "contig",
    "start",
    "end",
    "variant_type",
    "allele_key",
    "target_support_n",
    "target_total_n",
    "offtarget_support_n",
    "offtarget_total_n",
    "target_missing_n",
    "offtarget_missing_n",
    "strictness_class",
    "discovery_score",
    "support_score",
    "penalty_score",
    "final_score",
]

REGIONS_COLUMNS: list[str] = [
    "region_id",
    "contig",
    "start",
    "end",
    "n_loci",
    "variant_types",
    "dominant_strictness_class",
    "target_consistency",
    "offtarget_exclusion",
    "final_score",
]

EVIDENCE_COLUMNS: list[str] = [
    "locus_id",
    "source_type",
    "sample_id",
    "evidence_class",
    "metric_name",
    "metric_value",
    "details",
]

SAMPLE_SUPPORT_COLUMNS: list[str] = [
    "locus_id",
    "sample_id",
    "cohort_role",
    "genotype",
    "allele_supported",
    "depth",
    "allele_fraction",
    "evidence_class",
]

GFA_SEGMENT_COLUMNS: list[str] = [
    "locus_id",
    "contig",
    "start",
    "end",
    "segment_name",
    "segment_length",
    "segment_length_class",
    "graph_signal_type",
    "target_traverse_n",
    "target_total_n",
    "target_coordinate_covered_n",
    "target_missing_n",
    "offtarget_same_segment_traverse_n",
    "offtarget_same_segment_absent_n",
    "offtarget_coordinate_covered_n",
    "offtarget_missing_n",
    "offtarget_total_n",
    "strictness_class",
    "interpretation",
]

QC_COLUMNS: list[str] = [
    "metric",
    "value",
    "description",
]

COMPARE_COLUMNS: list[str] = [
    "compare_id",
    "locus_id_a",
    "locus_id_b",
    "source_a",
    "source_b",
    "contig",
    "start_a",
    "end_a",
    "start_b",
    "end_b",
    "match_class",
    "coordinate_overlap",
    "state_compatibility",
    "strictness_a",
    "strictness_b",
    "support_summary",
    "contradiction_summary",
    "comparison_score",
]

COMPARE_SUMMARY_COLUMNS: list[str] = [
    "match_class",
    "n_loci",
    "pct_total",
    "mean_overlap",
    "mean_score",
]

ANNOTATED_HITS_COLUMNS: list[str] = [
    *HITS_COLUMNS,
    "annotation_class",
    "gene_id",
    "gene_strand",
    "gene_start",
    "gene_end",
]

ANNOTATION_SUMMARY_COLUMNS: list[str] = [
    "annotation_class",
    "n_loci",
    "pct_total",
]

# Report output column schemas
RANKED_HITS_COLUMNS: list[str] = ["rank", *HITS_COLUMNS]

STRICTNESS_SUMMARY_COLUMNS: list[str] = [
    "strictness_class",
    "n_loci",
    "pct_hits",
]

SUPPORT_SUMMARY_COLUMNS: list[str] = [
    "source_type",
    "evidence_class",
    "n_records",
    "pct_of_source",
]


class TsvWriter:
    """Context-manager TSV writer with explicit column validation.

    Example::

        with TsvWriter(outdir / "hits.tsv", HITS_COLUMNS) as w:
            w.write_row({"locus_id": "PPX000001", "contig": "chr1", ...})
    """

    def __init__(self, path: Path, columns: list[str]) -> None:
        self.path = path
        self.columns = columns
        self._fh: Any = None
        self._writer: Any = None

    def __enter__(self) -> TsvWriter:
        self._fh = open(self.path, "w", newline="", encoding="utf-8")
        self._writer = csv.DictWriter(
            self._fh,
            fieldnames=self.columns,
            delimiter="\t",
            extrasaction="raise",
            lineterminator="\n",
        )
        self._writer.writeheader()
        return self

    def __exit__(self, *args: Any) -> None:
        if self._fh is not None:
            self._fh.close()
            self._fh = None
            self._writer = None

    def write_row(self, row: dict[str, Any]) -> None:
        """Write a single row.  Keys must match the declared column list."""
        if self._writer is None:
            raise RuntimeError("TsvWriter is not open.  Use as a context manager.")
        self._writer.writerow(row)

    def write_rows(self, rows: Iterable[dict[str, Any]]) -> None:
        """Write multiple rows."""
        for row in rows:
            self.write_row(row)

    def flush(self) -> None:
        """Flush buffered rows to disk."""
        if self._fh is not None:
            self._fh.flush()


def read_tsv(path: Path) -> list[dict[str, str]]:
    """Read a tab-separated file and return a list of row dicts."""
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return list(reader)
