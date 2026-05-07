"""Runtime statistics accumulators for Panex Privus.

Tracks per-run counts (records seen, loci emitted, loci merged, etc.)
for writing into qc.tsv and run.json.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class ScanStats:
    """Accumulated statistics from a single ``privy scan`` run.

    Populated by the VCF scan backend and written to ``qc.tsv`` and
    ``run.json``.
    """

    records_evaluated: int = 0
    """Total VCF records processed (before any filtering)."""

    records_skipped_filter: int = 0
    """Records skipped because FILTER != PASS (when pass_only=True)."""

    records_skipped_qual: int = 0
    """Records skipped because QUAL < min_qual."""

    records_skipped_multiallelic: int = 0
    """Records skipped because allow_multiallelic=False."""

    alleles_evaluated: int = 0
    """Total alternate alleles submitted to classify_strictness."""

    alleles_passed: int = 0
    """Alleles that passed discovery criteria (pattern_pass=True)."""

    alleles_contradicted: int = 0
    """Alleles classified as contradicted."""

    alleles_missing_only: int = 0
    """Alleles where all samples in one cohort are missing."""

    loci_emitted: int = 0
    """Loci written to hits.tsv."""

    regions_emitted: int = 0
    """Regions written to regions.tsv."""

    strictness_counts: dict[str, int] = field(default_factory=dict)
    """Per-strictness-class counts of passing loci."""

    n_target_samples: int = 0
    """Number of target samples found in VCF."""

    n_offtarget_samples: int = 0
    """Number of off-target samples found in VCF."""

    n_contigs_scanned: int = 0
    """Number of contigs visited."""

    def increment_strictness(self, class_name: str) -> None:
        """Increment the counter for *class_name*."""
        self.strictness_counts[class_name] = (
            self.strictness_counts.get(class_name, 0) + 1
        )

    def as_qc_rows(self, source: str = "vcf") -> list[dict[str, str]]:
        """Return a list of row dicts for ``qc.tsv``."""
        if source == "gfa":
            record_desc = "GFA graph segments processed for target-private status"
            skipped_filter_desc = "GFA records skipped by filter"
            skipped_qual_desc = "GFA records skipped by quality threshold"
            skipped_multiallelic_desc = "GFA records skipped as multiallelic"
            allele_eval_desc = "Graph segments evaluated for private-segment status"
            allele_pass_desc = "Graph segments passing discovery criteria"
            allele_contradicted_desc = "Graph segments classified as contradicted"
            target_desc = "Target samples found in GFA"
            offtarget_desc = "Off-target samples found in GFA"
        else:
            record_desc = "Total VCF records processed"
            skipped_filter_desc = "Records skipped: FILTER != PASS"
            skipped_qual_desc = "Records skipped: QUAL below threshold"
            skipped_multiallelic_desc = (
                "Records skipped: multiallelic (when allow_multiallelic=False)"
            )
            allele_eval_desc = "Alternate alleles evaluated for private-allele status"
            allele_pass_desc = "Alleles passing discovery criteria"
            allele_contradicted_desc = "Alleles classified as contradicted"
            target_desc = "Target samples found in VCF header"
            offtarget_desc = "Off-target samples found in VCF header"

        rows: list[dict[str, str]] = [
            {"metric": "records_evaluated",       "value": str(self.records_evaluated),
             "description": record_desc},
            {"metric": "records_skipped_filter",  "value": str(self.records_skipped_filter),
             "description": skipped_filter_desc},
            {"metric": "records_skipped_qual",    "value": str(self.records_skipped_qual),
             "description": skipped_qual_desc},
            {
                "metric": "records_skipped_multiallelic",
                "value": str(self.records_skipped_multiallelic),
                "description": skipped_multiallelic_desc,
            },
            {"metric": "alleles_evaluated",       "value": str(self.alleles_evaluated),
             "description": allele_eval_desc},
            {"metric": "alleles_passed",          "value": str(self.alleles_passed),
             "description": allele_pass_desc},
            {"metric": "alleles_contradicted",    "value": str(self.alleles_contradicted),
             "description": allele_contradicted_desc},
            {"metric": "loci_emitted",            "value": str(self.loci_emitted),
             "description": "Loci written to hits.tsv"},
            {"metric": "regions_emitted",         "value": str(self.regions_emitted),
             "description": "Regions written to regions.tsv"},
            {"metric": "n_target_samples",        "value": str(self.n_target_samples),
             "description": target_desc},
            {"metric": "n_offtarget_samples",     "value": str(self.n_offtarget_samples),
             "description": offtarget_desc},
            {"metric": "n_contigs_scanned",       "value": str(self.n_contigs_scanned),
             "description": "Contigs visited during scan"},
        ]
        # Add per-strictness-class counts
        for class_name, count in sorted(self.strictness_counts.items()):
            rows.append({
                "metric": f"strictness_{class_name}",
                "value": str(count),
                "description": f"Passing loci classified as {class_name}",
            })
        return rows

    def as_summary_dict(self) -> dict[str, object]:
        """Return a JSON-serialisable summary dict for ``run.json``."""
        return {
            "records_evaluated": self.records_evaluated,
            "records_skipped_filter": self.records_skipped_filter,
            "records_skipped_qual": self.records_skipped_qual,
            "alleles_evaluated": self.alleles_evaluated,
            "alleles_passed": self.alleles_passed,
            "alleles_contradicted": self.alleles_contradicted,
            "loci_emitted": self.loci_emitted,
            "regions_emitted": self.regions_emitted,
            "n_target_samples": self.n_target_samples,
            "n_offtarget_samples": self.n_offtarget_samples,
            "n_contigs_scanned": self.n_contigs_scanned,
            "strictness_counts": dict(self.strictness_counts),
        }
