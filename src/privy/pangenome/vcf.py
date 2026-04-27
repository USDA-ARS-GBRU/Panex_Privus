"""VCF adapter for the shared pangenome feature matrix."""

from __future__ import annotations

from pathlib import Path

from privy.io.vcf import (
    classify_variant_type,
    format_allele_key,
    get_vcf_samples,
    has_alt_allele,
    is_missing_genotype,
    stream_vcf_records,
)
from privy.pangenome.model import FeatureMatrix, FeatureRecord


def build_vcf_feature_matrix(vcf_path: Path) -> FeatureMatrix:
    """Build an allele-level feature matrix from a multisample VCF.

    Each alternate allele is treated as one pangenome feature.  Presence means
    a sample has at least one copy of that specific ALT allele in its genotype.
    Missing genotypes and reference/different-ALT genotypes are not present.
    """
    samples = tuple(get_vcf_samples(vcf_path))
    features: list[FeatureRecord] = []
    presence: dict[str, frozenset[str]] = {}

    for record in stream_vcf_records(vcf_path):
        if record.alts is None:
            continue
        start = record.pos - 1
        end = start + max(len(record.ref), 1)
        for alt_index, alt in enumerate(record.alts):
            feature_id = format_allele_key(record.chrom, record.pos, record.ref, alt)
            present_samples: set[str] = set()
            for sample in samples:
                try:
                    gt = record.samples[sample]["GT"]
                except (KeyError, TypeError):
                    continue
                if is_missing_genotype(gt):
                    continue
                if gt is not None and has_alt_allele(gt, alt_index):
                    present_samples.add(sample)

            features.append(
                FeatureRecord(
                    feature_id=feature_id,
                    source_type="vcf",
                    feature_type=classify_variant_type(record.ref, alt),
                    length=max(len(record.ref), len(alt), 1),
                    contig=record.chrom,
                    start=start,
                    end=end,
                )
            )
            presence[feature_id] = frozenset(present_samples)

    return FeatureMatrix(
        source_type="vcf",
        features=tuple(features),
        samples=samples,
        presence=presence,
    )
