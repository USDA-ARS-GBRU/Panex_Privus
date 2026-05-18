"""Extract focus-region genotype tables for interactive dashboards."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal

from privy.interactive.models import FocusRegion
from privy.io.gff import parse_gff3

VariantFilter = Literal["all", "snp", "indel", "sv"]


SITE_COLUMNS = [
    "contig",
    "pos",
    "ref",
    "alt",
    "variant_type",
    "qual",
    "filter",
    "{offtarget}_gt",
    "{offtarget}_alleles",
    "{offtarget}_state",
    "{offtarget}_carries_alt",
    "{derived}_gt",
    "{derived}_alleles",
    "{derived}_state",
    "{derived}_carries_alt",
    "{donor}_gt",
    "{donor}_alleles",
    "{donor}_state",
    "{donor}_carries_alt",
    "called_n",
    "alt_carrier_n",
    "derived_matches_donor_gt",
    "derived_matches_offtarget_gt",
    "donor_alt_absent_offtarget",
    "target_private_alt_pattern",
    "target_private_alt_indices",
    "target_private_alt_alleles",
    "overlapping_gene_ids",
    "overlapping_gene_names",
    "nearest_gene_id",
    "nearest_gene_name",
    "nearest_gene_distance_bp",
]


@dataclass(frozen=True)
class GeneInterval:
    gene_id: str
    name: str
    start: int
    end: int


@dataclass(frozen=True)
class ExtractionSummary:
    path: Path
    records_seen: int
    records_written: int
    skipped_filter: int
    skipped_variant_type: int
    skipped_missing: int
    skipped_biallelic: int


def extract_focus_sites_from_vcf(
    *,
    vcf: Path,
    gff3: Path | None,
    focus: FocusRegion,
    samples: tuple[str, str, str],
    out_tsv: Path,
    pass_only: bool = True,
    require_all_called: bool = True,
    variant_filter: VariantFilter = "all",
    biallelic_only: bool = False,
) -> ExtractionSummary:
    """Extract a focus-region sites TSV from a multisample VCF/BCF.

    The sample order is ``OFFTARGET, DERIVED, DONOR``. Coordinates in the output
    are 1-based VCF/GFF-style coordinates for direct display in the browser.
    """
    if variant_filter not in {"all", "snp", "indel", "sv"}:
        raise ValueError("--variant-type must be one of: all, snp, indel, sv.")

    genes = _read_gene_intervals(gff3, focus) if gff3 is not None else []
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    header = _site_columns(samples)

    records_seen = 0
    records_written = 0
    skipped_filter = 0
    skipped_variant_type = 0
    skipped_missing = 0
    skipped_biallelic = 0

    import pysam  # noqa: PLC0415

    with pysam.VariantFile(str(vcf)) as vf, out_tsv.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        _validate_samples(vf, samples)
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=header)
        writer.writeheader()
        for record in _iter_region_records(vf, focus):
            records_seen += 1
            alts = tuple(record.alts or ())
            if biallelic_only and len(alts) != 1:
                skipped_biallelic += 1
                continue
            if pass_only and not _record_passes(record):
                skipped_filter += 1
                continue
            record_type = _record_variant_type(record.ref, alts)
            if variant_filter != "all" and record_type != variant_filter:
                skipped_variant_type += 1
                continue
            sample_infos = {
                sample: _sample_info(record, sample)
                for sample in samples
            }
            if require_all_called and any(
                info["state"] == "no_call" for info in sample_infos.values()
            ):
                skipped_missing += 1
                continue
            writer.writerow(
                _site_row(
                    record=record,
                    alts=alts,
                    record_type=record_type,
                    samples=samples,
                    sample_infos=sample_infos,
                    genes=genes,
                )
            )
            records_written += 1

    return ExtractionSummary(
        path=out_tsv,
        records_seen=records_seen,
        records_written=records_written,
        skipped_filter=skipped_filter,
        skipped_variant_type=skipped_variant_type,
        skipped_missing=skipped_missing,
        skipped_biallelic=skipped_biallelic,
    )


def _site_columns(samples: tuple[str, str, str]) -> list[str]:
    offtarget, derived, donor = samples
    return [
        column.format(offtarget=offtarget, derived=derived, donor=donor)
        for column in SITE_COLUMNS
    ]


def _read_gene_intervals(path: Path, focus: FocusRegion) -> list[GeneInterval]:
    genes = []
    for record in parse_gff3(path, feature_types=frozenset({"gene"})):
        if record.seqid != focus.contig:
            continue
        start = record.start + 1
        end = record.end
        if end < focus.start or start > focus.end:
            continue
        gene_id = record.attrs.get("ID", record.gene_id)
        name = record.attrs.get("Name", record.gene_id or gene_id)
        genes.append(GeneInterval(gene_id=gene_id, name=name, start=start, end=end))
    return sorted(genes, key=lambda gene: (gene.start, gene.end, gene.name))


def _validate_samples(vf: Any, samples: tuple[str, str, str]) -> None:
    available = set(vf.header.samples)
    missing = [sample for sample in samples if sample not in available]
    if missing:
        raise ValueError(f"VCF is missing sample(s): {', '.join(missing)}")


def _iter_region_records(vf: Any, focus: FocusRegion):
    start0 = focus.start - 1
    try:
        yield from vf.fetch(focus.contig, start0, focus.end)
    except (OSError, ValueError):
        for record in vf:
            if record.chrom == focus.contig and focus.start <= record.pos <= focus.end:
                yield record


def _record_passes(record: Any) -> bool:
    filters = list(record.filter.keys())
    return not filters or filters == ["PASS"] or "PASS" in filters


def _filter_text(record: Any) -> str:
    filters = list(record.filter.keys())
    return ";".join(filters) if filters else "PASS"


def _record_variant_type(ref: str, alts: tuple[str, ...]) -> str:
    if any(_is_symbolic_alt(alt) for alt in alts):
        return "sv"
    if alts and len(ref) == 1 and all(len(alt) == 1 for alt in alts):
        return "snp"
    return "indel"


def _is_symbolic_alt(alt: str) -> bool:
    return alt.startswith("<") or "[" in alt or "]" in alt


def _sample_info(record: Any, sample: str) -> dict[str, Any]:
    call = record.samples[sample]
    gt = call.get("GT")
    phased = bool(getattr(call, "phased", False))
    alleles = _allele_strings(gt, record.ref, tuple(record.alts or ()))
    carries_alt = _carries_alt(gt)
    return {
        "gt": _format_gt(gt, phased),
        "gt_tuple": tuple(gt or ()),
        "alleles": "/".join(alleles) if alleles else ".",
        "state": _genotype_state(gt),
        "carries_alt": carries_alt,
        "alt_indices": _alt_indices(gt),
    }


def _format_gt(gt: tuple[int | None, ...] | None, phased: bool) -> str:
    if gt is None:
        return "."
    sep = "|" if phased else "/"
    return sep.join("." if allele is None else str(allele) for allele in gt)


def _allele_strings(
    gt: tuple[int | None, ...] | None,
    ref: str,
    alts: tuple[str, ...],
) -> list[str]:
    if gt is None:
        return []
    alleles = [ref, *alts]
    out = []
    for allele in gt:
        if allele is None:
            out.append(".")
        elif 0 <= allele < len(alleles):
            out.append(alleles[allele])
        else:
            out.append(f"<allele{allele}>")
    return out


def _genotype_state(gt: tuple[int | None, ...] | None) -> str:
    if gt is None or any(allele is None for allele in gt):
        return "no_call"
    unique = set(gt)
    if unique == {0}:
        return "hom_ref"
    if len(unique) == 1:
        return "hom_alt"
    if 0 in unique:
        return "het"
    return "alt_mixed"


def _carries_alt(gt: tuple[int | None, ...] | None) -> bool:
    return gt is not None and any(allele is not None and allele > 0 for allele in gt)


def _alt_indices(gt: tuple[int | None, ...] | None) -> set[int]:
    if gt is None:
        return set()
    return {allele for allele in gt if allele is not None and allele > 0}


def _site_row(
    *,
    record: Any,
    alts: tuple[str, ...],
    record_type: str,
    samples: tuple[str, str, str],
    sample_infos: dict[str, dict[str, Any]],
    genes: list[GeneInterval],
) -> dict[str, Any]:
    offtarget, derived, donor = samples
    offtarget_info = sample_infos[offtarget]
    derived_info = sample_infos[derived]
    donor_info = sample_infos[donor]
    shared_private = sorted(
        derived_info["alt_indices"]
        & donor_info["alt_indices"]
        - offtarget_info["alt_indices"]
    )
    donor_private = donor_info["alt_indices"] - offtarget_info["alt_indices"]
    overlapping, nearest = _gene_context(
        genes,
        start=record.pos,
        end=record.pos + max(1, len(record.ref)) - 1,
    )
    row: dict[str, Any] = {
        "contig": record.chrom,
        "pos": record.pos,
        "ref": record.ref,
        "alt": ",".join(alts),
        "variant_type": record_type,
        "qual": "." if record.qual is None else record.qual,
        "filter": _filter_text(record),
        "called_n": sum(info["state"] != "no_call" for info in sample_infos.values()),
        "alt_carrier_n": sum(bool(info["carries_alt"]) for info in sample_infos.values()),
        "derived_matches_donor_gt": _bool_text(
            derived_info["gt_tuple"] == donor_info["gt_tuple"]
            and derived_info["state"] != "no_call"
        ),
        "derived_matches_offtarget_gt": _bool_text(
            derived_info["gt_tuple"] == offtarget_info["gt_tuple"]
            and derived_info["state"] != "no_call"
        ),
        "donor_alt_absent_offtarget": _bool_text(bool(donor_private)),
        "target_private_alt_pattern": _bool_text(bool(shared_private)),
        "target_private_alt_indices": (
            ",".join(str(index) for index in shared_private) if shared_private else "NA"
        ),
        "target_private_alt_alleles": (
            ",".join(alts[index - 1] for index in shared_private) if shared_private else "NA"
        ),
        "overlapping_gene_ids": (
            ",".join(gene.gene_id for gene in overlapping) if overlapping else "NA"
        ),
        "overlapping_gene_names": (
            ",".join(gene.name for gene in overlapping) if overlapping else "NA"
        ),
        "nearest_gene_id": nearest.gene_id if nearest else "NA",
        "nearest_gene_name": nearest.name if nearest else "NA",
        "nearest_gene_distance_bp": (
            _distance_to_gene(nearest, record.pos, record.pos + max(1, len(record.ref)) - 1)
            if nearest
            else "NA"
        ),
    }
    for sample in samples:
        info = sample_infos[sample]
        row[f"{sample}_gt"] = info["gt"]
        row[f"{sample}_alleles"] = info["alleles"]
        row[f"{sample}_state"] = info["state"]
        row[f"{sample}_carries_alt"] = _bool_text(bool(info["carries_alt"]))
    return row


def _gene_context(
    genes: list[GeneInterval],
    start: int,
    end: int,
) -> tuple[list[GeneInterval], GeneInterval | None]:
    overlapping = [gene for gene in genes if gene.start <= end and gene.end >= start]
    if overlapping:
        return overlapping, sorted(overlapping, key=lambda gene: (gene.start, gene.end))[0]
    nearest = min(
        genes,
        key=lambda gene: _distance_to_gene(gene, start, end),
        default=None,
    )
    return [], nearest


def _distance_to_gene(gene: GeneInterval, start: int, end: int) -> int:
    if gene.start <= end and gene.end >= start:
        return 0
    if end < gene.start:
        return gene.start - end
    return start - gene.end


def _bool_text(value: bool) -> str:
    return "true" if value else "false"
