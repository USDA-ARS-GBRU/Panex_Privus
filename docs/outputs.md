---
title: Output Files
description: Panex Privus output files and important result columns.
---

# Output Files

`privy scan` writes source-specific output directories. VCF discovery results
go under `results/vcf/`, GFA discovery results go under `results/gfa/`, and a
combined VCF+GFA run also writes reconciliation files under `results/compare/`.
Each scan directory contains six primary outputs.

| File | Purpose |
|------|---------|
| `hits.tsv` | One row per candidate private locus, sorted by confidence score |
| `regions.tsv` | Nearby hits merged into candidate regions |
| `evidence.tsv` | Evidence records from VCF/GFA and optional BAM |
| `sample_support.tsv` | Per-sample support/missingness table |
| `qc.tsv` | Scan metrics |
| `run.json` | Run metadata and resolved configuration |

## Key `hits.tsv` Columns

| Column | Meaning |
|--------|---------|
| `locus_id` | Unique locus ID (`PPX...` for VCF, `GPX...` for GFA) |
| `contig`, `start`, `end` | 0-based half-open coordinates |
| `variant_type` | `snp`, `indel`, `sv`, or `graph_region` |
| `allele_key` | VCF allele key or GFA segment key |
| `target_support_n` | Number of target samples supporting the allele/segment |
| `offtarget_support_n` | Number of off-target samples supporting it |
| `target_missing_n` | Target samples with missing/uninformative data |
| `offtarget_missing_n` | Off-target samples with missing/uninformative data |
| `strictness_class` | Missingness-aware confidence class |
| `final_score` | Combined ranking score |

Coordinates use the BED/pysam convention: 0-based, half-open `[start, end)`.
For a SNP at VCF position 12345, `start=12344` and `end=12345`.

## Report Outputs

`privy report` can write:

- `summary.tsv`
- `ranked_hits.tsv`
- `strictness_summary.tsv`
- `support_summary.tsv`
- `contradiction_summary.tsv`
- `report.md`
- `report.html`

## Compare Outputs

`privy compare` writes:

- `compare.tsv`
- `compare_summary.tsv`
- `compare.json`

Match classes include `supported`, `partially_supported`, `contradicted`,
`source_specific`, `uninformative`, and `missing_data`.

## Annotate Outputs

`privy annotate` writes:

- `annotated_hits.tsv`
- `annotation_summary.tsv`
- `annotate.json`

Annotation classes are `CDS`, `UTR`, `exonic`, `intronic`, and `intergenic`.

## Export Outputs

`privy export` writes:

- `hits.bed` / `regions.bed`
- `hits.gff3` / `regions.gff3`
- `export.json`

BED scores are scaled to the standard 0-1000 range. GFF3 exports convert Privy
coordinates to 1-based closed feature coordinates.
