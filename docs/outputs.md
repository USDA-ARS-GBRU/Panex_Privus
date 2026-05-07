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

`privy pangenome` writes pangenome-wide and sub-pangenome summaries to the
directory you choose with `--outdir`.

| File | Purpose |
|------|---------|
| `feature_summary.tsv` | One row per pangenome feature with full, target, and off-target counts |
| `coverage_histogram.tsv` | Coverage histogram by feature count and bp |
| `composition.tsv` | Core, accessory, private, and absent feature totals |
| `growth_curves.tsv` | Permutation-based pangenome growth data |
| `pangenome_growth.png` | Growth curves for full, target, and off-target groups |
| `pangenome_coverage.png` | Feature coverage distribution |
| `pangenome_composition.png` | Group composition plot |
| `pangenome.json` | Run metadata, resolved groups, parameters, and output list |

`privy landscape` writes VCF sliding-window summaries to the directory you
choose with `--outdir`.

| File | Purpose |
|------|---------|
| `sample_windows.tsv` | One row per sample per window with missingness, genotype burden, private/rare ALT burden, and nearest local background |
| `windows.tsv` | One row per window with target/off-target summary metrics |
| `background_blocks.tsv` | Adjacent sample windows merged by nearest local background assignment |
| `similarity.tsv` | Pairwise sample genotype similarity for each window |
| `missingness_heatmap.png` | Sample-by-window missingness heatmap |
| `private_burden_heatmap.png` | Sample-by-window private ALT burden heatmap |
| `local_background_map.png` | Sample-by-window nearest-background map |
| `similarity_cluster_map.png` | Clustered mean sample-similarity heatmap |
| `landscape.json` | Run metadata, resolved groups, parameters, and output list |

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

## Pangenome Outputs

`feature_summary.tsv` uses a shared feature model. In the GFA adapter, each
feature is a graph segment. In the VCF adapter, each feature is one alternate
allele from a variant record.

For example snippets and publication-style captions for the pangenome tables
and plots, see [Figures and Tables]({{ '/figures-and-tables/' | relative_url }}).

| Column | Meaning |
|--------|---------|
| `feature_id` | Source-specific feature identifier |
| `source_type` | Input source, `gfa` or `vcf` |
| `feature_type` | Feature kind, such as `segment`, `snp`, `indel`, or `sv` |
| `contig`, `start`, `end` | Coordinates when available |
| `length` | Feature length used for bp-weighted summaries |
| `total_present_n` | Number of active samples containing the feature |
| `target_present_n` | Number of target samples containing the feature |
| `offtarget_present_n` | Number of off-target samples containing the feature |
| `full_category`, `target_category`, `offtarget_category` | `absent`, `private`, `accessory`, or `core` |
| `target_private` | `True` when present in targets and absent from off-targets |
| `offtarget_private` | `True` when present in off-targets and absent from targets |

## Landscape Outputs

`privy landscape` uses 0-based half-open coordinates and reports both window
mode and window size parameters in `landscape.json`.

### `sample_windows.tsv`

| Column | Meaning |
|--------|---------|
| `window_id` | Stable ID for the emitted landscape window |
| `contig`, `start`, `end`, `midpoint` | Window coordinates |
| `window_mode` | `records` or `bp` |
| `n_variants` | Number of VCF records in the window |
| `sample` | Sample name |
| `cohort_role` | `target` or `off_target` |
| `missing_rate` | Fraction of records where the sample genotype is missing |
| `het_rate` | Fraction of called records with heterozygous genotype |
| `nonref_rate` | Fraction of called records carrying any ALT allele |
| `minor_genotype_rate` | Fraction of called records where the sample has a minor genotype class |
| `rare_alt_rate` | Rare ALT allele events carried by the sample per called record |
| `private_alt_rate` | Cohort-private ALT allele events carried by the sample per called record |
| `median_call_freq` | Median frequency of the sample's genotype class in the window |
| `nearest_background` | Most similar sample in the same window |
| `nearest_similarity` | Genotype-match fraction to the nearest background sample |

### `windows.tsv`

| Column | Meaning |
|--------|---------|
| `span_bp` | Physical span of the emitted window |
| `density_variants_per_kb` | VCF record density across the window span |
| `target_mean_missing_rate` | Mean missingness across target samples |
| `offtarget_mean_missing_rate` | Mean missingness across off-target samples |
| `target_mean_nonref_rate` | Mean non-reference burden across target samples |
| `offtarget_mean_nonref_rate` | Mean non-reference burden across off-target samples |
| `target_private_alt_n` | Number of ALT alleles carried by targets and absent from off-targets |
| `offtarget_private_alt_n` | Number of ALT alleles carried by off-targets and absent from targets |
| `target_private_alt_rate` | Target-private ALT events per VCF record in the window |
| `offtarget_private_alt_rate` | Off-target-private ALT events per VCF record in the window |
| `top_nearest_background` | Most frequent nearest-background assignment in the window |

### `background_blocks.tsv`

| Column | Meaning |
|--------|---------|
| `block_id` | Stable local-background block ID |
| `sample` | Sample being assigned |
| `cohort_role` | Sample role |
| `contig`, `start`, `end` | Merged block coordinates |
| `n_windows` | Number of windows merged into the block |
| `nearest_background` | Assigned nearest local background, or `unassigned` |
| `mean_similarity` | Mean nearest-background similarity across merged windows |

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
