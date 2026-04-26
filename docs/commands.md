---
title: Command Reference
description: CLI reference for Panex Privus commands.
---

# Command Reference

## Global Options

```bash
privy --help
```

Common options:

| Option | Description |
|--------|-------------|
| `--config PATH` | YAML configuration file |
| `--project-name TEXT` | Project name written into outputs |
| `--outdir PATH` | Default output directory |
| `--threads INT` | Reserved for supported parallel paths |
| `--log-level TEXT` | `debug`, `info`, `warning`, or `error` |
| `--quiet` | Reduce console output |
| `--version` | Show package version |

## `privy scan`

Primary discovery command.

```bash
privy scan --vcf variants.vcf.gz --targets T1 T2 --off-targets O1 O2 --outdir results/
privy scan --gfa pangenome.gfa --targets T1 T2 --off-targets O1 O2 --outdir results/
privy scan --vcf variants.vcf.gz --gfa pangenome.gfa \
  --targets T1 T2 --off-targets O1 O2 --outdir results/
```

Scan outputs are written to source subdirectories: `results/vcf/` for VCF
discovery, `results/gfa/` for GFA discovery, and `results/compare/` when both
inputs are provided.

Key options:

| Option | Description |
|--------|-------------|
| `--vcf PATH` | Indexed multisample VCF |
| `--gfa PATH` | GFA graph file |
| `--bam PATH` | BAM support file; repeat for multiple |
| `--bam-manifest PATH` | TSV mapping BAM files to sample IDs |
| `--targets TEXT` | Target samples |
| `--off-targets TEXT` | Off-target samples |
| `--cohort-file PATH` | YAML or TSV cohort definition |
| `--region TEXT` | Restrict to `contig:start-end` |
| `--contig TEXT` | Restrict to one contig |
| `--merge-distance INT` | Merge nearby hits into regions |
| `--min-segment-length INT` | Minimum GFA segment length |

## `privy compare`

Compare two scan result sets.

```bash
privy compare --hits-a results/vcf/hits.tsv --hits-b results/gfa/hits.tsv \
  --outdir results/compare/
```

Key options:

| Option | Description |
|--------|-------------|
| `--hits-a PATH` | First `hits.tsv` |
| `--hits-b PATH` | Second `hits.tsv` |
| `--source-a TEXT` | Optional label for source A |
| `--source-b TEXT` | Optional label for source B |
| `--min-reciprocal-overlap FLOAT` | Minimum overlap for interval matching |
| `--breakpoint-tolerance-bp INT` | Gap tolerance for near misses |
| `--require-state-compatibility` | Require strictness compatibility |

## `privy report`

Generate Markdown and/or HTML reports from existing outputs.

```bash
privy report --hits results/vcf/hits.tsv --regions results/vcf/regions.tsv \
  --qc results/vcf/qc.tsv --format both --outdir report/
```

## `privy plot`

Generate diagnostic figures.

```bash
privy plot --hits results/vcf/hits.tsv --evidence results/vcf/evidence.tsv --outdir plots/
```

Plot types:

- `locus_panel`
- `strictness_bar`
- `score_distribution`
- `support_bar`
- `compare_summary`
- `all`

## `privy annotate`

Intersect hits with GFF3 annotation.

```bash
privy annotate --hits results/vcf/hits.tsv --gff annotation.gff3.gz --outdir annotated/
```

## `privy export`

Export hits and regions to BED or GFF3.

```bash
privy export --hits results/vcf/hits.tsv --regions results/vcf/regions.tsv \
  --format bed --outdir exported/

privy export --hits results/vcf/hits.tsv --regions results/vcf/regions.tsv \
  --format gff3 --outdir exported_gff3/
```
