---
title: Quickstart
description: First commands for scanning, comparing, reporting, plotting, annotating, and exporting Panex Privus results.
---

# Quickstart

This walkthrough takes you from an input file to a ranked list of candidate
private regions.

## VCF Scan

Make sure your VCF is compressed and indexed:

```bash
bgzip -c variants.vcf > variants.vcf.gz
tabix -p vcf variants.vcf.gz
```

Run a scan:

```bash
privy scan \
  --vcf variants.vcf.gz \
  --targets T1 T2 T3 \
  --off-targets O1 O2 O3 \
  --outdir results/vcf/
```

Important outputs:

- `hits.tsv`: ranked candidate private alleles
- `regions.tsv`: merged candidate regions
- `sample_support.tsv`: per-sample genotype/evidence table
- `qc.tsv`: scan metrics
- `run.json`: run metadata and resolved configuration

## GFA Scan

Run a graph scan:

```bash
privy scan \
  --gfa pangenome.gfa \
  --targets T1 T2 T3 \
  --off-targets O1 O2 O3 \
  --outdir results/gfa/
```

GFA segments must have coordinate tags such as `SN:Z:chr1`, `SO:i:1000`, and
`LN:i:500`. Minigraph-cactus output usually includes these tags.

## Add BAM Evidence to a VCF Scan

BAM files must be coordinate-sorted and indexed:

```bash
samtools sort sample.bam -o sample.sorted.bam
samtools index sample.sorted.bam
```

Run with repeated `--bam` flags:

```bash
privy scan \
  --vcf variants.vcf.gz \
  --targets T1 T2 \
  --off-targets O1 O2 \
  --bam T1.sorted.bam \
  --bam O1.sorted.bam \
  --outdir results/vcf_bam/
```

Or use a manifest:

```text
bam_path	sample_id
T1.sorted.bam	T1
O1.sorted.bam	O1
```

```bash
privy scan \
  --vcf variants.vcf.gz \
  --targets T1 T2 \
  --off-targets O1 O2 \
  --bam-manifest manifest.tsv \
  --outdir results/vcf_bam/
```

## Compare VCF and GFA Scans

```bash
privy compare \
  --hits-a results/vcf/hits.tsv \
  --hits-b results/gfa/hits.tsv \
  --outdir results/compare/
```

## Report, Plot, Annotate, Export

```bash
privy report --hits results/vcf/hits.tsv --regions results/vcf/regions.tsv \
  --qc results/vcf/qc.tsv --format both --outdir results/report/

privy plot --hits results/vcf/hits.tsv --evidence results/vcf/evidence.tsv \
  --outdir results/plots/

privy annotate --hits results/vcf/hits.tsv --gff annotation.gff3.gz \
  --outdir results/annotated/

privy export --hits results/vcf/hits.tsv --regions results/vcf/regions.tsv \
  --format gff3 --outdir results/exported/
```
