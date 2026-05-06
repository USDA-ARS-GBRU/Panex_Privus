---
title: Troubleshooting
description: Common Panex Privus setup and analysis issues.
---

# Troubleshooting

## "VCF index not found"

Your VCF must have a `.tbi` or `.csi` index next to it.

```bash
tabix -p vcf variants.vcf.gz
```

## "No target samples from the cohort definition were found in the VCF header"

The sample names passed to `--targets` do not match the VCF header.

List sample names:

```bash
bcftools query -l variants.vcf.gz
```

Sample names are case-sensitive.

## `hits.tsv` Is Empty

Check `qc.tsv` first. Common causes:

- all records were skipped by FILTER
- the target support threshold is too strict
- off-target samples carry the same allele
- cohort labels do not match the input sample names

Try:

```bash
privy scan --vcf variants.vcf.gz --no-pass-only ...
```

or lower `--min-target-support` if partial target support is biologically
acceptable.

## GFA Scan Produces No Hits

The examples below assume a compressed `pangenome.gfa.gz`. For plain `.gfa`,
replace `gzip -cd pangenome.gfa.gz` with `cat pangenome.gfa`.

Check that GFA segments have coordinate tags:

```bash
gzip -cd pangenome.gfa.gz | grep "^S" | head -3
```

Expected tags look like:

```text
SN:Z:chr1  SO:i:1000  LN:i:500
```

Also verify sample names from W-lines:

```bash
gzip -cd pangenome.gfa.gz | grep "^W" | awk '{print $2}' | sort -u
```

For P-lines:

```bash
gzip -cd pangenome.gfa.gz | grep "^P" | awk '{print $2}' | cut -d'#' -f1 | sort -u
```

## GFA Scan Is Killed While Parsing

Update to the newest development version first. Current `privy scan --gfa`
builds a single-pass scan-specific streaming index, logs indexing progress, and
does not retain full GFA sequences, links, walks, or paths in memory.

```bash
git pull origin main
python -m pip install -U .
```

Then run directly on the compressed minigraph-cactus graph:

```bash
privy scan --gfa pangenome.gfa.gz --targets T1 T2 --off-targets O1 O2 --outdir results/
```

## BAM Support Is Uninformative

Common reasons:

- BAM file is not indexed
- contig names do not match the VCF
- depth is below `--bam-min-depth`
- the locus is an indel or structural allele where SNP-style pileup cannot
  confidently count the alternate allele

## pysam Installation Errors

Install `pysam` with conda first:

```bash
conda install -c bioconda pysam
pip install .
```
