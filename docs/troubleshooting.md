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

Check that GFA segments have coordinate tags:

```bash
grep "^S" pangenome.gfa | head -3
```

Expected tags look like:

```text
SN:Z:chr1  SO:i:1000  LN:i:500
```

Also verify sample names from W-lines:

```bash
grep "^W" pangenome.gfa | awk '{print $2}' | sort -u
```

For P-lines:

```bash
grep "^P" pangenome.gfa | awk '{print $2}' | cut -d'#' -f1 | sort -u
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
