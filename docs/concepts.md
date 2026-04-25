---
title: Core Concepts
description: Key Panex Privus ideas including target-private signal, strictness, and evidence classes.
---

# Core Concepts

## Target-Private Signal

A target-private allele or segment is present in the target cohort and absent
from the off-target cohort.

For example, if five soybean lines share a trait and ten other lines do not, a
natural question is:

> What DNA variants are present in the five trait-positive lines and absent
> from the ten reference lines?

Panex Privus automates that search from VCF genotype calls or GFA graph
traversal.

## Primary Discovery Inputs

VCF and GFA are independent primary discovery backends.

| Input | Role | Requirement |
|-------|------|-------------|
| VCF | Genotype-call-based allele discovery | bgzip-compressed `.vcf.gz` plus `.tbi` or `.csi` index |
| GFA | Graph-traversal-based segment discovery | Plain text GFA with coordinate tags |

If both `--vcf` and `--gfa` are supplied to one scan, the VCF backend runs and
the GFA argument is ignored. To compare VCF and GFA evidence, run two separate
scans and use `privy compare`.

## BAM Support

BAM is a support layer, not a discovery caller. It queries read-level evidence at
loci already discovered from a VCF scan.

For each candidate locus, BAM support can add:

- read depth
- reference and alternate allele counts
- allele fraction
- support, absence, contradiction, ambiguous, or uninformative evidence classes

## Missingness and Strictness

Privy never silently treats missing data as absence. Every passing hit receives a
`strictness_class`.

| Class | Meaning |
|-------|---------|
| `strict_complete` | All targets support; all off-targets are confidently absent; no missing calls |
| `strict_target_missing` | Off-target exclusion holds, but at least one target is missing |
| `strict_offtarget_missing` | Target support holds, but at least one off-target is missing |
| `strict_both_missing` | Pattern is consistent, but both groups have missing data |
| `relaxed_threshold` | Passes user thresholds but not complete strict logic |
| `contradicted` | Target-private model fails |

In practice, start with `strict_complete`, then inspect missingness-aware classes
for biologically interesting candidates that need follow-up validation.
