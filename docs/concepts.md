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

If both `--vcf` and `--gfa` are supplied to one scan, Panex Privus writes
source-specific VCF and GFA result directories plus comparison outputs. You can
also run the two scans separately and use `privy compare` later.

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

## Pangenome Summaries

`privy pangenome` describes the feature space behind discovery. In the GFA
adapter, each feature is a graph segment. In the VCF adapter, each feature is
one alternate allele. The command summarizes how many features are core,
accessory, private, or absent in the full cohort, target cohort, and off-target
cohort.

Use pangenome summaries when you want to ask:

> How different are the target and off-target sub-pangenomes overall?

This is different from `privy scan`, which asks which individual loci match the
target-private pattern strongly enough to become candidates.

## VCF Landscapes

`privy landscape` is a windowed VCF context layer. It does not replace
population-genetic tools such as VCFtools, pixy, scikit-allel, PLINK, or
R/qtl2. Its novelty is Panex Privus' target/off-target framing: the same cohort
definition used for discovery is applied to genome-wide windows.

Landscape windows answer questions such as:

- Are candidate regions surrounded by high missingness?
- Are target-private alleles concentrated in particular chromosome intervals?
- Which samples are locally most similar to each other?
- Do adjacent windows form local background blocks that suggest shared genomic
  background?

The first implementation supports fixed-record windows by default and
base-pair windows when requested. Fixed-record windows keep variant counts more
stable across uneven SNP density. Base-pair windows are often easier to explain
on chromosome-scale figures.

## Local Background Blocks

A local background block is a run of adjacent windows where a sample's nearest
genotypic neighbor stays the same, subject to a minimum similarity threshold.
These blocks are useful for exploratory maps of shared genomic background.

They should not automatically be interpreted as a formal recombination-rate map.
Formal genetic maps usually need a cross design, progeny, marker order, and a
model such as those used by QTL/genetic-map software. In Panex Privus, local
background blocks are best read as:

> This sample looks locally most similar to that sample or group across this
> chromosome interval.

For controlled crosses, MAGIC populations, or founder-aware designs, these
blocks can become a bridge to more formal recombination or founder-haplotype
analyses.
