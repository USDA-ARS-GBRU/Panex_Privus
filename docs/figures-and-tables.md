---
title: Figures and Tables
description: How to interpret Panex Privus pangenome outputs and adapt them for publications.
---

# Figures and Tables

Panex Privus writes pangenome outputs as both machine-readable tables and
publication-starting plots. The tables are the analysis record. The plots are
interpretive summaries that help you explain full, target, and off-target
pangenomes in a manuscript, report, or slide deck.

The examples below were generated from the small test fixtures included in the
repository:

```bash
privy pangenome \
  --gfa tests/data/small_cohort.gfa \
  --targets T1 \
  --targets T2 \
  --permutations 25 \
  --outdir docs/assets/examples/pangenome-gfa

privy pangenome \
  --vcf tests/data/small_cohort.vcf \
  --targets T1 \
  --targets T2 \
  --permutations 25 \
  --outdir docs/assets/examples/pangenome-vcf
```

In both examples, `T1` and `T2` are the target samples. Because no off-targets
are specified, every other sample in the input becomes off-target.

## Feature Summary

`feature_summary.tsv` is the main feature-level table. For GFA input, each row
is a graph segment. For VCF input, each row is one alternate allele.

Small GFA example:

| feature_id | source_type | feature_type | total_present_n | target_present_n | offtarget_present_n | target_category | offtarget_category | target_private |
|------------|-------------|--------------|-----------------|------------------|---------------------|-----------------|--------------------|----------------|
| `s1` | `gfa` | `segment` | 5 | 2 | 3 | `core` | `core` | `False` |
| `s2_target` | `gfa` | `segment` | 2 | 2 | 0 | `core` | `absent` | `True` |
| `s4_target` | `gfa` | `segment` | 1 | 1 | 0 | `private` | `absent` | `True` |
| `s5` | `gfa` | `segment` | 4 | 1 | 3 | `private` | `core` | `False` |

Small VCF example:

| feature_id | source_type | feature_type | total_present_n | target_present_n | offtarget_present_n | target_category | offtarget_category | target_private |
|------------|-------------|--------------|-----------------|------------------|---------------------|-----------------|--------------------|----------------|
| `chr1:100:A:T` | `vcf` | `snp` | 2 | 2 | 0 | `core` | `absent` | `True` |
| `chr1:500:A:T` | `vcf` | `snp` | 3 | 2 | 1 | `core` | `private` | `False` |
| `chr1:800:A:G` | `vcf` | `snp` | 0 | 0 | 0 | `absent` | `absent` | `False` |
| `chr1:900:AGG:A` | `vcf` | `indel` | 2 | 2 | 0 | `core` | `absent` | `True` |

Use this table when you need to report specific graph segments or VCF alleles.
For publication methods, state the feature type: GFA segments and VCF alleles
are analyzed through the same matrix, but they are not biologically identical
objects.

## Composition

`composition.tsv` summarizes feature categories for each group.

| group | category | n_features | n_bp |
|-------|----------|------------|------|
| `full` | `private` | 1 | 7 |
| `full` | `accessory` | 4 | 35 |
| `full` | `core` | 2 | 16 |
| `target` | `absent` | 2 | 17 |
| `target` | `private` | 2 | 15 |
| `target` | `core` | 3 | 26 |

Interpretation:

- `core`: present in every sample in that group.
- `private`: present in exactly one sample in that group.
- `accessory`: present in more than one sample but not all samples.
- `absent`: not present in that group.

For target/off-target analysis, `target_private=True` in `feature_summary.tsv`
is often the most direct flag. The composition plot is broader: it describes
the shape of each group’s pangenome, not only candidate target-private signal.

## Pangenome Growth

![Pangenome growth example]({{ '/assets/examples/pangenome-gfa/pangenome_growth.png' | relative_url }})

**Figure caption example.** Pangenome growth curves for the full cohort, target
sub-pangenome, and off-target sub-pangenome. Curves show the mean number of
observed GFA segment features as samples are added across 25 deterministic
permutations. Shaded intervals show the 2.5th to 97.5th percentile range across
permutations. The target and off-target curves are computed independently from
the same graph-derived feature matrix.

Use this figure to show whether the feature set is still expanding as more
samples are added. A curve that continues rising steeply suggests additional
sampling may reveal more features. A curve approaching a plateau suggests the
observed cohort is closer to saturating the feature space under the chosen input
representation.

## Coverage Distribution

![Pangenome coverage example]({{ '/assets/examples/pangenome-gfa/pangenome_coverage.png' | relative_url }})

**Figure caption example.** Feature coverage distribution for full, target, and
off-target groups. The x-axis gives the number of samples containing a feature,
and the y-axis gives the number of features at that coverage level. In GFA
analysis, features are graph segments; in VCF analysis, features are alternate
alleles.

Use this figure to distinguish singleton/private-rich datasets from datasets
dominated by features shared across many samples. For larger cohorts, this plot
is useful for spotting heavy private-feature tails, highly conserved cores, or
unexpected representation imbalance between groups.

## Pangenome Composition

![Pangenome composition example]({{ '/assets/examples/pangenome-gfa/pangenome_composition.png' | relative_url }})

**Figure caption example.** Core, accessory, private, and absent feature counts
for the full cohort, target sub-pangenome, and off-target sub-pangenome.
Categories are assigned independently within each group. A feature can therefore
be target-core and off-target-absent, which is the strongest pangenome-level
pattern for target-private signal.

Use this figure to communicate the structure of each group’s pangenome. It is
especially useful when the target group has a distinct sub-pangenome profile
relative to the background group.

## Reporting In Methods

A reproducible methods paragraph should include:

- input source: GFA graph segments, VCF alternate alleles, or both
- sample counts and names or cohort-definition file
- whether off-targets were provided explicitly or inferred from remaining input
  samples
- number of permutations and random seed used for growth curves
- software version or commit
- any upstream filtering performed before creating the GFA or VCF

Example language:

> We summarized pangenome composition with Panex Privus using graph segments as
> features. Target samples were `T1` and `T2`; all remaining graph samples were
> treated as off-targets. Growth curves were estimated from 25 deterministic
> sample-order permutations with seed 42. Features were classified as core,
> accessory, private, or absent independently for the full, target, and
> off-target groups.

For VCF analyses, replace "graph segments" with "alternate alleles from the
multisample VCF."

## Interpretation Caveats

- Target-private pangenome features are hypotheses, not causal claims.
- GFA segment features and VCF allele features are comparable through the Privy
  matrix model, but they represent different upstream data products.
- Missing samples, graph construction choices, variant normalization, and
  upstream filtering can change core/accessory/private counts.
- Publication figures should report the input representation and group
  definitions directly in the caption or methods.
