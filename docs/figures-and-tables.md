---
title: Figures and Tables
description: How to interpret Panex Privus pangenome outputs and adapt them for publications.
---

# Figures and Tables

Panex Privus writes analysis outputs as machine-readable tables, Markdown/HTML
reports, and publication-starting plots. The tables are the analysis record.
The plots and reports are interpretive summaries that help you explain
target-private signal, VCF/GFA agreement, read-level support, annotations, and
pangenome structure in a manuscript, report, or slide deck.

## Command Outputs At A Glance

| Command | Main tables | Main figures or reports | Research use |
|---------|-------------|-------------------------|--------------|
| `privy scan` | `hits.tsv`, `regions.tsv`, `sample_support.tsv`, `evidence.tsv`, `qc.tsv` | Input to `privy report` and `privy plot` | Candidate discovery and ranking |
| `privy pangenome` | `feature_summary.tsv`, `coverage_histogram.tsv`, `composition.tsv`, `growth_curves.tsv` | `pangenome_growth.png`, `pangenome_coverage.png`, `pangenome_composition.png` | Full, target, and off-target pangenome summaries |
| `privy compare` | `compare.tsv`, `compare_summary.tsv` | `compare_summary.png` through `privy plot` | VCF/GFA or run-to-run concordance |
| `privy report` | `summary.tsv`, `ranked_hits.tsv`, `strictness_summary.tsv` | `report.md`, `report.html` | Shareable candidate summaries |
| `privy plot` | Reads scan/compare TSVs | `locus_panel.png`, `strictness_bar.png`, `score_distribution.png`, `support_bar.png`, `compare_summary.png` | Diagnostic and presentation figures |
| `privy annotate` | `annotated_hits.tsv`, `annotation_summary.tsv` | Input to downstream figures/tables | Gene-context interpretation |
| `privy export` | `hits.bed/gff3`, `regions.bed/gff3` | Genome browser tracks | IGV, JBrowse, bedtools, and genome-browser workflows |

## Scan Tables

`privy scan` is the discovery step. It asks which alleles or graph segments are
present in the target group and absent from off-target genomes.

The most publication-relevant scan table is `hits.tsv`. Each row is one
candidate private locus. For a manuscript or supplement, useful columns include:

| Column | How to use it |
|--------|---------------|
| `locus_id` | Stable identifier for cross-referencing figures and supplements |
| `contig`, `start`, `end` | Genomic interval; coordinates are 0-based half-open |
| `variant_type` | Distinguishes SNP, indel, SV, or graph-region candidates |
| `allele_key` | Source-specific allele or graph-segment identifier |
| `target_support_n`, `offtarget_support_n` | Cohort support counts |
| `target_missing_n`, `offtarget_missing_n` | Missingness counts to report uncertainty |
| `strictness_class` | Missingness-aware confidence category |
| `final_score` | Ranking score for prioritization |

**Table caption example.** Ranked target-private candidate loci discovered by
Panex Privus. Candidate intervals are reported in 0-based half-open
coordinates. `target_support_n` and `offtarget_support_n` give the number of
samples supporting the allele or graph segment in each cohort, while
`target_missing_n` and `offtarget_missing_n` preserve missing or uninformative
data separately from biological absence. Candidates are sorted by `final_score`.

`regions.tsv` is better for locus clusters. Use it when nearby candidate hits
should be discussed as one candidate interval rather than many individual
alleles or graph segments.

**Table caption example.** Candidate private regions produced by merging nearby
passing loci. Each region records the number of constituent loci, dominant
strictness class, target consistency, off-target exclusion, and maximum or
representative ranking score.

## Compare Tables

`privy compare` reconciles two scan outputs, usually one VCF run and one GFA
run. The goal is not to force the two sources to agree perfectly, but to make
agreement, partial agreement, and source-specific discoveries explicit.

Key `compare.tsv` columns:

| Column | How to use it |
|--------|---------------|
| `locus_id_a`, `locus_id_b` | Candidate IDs from the two input scans |
| `source_a`, `source_b` | Labels such as `vcf` and `gfa` |
| `coordinate_overlap` | How much the two intervals overlap |
| `match_class` | `supported`, `partially_supported`, `contradicted`, `source_specific`, `uninformative`, or `missing_data` |
| `state_compatibility` | Whether the cohort-support states are compatible |
| `comparison_score` | Summary score for agreement |

**Table caption example.** Concordance between VCF-derived and graph-derived
target-private candidates. Candidate pairs were matched by reciprocal
coordinate overlap and classified by match class. Source-specific rows indicate
candidates detected in one input representation but not matched in the other.

`compare_summary.tsv` is useful in the main text because it reduces many locus
pairs to counts and percentages by match class.

## Plot Diagnostics

`privy plot` turns scan and compare outputs into quick figures. These are meant
as diagnostic and presentation-ready starting points; for formal publication,
export SVG or PDF and adjust labels, typography, or panel layout as needed.

Common plot outputs:

| Plot | Best use | Caption focus |
|------|----------|---------------|
| `locus_panel.png` | Show top-ranked candidates | Ranking by `final_score`; color by strictness class |
| `strictness_bar.png` | Summarize confidence classes | Counts of strict, relaxed, missing, and contradicted candidates |
| `score_distribution.png` | Inspect ranking distribution | Score range and class-specific score patterns |
| `support_bar.png` | Summarize evidence records | BAM/VCF/GFA support, absence, ambiguity, contradiction |
| `compare_summary.png` | Present VCF/GFA concordance | Match-class distribution between sources |

**Figure caption example.** Diagnostic summary of Panex Privus candidate loci.
The locus panel shows top-ranked candidates by final score, with colors
indicating strictness class. Strictness and score-distribution panels summarize
confidence and ranking behavior across all emitted loci.

## Report Tables

`privy report` packages scan outputs into collaborator-friendly Markdown or HTML.
Use it for lab notebooks, supplemental summaries, and review handoffs. The
report is not a replacement for the raw TSVs; it is a readable layer over them.

Publication-useful report outputs:

- `summary.tsv`: run-level counts and top-locus metadata
- `ranked_hits.tsv`: top candidates with explicit rank values
- `strictness_summary.tsv`: counts and percentages by strictness class
- `support_summary.tsv`: evidence class counts when evidence is supplied
- `contradiction_summary.tsv`: contradiction metrics from QC and compare inputs

**Table caption example.** Summary of Panex Privus scan results, including the
number of candidate loci, merged candidate regions, evaluated records, top-ranked
locus, and strictness-class distribution.

## Annotation And Export Tables

`privy annotate` connects candidate loci to a GFF3 annotation. Use
`annotated_hits.tsv` when you need to prioritize candidates by gene context.

**Table caption example.** Target-private candidate loci annotated against the
reference gene model. Annotation classes distinguish coding sequence, UTR,
exonic, intronic, and intergenic candidates, with gene identifiers reported for
overlapping genic features.

`privy export` writes BED or GFF3 files for genome browsers and interval tools.
These are usually not final manuscript tables, but they are useful for manual
inspection, figure panel preparation, and downstream intersection analyses.

**Track caption example.** BED track of Panex Privus candidate private loci,
scaled by final score and displayed against the reference genome annotation.

## Pangenome Tables And Plots

`privy pangenome` describes the full pangenome and the target/off-target
sub-pangenomes. It complements `privy scan`: the scan finds candidate
target-private loci, while pangenome analysis describes the feature space those
candidates come from.

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

![Pangenome growth example](assets/examples/pangenome-gfa/pangenome_growth.png)

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

![Pangenome coverage example](assets/examples/pangenome-gfa/pangenome_coverage.png)

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

![Pangenome composition example](assets/examples/pangenome-gfa/pangenome_composition.png)

**Figure caption example.** Core, accessory, private, and absent feature counts
for the full cohort, target sub-pangenome, and off-target sub-pangenome.
Categories are assigned independently within each group. A feature can therefore
be target-core and off-target-absent, which is the strongest pangenome-level
pattern for target-private signal.

Use this figure to communicate the structure of each group’s pangenome. It is
especially useful when the target group has a distinct sub-pangenome profile
relative to the background group.

## Reporting Pangenome Analyses In Methods

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

## General Interpretation Caveats

- Target-private loci and pangenome features are hypotheses, not causal claims.
- Score rank is a prioritization aid, not a probability of causality.
- GFA segment features and VCF allele features are comparable through the Privy
  matrix model, but they represent different upstream data products.
- Missing samples, graph construction choices, variant normalization, and
  upstream filtering can change core/accessory/private counts.
- Publication figures should report the command, input representation, cohort
  definitions, and important thresholds directly in the caption or methods.
