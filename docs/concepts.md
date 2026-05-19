---
title: Core Concepts
description: Definitions for Panex Privus cohorts, discovery logic, scores, modules, landscape metrics, and output interpretation.
---

# Core Concepts

Panex Privus is designed to make target-private genomics analyses explicit.
This page defines the biological question, the evidence model, and what each
module does internally. It is written for readers who want to understand the
tables and figures rather than treat the command line as a black box.

## Analysis Contract

Panex Privus asks a comparative question:

> Which alleles, graph segments, regions, or genome windows distinguish a
> target cohort from an off-target cohort?

The answer is evidence for an association pattern, not proof of causality. A
high-scoring locus may be biologically interesting because it is shared by the
target cohort and absent from the comparison cohort, but follow-up validation is
still needed before claiming mechanism.

The main design rules are:

- Cohort roles are explicit: every active sample is a target or an off-target.
- Missing data is not silently treated as absence.
- Coordinates are reported consistently as 0-based half-open intervals.
- Scores are decomposed into discovery, support, and penalty components.
- Context modules such as `privy pangenome` and `privy landscape` explain the
  data around discovery; they do not replace the discovery decision itself.
- Interactive dashboards make existing analysis outputs easier to browse and
  share; they are an interpretive layer over auditable TSV and JSON files.

## Command Map

| Module | Main question | Primary inputs | Main outputs |
|--------|---------------|----------------|--------------|
| `privy scan` | Which loci or graph segments match the target-private pattern? | VCF or GFA, optional BAM support | `hits.tsv`, `regions.tsv`, evidence, QC, run metadata |
| `privy index gfa` | Can a large GFA be pre-indexed for faster repeated scans? | GFA | `<GFA>.privy.gfaidx` or a user-chosen index path |
| `privy pangenome` | What is the full, target, and off-target feature composition? | GFA segments or VCF alternate alleles | Feature summary, coverage histogram, composition, growth curves |
| `privy landscape` | How do VCF metrics change across chromosomes in sliding windows? | Multisample VCF or BCF | Window tables, similarity tables, local background blocks, candidate introgression blocks |
| `privy compare` | Do two scan result sets support the same loci? | Two `hits.tsv` files | `compare.tsv`, summary table, metadata |
| `privy report` | How should scan outputs be summarized for collaborators? | Scan and optional compare outputs | Ranked tables, summaries, Markdown/HTML report |
| `privy plot` | What figures summarize an existing run? | Scan, landscape, or pangenome outputs | Scan diagnostics, landscape maps, and pangenome plots |
| `privy interactive` | How can results be reviewed as a shareable browser or dashboard? | Existing scan, landscape, or pangenome outputs; or focus-region VCF/GFF3 inputs | Self-contained HTML dashboards and JSON metadata |
| `privy annotate` | Which candidates overlap gene models? | `hits.tsv`, GFF3 | Annotated hits and annotation summary |
| `privy export` | How can candidates move into genome browsers or downstream tools? | `hits.tsv`, `regions.tsv` | BED or GFF3 tracks |

## Key Vocabulary

These terms appear throughout the documentation and output files.

| Term | Meaning in Panex Privus |
|------|-------------------------|
| Panex-native | Computed inside Panex Privus using its cohort definitions, coordinate conventions, explicit missingness model, and target/off-target framing |
| Cohort definition | The resolved target, off-target, and ignored sample lists used by a run |
| Target/off-target-aware | A metric is calculated or interpreted with sample roles preserved instead of treating all samples as one undifferentiated panel |
| Missingness | Genotype, graph-path, read, or source evidence that is absent or uninformative; missing data is reported separately from confirmed absence |
| Missing rate | In a landscape window, `missing_n / n_variants` for one sample |
| Non-reference burden | The amount of called VCF signal where a sample carries any ALT allele; in landscape output, `nonref_rate = nonref_n / called_n` |
| Burden | A count or rate of events in a region or window; it is descriptive, not a disease burden or causal effect estimate |
| ALT allele | A non-reference allele listed in a VCF record. Multiallelic records have more than one ALT allele, and Panex Privus can count each ALT separately |
| Rare ALT burden | Rare ALT allele events carried by a sample per called record, using the configured carrier-count or carrier-frequency thresholds |
| Private ALT burden | Cohort-private ALT allele events carried by a sample per called record |
| Local sample similarity | Pairwise genotype-match fraction between two samples within one landscape window, calculated only where both samples have called genotypes |
| Similarity output mode | The amount of pairwise similarity written to disk: compact genome summaries, full per-window pairs, or no similarity table |
| Local background | The sample most similar to a focal sample in a landscape window |
| Local background block | Adjacent landscape windows merged when the nearest local background stays the same and passes the similarity threshold |
| Local PCA coordinates | Optional two-axis embedding of each window's local similarity matrix for exploratory local-structure scans |
| Candidate introgression block | Adjacent target-sample windows where the nearest local background is an off-target sample and configured similarity, delta, missingness, and minimum-window filters pass |
| Sliding window | A fixed-record or fixed-base-pair interval moved along each contig to summarize local VCF patterns |
| Focus region | A user-selected genomic interval, such as `Gm15:1-4000000`, rendered by `privy interactive --focus` as one shareable region browser |
| Interactive dashboard | A self-contained HTML file that embeds bounded tables, summaries, JavaScript, and provenance for local review or collaborator sharing |

In plain language, the `privy landscape` sentence:

> Panex-native view of how missingness, non-reference burden, private/rare ALT
> burden, and local sample similarity change along chromosomes under the same
> target/off-target cohort definition used by `privy scan`

means:

> Use the same target and off-target groups from discovery, move across the VCF
> in windows, and report where calls are missing, where samples carry ALT
> alleles, where rare or cohort-private ALT alleles are concentrated, and which
> samples look most similar to each other locally.

## Cohorts and Sample Roles

A **target** sample belongs to the focal group. In a trait study, this might be
the trait-positive group. In a breeding or population comparison, it may be the
line, clade, panel, or treatment group of interest.

An **off-target** sample belongs to the comparison group. Off-target samples are
used to test whether a signal is absent outside the focal group.

An **ignored** sample is present in the input but excluded from the analysis.
Ignored samples do not contribute to support counts, absence counts, pangenome
categories, or landscape windows.

An **active** sample is any target or off-target sample after ignored samples
are removed. Active samples are the denominator for most cohort-aware summaries.

In `privy scan`, targets and off-targets define the discovery pattern. In
`privy pangenome` and `privy landscape`, if targets are provided but
off-targets are omitted, Panex Privus can infer the off-target group as every
non-target, non-ignored sample in the input.

`privy scan`, `privy pangenome`, and `privy landscape` all accept the same
cohort inputs: grouped `--targets`/`--off-targets` flags, role-specific sample
list files, or a YAML/TSV file passed with `--cohort-file`.

## What "Private" Means

The word **private** appears in several contexts. These meanings are related
but not identical.

| Term | Meaning |
|------|---------|
| Target-private allele or segment | Present in target samples and absent from off-target samples under the configured support thresholds |
| Off-target-private allele or segment | Present in off-target samples and absent from target samples |
| Pangenome `private` category | Present in exactly one sample within the group being summarized |
| Landscape private ALT event | An ALT allele carried by one cohort and absent from the opposite cohort within a window |

The pangenome `private` category is a singleton-within-group category. It does
not automatically mean target-private. For example, a feature present in one
off-target sample is `private` in the off-target pangenome and
`offtarget_private=True` in the feature summary.

## Coordinates and Intervals

Panex Privus output coordinates use the BED and pysam convention:

```text
0-based, half-open: [start, end)
```

For a SNP at VCF position `12345`, Panex Privus reports:

```text
start = 12344
end   = 12345
```

For GFF3 export, coordinates are converted back to the GFF3 convention:
1-based, closed intervals. BED export keeps 0-based half-open coordinates.

VCF-derived locus IDs use a `PPX` prefix. GFA-derived locus IDs use a `GPX`
prefix. Region IDs are produced when nearby passing loci are merged.

## Primary Inputs and Evidence Layers

VCF and GFA are primary discovery backends.

| Input | Role | Requirement |
|-------|------|-------------|
| VCF | Genotype-call-based allele discovery | bgzip-compressed `.vcf.gz` plus `.tbi` or `.csi` index |
| GFA | Graph-traversal-based segment discovery | GFA with coordinate tags for segments |

BAM is a support layer, not a discovery caller. BAM evidence is only queried at
loci that were already found by a scan.

GFF3 is an annotation layer. It explains where candidates fall relative to gene
models; it does not alter discovery scores.

If both `--vcf` and `--gfa` are supplied to one scan, Panex Privus writes
source-specific VCF and GFA result directories plus comparison outputs. You can
also run the two scans separately and compare them later with `privy compare`.

## Missingness and Strictness

Missingness is one of the central Panex Privus concepts. A missing genotype,
missing path, absent sample, or uninformative evidence record is not the same
thing as confirmed absence.

Every passing scan hit receives a `strictness_class`.

| Class | Meaning |
|-------|---------|
| `strict_complete` | All called targets support the signal, all called off-targets are absent, and no required samples are missing |
| `strict_target_missing` | Off-target exclusion holds, but one or more target samples are missing or uninformative |
| `strict_offtarget_missing` | Target support holds, but one or more off-target samples are missing or uninformative |
| `strict_both_missing` | The target-private pattern is consistent, but both cohorts contain missing data |
| `relaxed_threshold` | The result passes a relaxed user threshold or is downgraded by a missingness tolerance |
| `contradicted` | The private model fails, usually because off-target support exceeds the allowed threshold or target support is too low |

Support fractions are based on called samples when calls exist. Missing samples
are reported separately and can affect strictness and penalties.

A useful review order is:

1. Start with `strict_complete` hits.
2. Inspect missingness-aware strict classes when candidate biology is strong.
3. Treat `relaxed_threshold` as lower-confidence and threshold-dependent.
4. Treat `contradicted` as evidence against the target-private model.

## `privy scan`

`privy scan` is the primary discovery module. It produces candidate loci and
candidate regions that match a target-private pattern.

The currently implemented discovery mode is `private_allele`. The CLI reserves
names for `private_genotype` and `private_sv_state`, but those modes are not
implemented yet.

### VCF Discovery

In a VCF scan, Panex Privus streams records from an indexed multisample VCF. It
can restrict the scan by contig or region, require `FILTER=PASS`, apply a
minimum QUAL threshold, and include or skip multiallelic records.

For each VCF record, Panex Privus evaluates each alternate allele separately.
A sample supports an ALT allele when its called genotype contains that ALT
allele index. A sample is missing when the genotype is missing or uninformative.
The scan then counts:

- target samples supporting the ALT allele
- off-target samples supporting the ALT allele
- target samples missing at the record
- off-target samples missing at the record

The default private-allele rule requires all called targets to support the ALT
allele and no called off-targets to support it. These defaults are controlled
by `min_target_support=1.0` and `max_off_target_support=0.0`.

### GFA Discovery

In a GFA scan, Panex Privus evaluates graph segments instead of VCF alleles. A
segment can become a candidate when target paths or walks traverse it and
off-target paths or walks do not traverse that same segment.

This is private graph-node evidence. It should not be described as a VCF-style
ALT allele. In a graph traversal, an off-target sample may take an alternate
path through the same genomic interval, or it may have no informative walk
coverage there. `graph_segments.tsv` reports same-segment traversal and
coordinate-overlap coverage separately so those cases remain visible.

GFA discovery depends on segment coordinates. Segments without usable coordinate
tags cannot be placed on chromosome-scale output intervals and are skipped for
scan output. Sample path or walk coverage is used to distinguish:

- **traverses**: the sample path/walk uses the segment
- **absent**: the sample covers the locus but does not traverse the segment
- **missing**: the sample has no informative path/walk coverage for that locus

For large GFA files, `privy index gfa` can build a reusable sidecar index so
later GFA scans do not need to re-parse all graph walks.

### Candidate Regions

`regions.tsv` is built from passing loci. Loci are merged when they are on the
same contig and within the configured `merge_distance`. If
`same_variant_class_only` is enabled, only loci with the same variant class are
merged together.

Regions are summary intervals. They are useful for follow-up because a cluster
of target-private loci may be easier to validate than a single marker.

### Scoring

Panex Privus scoring is additive and transparent:

```text
final_score = discovery_score + support_score - penalty_score
```

`discovery_score` is based on the cohort pattern:

- target support fraction
- off-target exclusion fraction
- a specificity bonus for `strict_complete`
- optional contribution from VCF QUAL

`support_score` is based on secondary evidence such as BAM support. When no
secondary evidence is available, this component is zero.

`penalty_score` reduces confidence for missingness and contradiction. Target
missingness and off-target missingness contribute separately. Contradicted
patterns receive the strongest penalty.

The score is a ranking aid, not a probability that the locus is causal.

### BAM Support

BAM support asks whether aligned reads agree with VCF-discovered candidate
loci. It does not discover new loci.

For SNP loci, Panex Privus queries read pileups and records reference counts,
alternate counts, other-base counts, depth, and alternate allele fraction.
BAM evidence classes are:

| Evidence class | Meaning |
|----------------|---------|
| `support` | A target BAM has enough ALT-supporting reads |
| `absence` | An off-target BAM has enough depth and lacks ALT-supporting reads |
| `contradiction` | An off-target BAM carries the candidate ALT allele |
| `ambiguous` | A target BAM has depth but does not clearly support the ALT allele |
| `uninformative` | Depth is too low or the locus type cannot be evaluated from the BAM pileup |

For non-SNP loci, the current BAM layer reports depth but treats allele support
as uninformative because simple pileup counts do not resolve indel, structural,
or graph states well enough.

## `privy compare`

`privy compare` reconciles two scan result sets, usually VCF hits and GFA hits.
It does not re-open the original VCF or GFA. It compares the coordinates and
state summaries already written in two `hits.tsv` files.

The comparison engine:

1. Reads source A and source B hits.
2. Optionally normalizes contig names such as `SAMPLE#HAP#CONTIG` to `CONTIG`.
3. Finds same-contig coordinate candidates.
4. Applies the selected overlap rule: `contained`, `reciprocal`, or `any`.
5. Uses breakpoint distance as a fallback for near-but-not-overlapping loci.
6. Checks strictness compatibility when requested.
7. Emits match classes and comparison scores.

Match classes are:

| Match class | Meaning |
|-------------|---------|
| `supported` | The two sources overlap and their states are compatible |
| `partially_supported` | The sources are near or overlapping, but support is incomplete, weak, or state compatibility is imperfect |
| `contradicted` | One source contradicts the other source's target-private state |
| `source_specific` | A locus appears in only one source |
| `uninformative` | The comparison source exists but does not provide usable signal |
| `missing_data` | One or both sources have no data at the locus |

`comparison_score` is a concordance score from 0 to 1. Supported and partially
supported matches are scaled by overlap, while source-specific and contradicted
classes receive fixed lower scores.

## `privy pangenome`

`privy pangenome` describes the feature space behind discovery. It answers a
cohort-scale composition question:

> What features exist in the full cohort, target sub-pangenome, and off-target
> sub-pangenome?

The module converts inputs into a feature-by-sample presence matrix. The shared
feature model lets the same summaries work for different input types.

| Adapter | Feature definition |
|---------|--------------------|
| GFA | One graph segment |
| VCF | One alternate allele from one VCF record |

For each feature and each group, Panex Privus counts how many samples contain
the feature and assigns a category:

| Category | Meaning within the group being summarized |
|----------|-------------------------------------------|
| `absent` | Present in zero samples |
| `private` | Present in exactly one sample |
| `accessory` | Present in more than one sample but not all samples |
| `core` | Present in every sample |

The `target_private` flag means the feature is present in at least one target
and absent from every off-target. The `offtarget_private` flag is the reverse.

Pangenome coverage histograms count features and base pairs at each sample
coverage level. Growth curves permute sample order and report how the observed
pangenome expands as samples are added. These curves are descriptive summaries,
not model-based estimates of all unsampled diversity.

## `privy landscape`

`privy landscape` is a windowed VCF context module. It is separate from
`privy scan`: it explains how genome-wide VCF properties vary along
chromosomes, but it does not decide which individual loci are target-private
candidates.

Landscape analysis is Panex-native because it applies the same target and
off-target cohort framing used by discovery to sliding windows across the VCF.

### Window Modes

Landscape windows can be made in two ways.

| Mode | Meaning | Best use |
|------|---------|----------|
| Fixed-record windows | Each window contains a configured number of VCF records, advanced by a record step | Keeps variant count per window more stable across uneven SNP density |
| Base-pair windows | Each window spans a configured physical interval, advanced by a bp step | Easier to interpret on chromosome coordinate plots |

Coordinates in landscape outputs are 0-based half-open intervals. Window IDs
use the `LW` prefix.

### Per-Sample Window Metrics

`sample_windows.tsv` contains one row per active sample per window.

| Metric | Definition |
|--------|------------|
| `called_n` | Number of records in the window where the sample has a called genotype |
| `missing_n` | Number of records in the window where the sample genotype is missing |
| `missing_rate` | `missing_n / n_variants` |
| `het_n` | Number of called records with a heterozygous genotype |
| `het_rate` | `het_n / called_n` |
| `nonref_n` | Number of called records where the genotype contains any ALT allele |
| `nonref_rate` | `nonref_n / called_n` |
| `minor_genotype_n` | Number of called records where the sample's normalized genotype is among the least common genotype classes in that window record |
| `minor_genotype_rate` | `minor_genotype_n / called_n` |
| `rare_alt_n` | Number of rare ALT allele events carried by the sample |
| `rare_alt_rate` | `rare_alt_n / called_n` |
| `private_alt_n` | Number of cohort-private ALT allele events carried by the sample |
| `private_alt_rate` | `private_alt_n / called_n` |
| `median_call_freq` | Median frequency of the sample's genotype class across informative records in the window |
| `nearest_background` | Sample with the highest local genotype similarity to this sample |
| `nearest_similarity` | Genotype-match fraction to the nearest background sample |

ALT-count metrics are event counts. A multiallelic record can contribute more
than one ALT event, so `rare_alt_rate` and `private_alt_rate` are best read as
events per called record rather than guaranteed 0-to-1 probabilities.

### Rare ALT Burden

A rare ALT event is an ALT allele whose carrier count or carrier frequency is
small among active samples. The defaults can be changed with
`--rare-max-count` and `--rare-max-freq`.

Rare ALT burden is cohort-agnostic. A rare ALT can occur in targets,
off-targets, or both. It is useful for identifying windows enriched for
low-frequency variation.

### Private ALT Burden

Private ALT burden is the landscape metric behind the
`private_burden_heatmap.png` figure.

For each window and ALT allele:

- The ALT is **target-private** if at least one target carries it and no
  off-target carries it.
- The ALT is **off-target-private** if at least one off-target carries it and
  no target carries it.
- A target sample's `private_alt_n` increases when that target carries a
  target-private ALT.
- An off-target sample's `private_alt_n` increases when that off-target carries
  an off-target-private ALT.

`private_alt_rate` is:

```text
private_alt_n / called_n
```

This is not the same as the `privy scan` hit count. Landscape private ALT
burden is a window-level descriptive metric. It does not apply scan scoring,
strictness classes, region merging, BAM support, or candidate ranking.

### Pairwise Similarity and Nearest Background

For every pair of active samples in a window, Panex Privus compares records
where both samples have called genotypes. Genotypes are normalized before
comparison, so allele order does not affect the match.

```text
similarity = matching normalized genotypes / compared variants
```

`nearest_background` is the sample with the highest similarity to the focal
sample in that window. If there is a tie, the lexicographically earlier sample
name is chosen. The nearest background can be a target or an off-target; it is
a local similarity assignment, not a claim about ancestry by itself.

`--similarity-output` controls how much of the pairwise similarity matrix is
written to disk. `full` writes every sample-pair-by-window row and is the
default because chromosome-level landscape plots use those local similarities.
`summary` writes one genome-wide mean row per sample pair for leaner runs.
`none` skips the similarity table while still computing the internal local
similarity values needed for nearest-background and candidate-introgression
calls.

### Local PCA Coordinates

When `--local-pca` is enabled, Panex Privus writes `local_pca.tsv`. For each
window, it converts the pairwise local genotype similarity matrix into a
distance matrix and embeds that matrix into two PCA-like axes. This is closest
in spirit to local PCA or local structure scans: it helps show whether samples
separate differently in different chromosome intervals.

These coordinates are exploratory. Axis sign can flip between windows, and the
coordinates are not a formal ancestry or recombination model. Use them to find
regions where the local sample relationship changes, then follow up with a
method suited to the study design.

### Local Background Blocks

A local background block is a run of adjacent windows where a sample's nearest
background assignment stays the same and passes the minimum similarity
threshold. Blocks use the `LB` prefix.

If a window has no nearest sample or the similarity is below
`--min-background-similarity`, that window is assigned to `unassigned`.

Local background blocks are exploratory shared-genomic-background segments.
They should not automatically be interpreted as a formal recombination-rate
map. Formal genetic maps usually need a cross design, progeny, marker order,
and a recombination model. In Panex Privus, a local background block should be
read as:

> This sample looks locally most similar to this other sample across this
> chromosome interval.

For controlled crosses, MAGIC populations, founder panels, or pedigrees, the
landscape outputs can help choose regions and samples for more formal
recombination or founder-haplotype analyses.

### Candidate Introgression Blocks

A candidate introgression block is a stricter target-focused view of local
background. Panex Privus looks for adjacent windows where a target sample is
locally closest to an off-target sample. The block must pass the configured
minimum target-to-off-target similarity, optional similarity advantage over the
nearest target sample, maximum missingness, and minimum-window filters.

These rows are best read as:

> This target sample has a local donor-like background similar to this
> off-target sample across this chromosome interval.

They are candidate intervals for follow-up. Shared ancestral variation,
incomplete lineage sorting, low recombination, selection, structural variation,
missingness, and VCF representation can all produce donor-like local similarity
without recent introgression.

## `privy report`

`privy report` turns existing scan and compare outputs into collaborator-ready
summaries. It reads output tables; it does not re-run discovery.

The report module can write:

- run-level summary metrics
- top ranked hits
- strictness class counts
- support evidence summaries
- contradiction summaries
- optional candidate region summaries
- Markdown and optional HTML reports

The report is an interpretation layer. If a value in the report looks
surprising, the source of truth is still the underlying TSV and JSON output.

## `privy plot`

`privy plot` creates figures from existing output tables. The default
`--plot-set scan` creates diagnostic figures from scan and compare tables.
`--plot-set landscape` and `--plot-set pangenome` render figures from existing
result directories after those analyses finish.

Typical figures include:

- top-locus panels ranked by `final_score`
- strictness class bar plots
- score distributions
- evidence class summaries
- compare match-class summaries
- landscape heatmaps and local background maps
- pangenome growth, coverage, and composition plots

Plotting does not alter hit calls, scores, regions, window metrics, or
pangenome summaries.

## `privy interactive`

`privy interactive` builds shareable HTML dashboards from existing Privy outputs
or from one or more user-selected focus regions. It is designed for review,
collaboration, and supplementary exploration, not for changing the underlying
analysis.

The command currently has four dashboard modes:

- `--focus`: extract or read site-level genotype states for a genomic region
  and render an interactive genome/gene/variant browser.
- `--scan`: summarize existing scan outputs, including ranked hits, regions,
  score distributions, strictness classes, QC metrics, and optional VCF/GFA
  comparison summaries.
- `--landscape`: summarize existing landscape outputs with window profiles,
  sample-by-window metrics, local background assignments, candidate
  introgression blocks, filtering provenance, and metadata.
- `--pangenome`: summarize existing pangenome outputs with source-aware feature
  counts, composition, coverage histograms, growth curves, target-private
  features, and searchable feature tables.

Interactive dashboards do not rerun discovery. The scan, landscape, and
pangenome modes read existing result directories. Focus-region mode can extract
site genotypes from an indexed VCF/BCF, and it writes the extracted
`*.sites.tsv` table beside the dashboard so the browser state remains
auditable.

For large datasets, the HTML embeds bounded table rows while the companion JSON
metadata records full row counts, source files, and run parameters. This keeps
the dashboard portable enough to email or open locally while preserving the
analysis record in the original TSV/JSON outputs.

## `privy annotate`

`privy annotate` intersects scan hits with a GFF3 annotation. It classifies each
locus using this hierarchy:

```text
CDS -> UTR -> exonic -> intronic -> intergenic
```

The first overlapping gene is used for gene context. A contig alias file can be
used when the hit table and GFF3 use different contig names.

Annotation answers:

> Where does this candidate fall relative to known gene models?

It does not answer:

> Does this gene cause the phenotype?

## `privy export`

`privy export` converts scan hits and regions into downstream genome-tool
formats. It does not re-open the original inputs or change the analysis.

BED export keeps Panex Privus' 0-based half-open coordinates and scales
`final_score` into the BED 0-to-1000 score field. GFF3 export converts
coordinates to 1-based closed intervals and writes useful details as GFF3
attributes.

Use export files for genome browsers, interval joins, annotation workflows, and
other downstream tools.

## Configuration and Reproducibility

Configuration is resolved in this order:

1. Package defaults
2. YAML config passed with `--config`
3. CLI flags

The resolved configuration and run metadata are written to JSON outputs such as
`run.json`, `pangenome.json`, `landscape.json`, `compare.json`,
`annotate.json`, and `export.json`. These metadata files are intended to make
analyses reproducible and reviewable.

## Reading Results Without Treating Them as a Black Box

A practical review workflow is:

1. Start with `qc.tsv` and `run.json` to confirm what was scanned.
2. Inspect `hits.tsv` by `strictness_class` and `final_score`.
3. Use `sample_support.tsv` to see which samples drive each candidate.
4. Use `evidence.tsv` when BAM support was included.
5. Use `regions.tsv` to move from single markers to candidate intervals.
6. Use `privy compare` to check whether VCF and GFA evidence agree.
7. Use `privy pangenome` to understand cohort-wide feature composition.
8. Use `privy landscape` to inspect missingness, private ALT burden, and local
   background around chromosomes or candidate regions.
9. Use `privy interactive` dashboards for browsable review with collaborators.
10. Use annotation, export, reports, and plots to prepare follow-up
    interpretation.

The strongest Panex Privus results are usually not just high-scoring rows. They
are candidates where the cohort pattern, strictness class, sample support,
cross-source comparison, local genomic context, and biological annotation all
tell a coherent story.
