---
title: Architecture
description: Panex Privus software architecture, formal data model, and algorithmic methods.
math: true
mermaid: true
---

# Architecture

Panex Privus is a cohort-aware comparative genomics toolkit for discovering
target-private genomic signal: alleles or graph segments that are observed in a
target cohort and absent from an explicitly defined off-target cohort. This
page describes the software architecture and the main algorithms at a level
intended for manuscript review, supplement preparation, and reproducible method
inspection.

The architecture follows three scientific constraints.

1. Evidence must remain source-aware. VCF, GFA, BAM, GFF3, and derived tables
   have different evidentiary meanings and are not collapsed into an
   undocumented black box.
2. Missingness must remain distinct from absence. A missing genotype, missing
   graph path, or low-depth BAM observation is not treated as confirmed lack of
   signal.
3. Outputs must be auditable. Every discovery row is decomposed into support
   counts, strictness class, score components, and reproducibility metadata.

## Design Principles

Panex Privus is logic-centered rather than format-centered. VCF and GFA are
primary discovery backends, while BAM, GFF3, reports, plots, and dashboards are
supporting or interpretive layers. A VCF allele and a GFA graph segment are not
the same biological object, but both can be evaluated against the same
target/off-target inference question.

The implementation favors streaming and bounded intermediate state. The VCF
scan consumes one `pysam.VariantRecord` at a time. The GFA scanner builds a
compact scan index containing coordinate-tagged segments, sample traversal
bitmasks, and per-sample coordinate coverage intervals. Landscape analysis can
write row streams directly while retaining only the summaries needed for local
background and candidate-introgression blocks.

The package also separates discovery from interpretation. `privy scan`
identifies candidate target-private loci. `privy pangenome` describes the
feature space from which candidates arise. `privy landscape` describes
windowed VCF context. `privy compare` reconciles two result sets by coordinate
and state compatibility. Reporting, plotting, annotation, export, and
interactive dashboards consume existing outputs rather than silently rerunning
discovery.

## System Map

<figure class="method-diagram">
<div class="mermaid">
flowchart LR
  cohort["Cohort definition<br/>targets, off-targets, ignored"]
  config["Resolved configuration<br/>defaults + YAML + CLI"]
  vcf["Indexed VCF/BCF"]
  gfa["GFA graph<br/>segments, paths, walks"]
  bam["BAM support layer"]
  scan["privy scan<br/>discovery kernel"]
  hits["hits.tsv<br/>regions.tsv<br/>evidence.tsv<br/>run.json"]
  pangenome["privy pangenome<br/>feature matrix summaries"]
  landscape["privy landscape<br/>windowed VCF context"]
  compare["privy compare<br/>cross-source reconciliation"]
  annotate["annotate/export/plot/report/interactive"]

  cohort --> scan
  config --> scan
  vcf --> scan
  gfa --> scan
  bam --> scan
  scan --> hits
  vcf --> landscape
  vcf --> pangenome
  gfa --> pangenome
  hits --> compare
  hits --> annotate
  pangenome --> annotate
  landscape --> annotate
</div>
<figcaption>Figure 1. Panex Privus separates cohort resolution, primary
discovery, contextual analyses, cross-source reconciliation, and downstream
interpretation.</figcaption>
</figure>

## Command Responsibilities

The command surface is organized by analytical phase rather than by file
format. This keeps the user-facing workflow aligned with the scientific
question being asked.

| Command | Architectural role | Primary output contract |
|---------|--------------------|-------------------------|
| `privy scan` | Primary discovery from VCF alleles or GFA graph segments | ranked hits, merged regions, evidence rows, QC, run metadata |
| `privy index gfa` | Reusable scan index construction for large GFA files | SQLite-backed `.privy.gfaidx` sidecar |
| `privy pangenome` | Full, target, and off-target feature-space summaries | feature summary, coverage, composition, growth curves |
| `privy landscape` | Windowed VCF context and local similarity analysis | sample/window metrics, similarity, local background, candidate blocks |
| `privy compare` | Cross-source or cross-run reconciliation | match classes, comparison scores, diagnostics |
| `privy report` | Human-readable synthesis of existing outputs | ranked summaries, strictness summaries, Markdown/HTML |
| `privy plot` | Static publication-oriented figures from existing tables | scan, landscape, and pangenome plots |
| `privy interactive` | Self-contained review dashboards | focus browsers and run-level HTML dashboards |
| `privy annotate` | Gene-model intersection | annotated hit tables and summary counts |
| `privy export` | Downstream interval conversion | BED and GFF3 tracks |

## Formal Data Model

Let $T$ denote the target sample set, $O$ the off-target sample set, and
$A = T \cup O$ the active cohort after ignored samples have been removed.
All scan outputs use 0-based half-open genomic intervals $[s,e)$.

The internal domain model is intentionally small.

| Object | Meaning | Principal fields |
|--------|---------|------------------|
| `CohortDefinition` | Biological grouping used by a run | targets, off-targets, ignored samples |
| `Locus` | A site or interval under evaluation | `contig`, `start`, `end`, type, primary source, source IDs |
| `AllelePattern` | Cohort-level support and missingness summary | target/off-target support, totals, missing counts, strictness class |
| `EvidenceRecord` | One normalized evidence statement | source type, sample or group, evidence class, metric, provenance |
| `ScoredHit` | Ranked discovery output | discovery, support, penalty, final score, strictness |
| `FeatureMatrix` | Sparse pangenome matrix | feature records, samples, feature-to-sample presence sets |

For each candidate allele or graph segment $a$, the discovery kernels compute
support indicators $S_{ia}$ and missing indicators $M_{ia}$ for sample
$i$. The central counts are:

$$
n_T^+(a) = \sum_{i \in T} S_{ia}, \quad
n_O^+(a) = \sum_{i \in O} S_{ia}
$$

$$
n_T^m(a) = \sum_{i \in T} M_{ia}, \quad
n_O^m(a) = \sum_{i \in O} M_{ia}
$$

Called-sample support fractions are then:

$$
p_T(a) = \frac{n_T^+(a)}{|T| - n_T^m(a)}, \quad
p_O(a) = \frac{n_O^+(a)}{|O| - n_O^m(a)}
$$

when the denominators are nonzero. Missingness fractions are tracked
separately:

$$
m_T(a) = \frac{n_T^m(a)}{|T|}, \quad
m_O(a) = \frac{n_O^m(a)}{|O|}
$$

This separation is the core statistical safeguard in the package: $p_O=0$
means no called off-target carries the signal, while $m_O>0$ means off-target
absence is not fully observed.

## Discovery Algorithms

`privy scan` currently implements the `private_allele` discovery mode for VCF
alleles and GFA graph segments. CLI names for additional modes are reserved,
but `private_genotype` and `private_sv_state` are not implemented discovery
kernels yet.

### VCF Private-Allele Kernel

The VCF backend streams an indexed multisample VCF by contig or requested
region. Each ALT allele in a record is evaluated independently. A sample
supports ALT $a$ if its called genotype contains the allele index for $a$;
the sample is missing if the genotype is missing or uninformative.

<figure class="method-diagram">
<div class="mermaid">
flowchart TD
  start["Indexed VCF record"]
  filters["Apply FILTER, QUAL,<br/>multiallelic policy"]
  alts["Enumerate ALT alleles"]
  counts["Count target/off-target<br/>support and missingness"]
  strict["Classify strictness"]
  pass{"pattern_pass?"}
  hit["Create Locus + HitRecord"]
  score["Score and rank"]
  regions["Merge nearby loci into regions"]
  outputs["Write TSV + JSON outputs"]

  start --> filters --> alts --> counts --> strict --> pass
  pass -- yes --> hit --> score --> regions --> outputs
  pass -- no --> outputs
</div>
<figcaption>Figure 2. VCF scanning is streaming at the record level. Passing
alleles are accumulated only after the target-private decision so they can be
ranked and merged into candidate regions.</figcaption>
</figure>

With thresholds $\tau_T$ (`min_target_support`) and $\tau_O$
(`max_off_target_support`), an allele can pass when:

$$
p_T(a) \ge \tau_T \quad \mathrm{and} \quad p_O(a) \le \tau_O
$$

The defaults are $\tau_T=1.0$ and $\tau_O=0.0$, meaning every called target
must carry the allele and no called off-target may carry it. The implementation
also requires at least one called target; when all targets are missing the
pattern is reported but not emitted as a passing hit.

### Strictness Classification

Strictness classes encode whether the target-private pattern is fully observed
or depends on incomplete data.

| Class | Emitted as pass? | Interpretation |
|-------|------------------|----------------|
| `strict_complete` | yes | Targets support, off-targets are absent, and no required sample is missing |
| `strict_target_missing` | yes | Off-target exclusion holds, but one or more targets are missing |
| `strict_offtarget_missing` | yes | Target support holds, but one or more off-targets are missing |
| `strict_both_missing` | yes | The pattern is otherwise consistent, but both cohorts contain missing calls |
| `relaxed_threshold` | context dependent | The pattern passes configured relaxed missingness handling or fails strict support thresholds |
| `contradicted` | no | Off-target support exceeds the allowed threshold or target support is insufficient |

The decision kernel gives contradiction priority over missingness. If a called
off-target carries the allele above $\tau_O$, the allele is contradicted even
if other off-target samples are missing. If target support is below
$\tau_T$, the allele is not emitted as a passing hit.

<figure class="method-diagram">
<div class="mermaid">
flowchart TD
  counts["support and missingness counts"]
  off{"p_O > tau_O?"}
  targets{"any called target?"}
  tpass{"p_T >= tau_T?"}
  missing{"missing in target<br/>or off-target?"}
  relax{"missingness exceeds<br/>relaxed tolerance?"}
  contradicted["contradicted"]
  no_pass["non-passing target-missing<br/>or relaxed-threshold state"]
  complete["strict_complete"]
  strict_missing["strict_target_missing /<br/>strict_offtarget_missing /<br/>strict_both_missing"]
  relaxed["relaxed_threshold"]

  counts --> off
  off -- yes --> contradicted
  off -- no --> targets
  targets -- no --> no_pass
  targets -- yes --> tpass
  tpass -- no --> no_pass
  tpass -- yes --> missing
  missing -- no --> complete
  missing -- yes --> relax
  relax -- yes --> relaxed
  relax -- no --> strict_missing
</div>
<figcaption>Figure 3. The strictness state is a named part of the result rather
than a hidden quality flag.</figcaption>
</figure>

### GFA Private-Segment Kernel

The GFA backend evaluates coordinate-tagged graph segments as graph-native
features. Segment coordinates come from `SN`, `SO`, and `LN` tags on GFA
`S`-lines. Sample traversal is read from `P` paths and `W` walks. W-line
coordinates and segment coordinates are 0-based half-open, matching the
`Locus` convention.

For each segment $g$, the scan index stores a traversal bitmask
$B_g$ over samples and per-sample coordinate coverage intervals
$\mathcal{C}_i$. At segment interval $I_g=[s_g,e_g)$:

$$
S_{ig}=1 \quad \mathrm{if} \quad i \in B_g
$$

$$
M_{ig}=1 \quad \mathrm{if} \quad i \notin B_g
\quad \mathrm{and} \quad
\mathcal{C}_i \cap I_g = \emptyset
$$

If a sample has coverage at $I_g$ but does not traverse $g$, it is counted
as absent, not missing. This distinction allows GFA scans to distinguish
alternative graph paths from unobserved graph coverage.

<figure class="method-diagram">
<div class="mermaid">
flowchart LR
  sline["S-lines<br/>SN/SO/LN coordinates"]
  paths["P/W traversals"]
  index["GFA scan index<br/>segment masks + coverage intervals"]
  segment["Coordinate-tagged segment"]
  present["Present mask at segment interval"]
  counts["Support, absence,<br/>missingness counts"]
  strict["Strictness classifier"]
  outputs["hits.tsv + graph_segments.tsv"]

  sline --> index
  paths --> index
  index --> segment --> present --> counts --> strict --> outputs
</div>
<figcaption>Figure 4. GFA discovery uses graph traversal as the support state
and coordinate coverage as the missingness state.</figcaption>
</figure>

GFA hits are private graph-segment evidence, not VCF ALT alleles. The
companion `graph_segments.tsv` therefore reports same-segment traversal,
coordinate-covered absence, missingness, and segment length summaries.

## Evidence and Scoring

Panex Privus uses transparent additive scoring as a ranking aid. Scores are
not probabilities of causality and should be interpreted with the component
columns retained.

For a passing hit $h$:

$$
F_h = D_h + S_h - P_h
$$

where $D_h$ is the discovery score, $S_h$ is the secondary support score,
and $P_h$ is the penalty score. The current VCF/GFA discovery score is:

$$
D_h =
w_D \min\left(
2,\,
\frac{p_T(h) + (1-p_O(h))}{2}
+ 0.2\,\mathbf{1}_{\mathrm{strict\_complete}}
+ 0.1\min(Q/60,1)
\right)
$$

where $Q$ is the VCF QUAL value when available; GFA scans use $Q=0$. The
support score is the weighted mean of normalized actionable secondary evidence:

$$
S_h = w_S \cdot \frac{1}{K}\sum_{k=1}^{K} e_k
$$

with $S_h=0$ when no actionable secondary evidence exists. The penalty is:

$$
P_h =
\begin{cases}
w_P, & \mathrm{if\ contradicted} \\
w_P \min(1,\ 0.4m_T + 0.3m_O), & \mathrm{otherwise}
\end{cases}
$$

The defaults are stored in the resolved configuration and written to
`run.json`. The score columns in `hits.tsv` are therefore reproducible from
the output table and metadata.

### BAM Support Layer

BAM evidence is queried only at previously discovered VCF loci. For SNPs,
Panex Privus counts reference, alternate, and other bases in each BAM pileup.
The ALT allele fraction is:

$$
\mathrm{AF}_{\mathrm{BAM}} = \frac{n_{\mathrm{ALT}}}{n_{\mathrm{REF}} +
n_{\mathrm{ALT}} + n_{\mathrm{other}}}
$$

Depth below `min_depth` is `uninformative`. For target samples, sufficient ALT
count and allele fraction are `support`; insufficient ALT evidence at adequate
depth is `ambiguous`. For off-target samples, sufficient ALT evidence is
`contradiction`; adequate depth without ALT evidence is `absence`. Non-SNP
loci currently record depth but are treated as uninformative for allele support
because simple pileup counts do not resolve indel, structural, or graph states
with sufficient specificity.

## Region Construction

Single markers are often too granular for biological interpretation, so
passing loci are merged into candidate regions after scoring. Loci are sorted
by `(contig, start)`. Adjacent loci $i$ and $j$ are merged if they are on
the same contig and:

$$
s_j - e_i \le d
$$

where $d$ is `merge_distance`. If `same_variant_class_only` is enabled, both
loci must also share the same `LocusType`. Region rows summarize the number of
constituent loci, variant-type composition, dominant strictness class, target
consistency, off-target exclusion, and mean final score.

## Pangenome Feature Architecture

`privy pangenome` converts VCF or GFA inputs into the same sparse feature
matrix:

$$
X_{fs} =
\begin{cases}
1, & \mathrm{feature\ } f \mathrm{\ is\ present\ in\ sample\ } s \\
0, & \mathrm{otherwise}
\end{cases}
$$

GFA features are graph segments. VCF features are individual ALT alleles.
For a group $G \in \lbrace A,T,O \rbrace$, feature coverage is:

$$
c_G(f) = \sum_{s \in G} X_{fs}
$$

Group categories are assigned as:

| Category | Rule |
|----------|------|
| `absent` | $c_G(f)=0$ |
| `private` | $c_G(f)=1$ |
| `accessory` | $1<c_G(f)<\lvert G \rvert$ |
| `core` | $c_G(f)=\lvert G \rvert$ |

The `target_private` flag is separate from the within-group `private`
category. It is true when $c_T(f)>0$ and $c_O(f)=0$.

Pangenome growth curves permute sample order. For permutation $\pi$, the
observed pangenome size after $n$ samples is:

$$
G_{\pi}(n) =
\left|\left\{f: \sum_{r=1}^{n} X_{f,\pi(r)} > 0 \right\}\right|
$$

The implementation also reports base-pair-weighted growth and singleton
features. These curves are descriptive finite-panel summaries, not model-based
extrapolations to unsampled diversity.

## Landscape Algorithms

`privy landscape` is a windowed VCF context engine. It does not alter scan
hits. Its purpose is to describe missingness, non-reference burden,
private/rare ALT burden, local sample similarity, and candidate donor-like
blocks under the same cohort definition used by discovery.

Windows can be fixed-record windows or fixed-base-pair windows. Record windows
stabilize the number of variants per window across heterogeneous variant
density; base-pair windows are easier to interpret on chromosome-scale axes.

For window $W$ and sample $i$, per-sample rates include:

$$
\mathrm{missing\_rate}_{iW} =
\frac{\#\{\mathrm{records\ in\ } W \mathrm{\ missing\ in\ } i\}}{|W|}
$$

$$
\mathrm{nonref\_rate}_{iW} =
\frac{\#\{\mathrm{called\ records\ in\ } W \mathrm{\ where\ } i
\mathrm{\ carries\ any\ ALT}\}}
{\#\{\mathrm{called\ records\ in\ } W \mathrm{\ for\ } i\}}
$$

For active sample pair $(i,j)$, local genotype similarity is the genotype
match fraction over records where both samples are called:

$$
\sigma_{ijW} =
\frac{1}{N_{ijW}}
\sum_{r \in W}
\mathbf{1}\{g_{ir}=g_{jr}\ \mathrm{and\ both\ are\ called}\}
$$

where $N_{ijW}$ is the number of compared records. Genotypes are normalized
before comparison, so allele ordering does not change the result. The nearest
local background for sample $i$ is:

$$
b(i,W) = \arg\max_{j \in A,\ j \ne i} \sigma_{ijW}
$$

with deterministic lexical tie-breaking.

<figure class="method-diagram">
<div class="mermaid">
flowchart TD
  vcf["Stream VCF records"]
  filters["Apply site filters<br/>PASS, QUAL, type, missingness, ALT frequency"]
  windows["Build record or bp windows"]
  metrics["Compute sample/window metrics"]
  sim["Compute pairwise local similarity"]
  nearest["Assign nearest background"]
  blocks["Merge adjacent background blocks"]
  intro["Filter target windows for<br/>candidate introgression blocks"]
  outputs["Write window, similarity,<br/>block, PCA, metadata tables"]

  vcf --> filters --> windows --> metrics --> sim --> nearest --> blocks --> outputs
  nearest --> intro --> outputs
</div>
<figcaption>Figure 5. Landscape analysis is a contextual, windowed analysis
that remains separate from target-private locus discovery.</figcaption>
</figure>

### Local Background and Candidate Introgression

Local background blocks merge adjacent windows for a sample when the nearest
background sample and nearest-background role remain constant and the
similarity exceeds `min_background_similarity`. A low-similarity window is
assigned to `unassigned`.

Candidate introgression blocks are stricter, target-only summaries. A target
sample window is eligible when:

$$
b(i,W) \in O
$$

$$
\sigma_{i,b(i,W),W} \ge \tau_{\mathrm{intro}}
$$

$$
\mathrm{missing\_rate}_{iW} \le \mu_{\max}
$$

and, when a nearest target similarity is available:

$$
\sigma_{i,b(i,W),W} - \max_{t \in T,\,t \ne i}\sigma_{itW}
\ge \delta_{\min}
$$

Adjacent eligible windows are merged when they have the same candidate donor
and consecutive window indices. Blocks shorter than `min_introgression_windows`
are discarded. These blocks are exploratory donor-like intervals; they are not
a formal local ancestry model, recombination map, or causal test.

### Local PCA

When `--local-pca` is enabled, each window's similarity matrix is transformed
into a distance matrix $D_{ijW}=1-\sigma_{ijW}$. Classical metric embedding is
then applied using:

$$
B_W = -\frac{1}{2}J D_W^2 J
$$

where $J=I-\frac{1}{n}\mathbf{1}\mathbf{1}^{\top}$. The leading positive
eigenvectors define two exploratory coordinates. Axis sign and scale should be
interpreted locally, not as globally aligned ancestry axes.

## Cross-Source Comparison

`privy compare` reconciles two `hits.tsv` files without reopening the original
VCF or GFA. Source B is indexed by canonical contig name. For each source A
locus, candidate source B loci are selected by coordinate proximity and tested
for overlap.

For intervals $A=[s_A,e_A)$ and $B=[s_B,e_B)$:

$$
I = \max(0,\ \min(e_A,e_B)-\max(s_A,s_B))
$$

$$
\rho_{\mathrm{reciprocal}} =
\frac{I}{\max(e_A,e_B)-\min(s_A,s_B)}
$$

$$
\rho_{\mathrm{containment}} =
\min\left(1,\ \max\left(\frac{I}{|A|}, \frac{I}{|B|}\right)\right)
$$

`overlap_mode` selects the overlap score used for matching. When no overlap
passes, a breakpoint-tolerance fallback can classify near but non-overlapping
loci as partially supported. State compatibility is based on strictness class:
contradicted states are incompatible, and optional strictness compatibility
requires both states to be in the same broad strict/relaxed category.

## Implementation Architecture

The package is organized around a CLI layer, reusable domain models, file IO,
analysis backends, and output-oriented modules.

```text
src/privy/
  cli/             Typer command definitions and CLI-to-config resolution
  core/            cohort, locus, strictness, evidence, scoring, intervals
  io/              VCF, BAM, GFA, GFF3, BED, TSV, JSON readers/writers
  backends/        command-level analysis engines
  compare/         reusable comparison helpers
  pangenome/       feature-matrix model and summaries
  landscape/       VCF window algorithms
  plot/            publication-oriented static figures
  report/          Markdown and HTML report generation
  interactive/     self-contained HTML dashboards
  utils/           logging, metrics, validation, parallel helpers
```

Configuration priority is:

1. Package defaults from Pydantic models.
2. YAML configuration supplied with `--config`.
3. CLI overrides from the subcommand invocation.

The resolved configuration and run metadata are written to JSON outputs so
reviewers can inspect sample sets, thresholds, scoring weights, input paths,
and command-specific options.

## Complexity and Scaling

The VCF scan is $O(R \cdot A_{\mathrm{ALT}} \cdot |A|)$ over retained VCF
records $R$, ALT alleles per record, and active samples. Memory use is
dominated by emitted hits and regions rather than input VCF size.

The GFA scan-index build is linear in GFA records plus cohort path/walk segment
references. Scan-time evaluation is linear in coordinate-tagged segments after
region and minimum-length filters. The SQLite sidecar index allows repeated
scans to stream indexed segments by contig rather than rebuilding traversal
state.

Landscape analysis is linear in retained variants for filtering and window
construction. Pairwise similarity within each window scales as
$O(|A|^2)$, which is the expected cost of a full local similarity matrix. The
record-window implementation uses rolling accumulators so overlapping windows
do not recompute every per-variant contribution from scratch.

## Validation Strategy

The test architecture mirrors the scientific risk points.

| Test layer | Main targets |
|------------|--------------|
| Unit tests | cohort validation, strictness classification, scoring, interval merging, GFA/GFF/VCF IO, pangenome summaries |
| Integration tests | VCF scan, GFA scan, BAM support, compare, landscape, pangenome, plot, report, annotate, export, interactive dashboards |
| Regression fixtures | missing data, multiallelic records, symbolic variants, off-target contradiction, coordinate normalization, region merging |

The most important invariant is that missingness must never be silently
converted into absence. The second most important invariant is coordinate
consistency: internal and tabular outputs use 0-based half-open intervals,
while GFF3 export converts to 1-based closed intervals.

## Non-Goals

Panex Privus is not a variant caller, assembler, pangenome graph constructor,
general genome browser, formal local ancestry tool, QTL mapper, or workflow
orchestrator. It assumes upstream tools have generated VCF, GFA, BAM, and GFF3
inputs and focuses on auditable target/off-target-aware interpretation.

## Reviewer's Reading Guide

For a methods review, the key claims to inspect are:

- The VCF and GFA discovery kernels use the same explicit support and
  missingness accounting.
- Strictness classes preserve uncertainty that would otherwise be hidden in a
  binary pass/fail result.
- Scores are decomposed and reproducible from output columns plus `run.json`.
- Pangenome and landscape modules are contextual analyses, not hidden changes
  to scan decisions.
- Candidate introgression blocks are deliberately framed as exploratory
  donor-like intervals, not formal local ancestry calls.
