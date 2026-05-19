---
title: Run Guide
description: End-to-end guide for running Panex Privus scans and downstream commands.
---

# Run Guide

This guide takes you from prepared input files to ranked private loci, reports,
plots, annotations, and exported intervals. The basic story is the same whether
you begin with a VCF or a graph: define the group you care about, define the
comparison group, ask Panex Privus what is private to the target group, then add
supporting evidence and summaries around the candidates.

For example tables, figure titles, and publication-style captions, see
[Figures and Tables](figures-and-tables.md).

## Where the Inputs Come From

Panex Privus does not build the VCF, graph, or read alignments itself. It starts
after you have generated comparative genomics inputs with upstream tools.

- A multisample **VCF** usually comes from variant calling against a reference
  genome, or from a pangenome workflow that emits genotypes for all samples.
- A **GFA** pangenome graph usually comes from a graph-building workflow such as
  minigraph-cactus, PGGB, or another pangenome assembler.
- **BAM** files come from mapping reads for each sample. For VCF support, reads
  are usually mapped to the same reference genome used for the VCF. For
  graph-centered projects, reads may also be mapped to a pangenome graph and
  projected or summarized in ways that preserve sample-level evidence.

The common pattern is: build or obtain VCF/GFA discovery inputs first, then use
BAMs as read-level support for candidate loci discovered from the VCF.

### Check Sample Names Before Choosing Cohorts

Use the sample names printed by the input file when choosing `--targets` and
`--off-targets`. Sample names are case-sensitive.

For a VCF:

```bash
bcftools query -l variants.vcf.gz
```

For a plain-text GFA:

```bash
awk -F'\t' '$1=="W"{print $2} $1=="P"{split($2,a,"#"); print a[1]}' pangenome.gfa | sort -u
```

For a compressed GFA:

```bash
gzip -cd pangenome.gfa.gz | awk -F'\t' '$1=="W"{print $2} $1=="P"{split($2,a,"#"); print a[1]}' | sort -u
```

GFA `W` lines store the sample name directly. GFA `P` path names often look like
`SAMPLE#HAP#CONTIG`, so the command above keeps the sample name before the first
`#`. If you run VCF and GFA scans together, check both lists and use names that
match the samples in the inputs you expect each backend to evaluate.

## Define Cohorts

All cohort-aware commands use the same cohort syntax:

- `privy scan`
- `privy pangenome`
- `privy landscape`

The most direct form is to list samples after grouped flags:

```bash
--targets T1 T2 T3 --off-targets O1 O2 O3
```

You can also put one sample per line in role-specific text files:

```bash
--targets-file targets.txt --off-targets-file off_targets.txt
```

Use `--ignore-samples I1 I2` or `--ignore-samples-file ignored.txt` for samples
present in the input that should be excluded from both groups.

For reproducible runs, prefer a single cohort file. The canonical TSV format
has a header with `sample_id` and `cohort_role` columns. Valid roles are
`target`, `off_target`, and `ignored`.

```text
sample_id	cohort_role
T1	target
T2	target
T3	target
O1	off_target
O2	off_target
O3	off_target
LowCoverage1	ignored
```

Save that as `cohort.tsv` and pass it to any cohort-aware command:

```bash
privy scan --vcf variants.vcf.gz --cohort-file cohort.tsv --outdir results/
privy pangenome --gfa pangenome.gfa.gz --cohort-file cohort.tsv --outdir results/pangenome/
privy landscape --vcf variants.vcf.gz --cohort-file cohort.tsv --outdir results/landscape/
```

YAML cohort files are also accepted:

```yaml
targets:
  - T1
  - T2
  - T3
off_targets:
  - O1
  - O2
  - O3
ignored_samples:
  - LowCoverage1
```

`privy scan` requires at least one target and one off-target sample. `privy
pangenome` and `privy landscape` can infer off-targets from every non-target,
non-ignored sample in the input when `off_targets` are omitted. If you combine
forms, explicit role flags and role-specific files override that same role from
the cohort file.

## Command Shape

Panex Privus commands have two layers of options:

1. global options for `privy` itself
2. subcommand options for `scan`, `pangenome`, `landscape`, `compare`, and the
   other workflow steps

Global options must come before the subcommand. Subcommand-specific options
come after it:

```bash
privy --config configs/privy.yaml --project-name soybean-protein scan --vcf variants.vcf.gz ...
```

Common global options:

| Option | Description |
|--------|-------------|
| `--config PATH` | YAML configuration file |
| `--project-name TEXT` | Project name written into outputs (default = `privy_run`) |
| `--outdir PATH` | Default output directory (default = `.`) |
| `--threads INT` | Reserved for supported parallel paths (default = 1) |
| `--log-level TEXT` | `debug`, `info`, `warning`, or `error` (default = `info`) |
| `--quiet` | Reduce console output (default = false) |
| `--version` | Show package version |

### Using `--config`

`--config PATH` points to a YAML file with reusable project settings. It is most
useful when you want the same cohort, thresholds, scoring weights, and compare
settings to be recorded and reused across runs.

The config file is optional. Any section you omit falls back to package
defaults. Explicit command-line flags override values from the config file.

Minimal example:

```yaml
project_name: soybean_privy_scan

cohorts:
  targets: [T1, T2, T3]
  off_targets: [O1, O2, O3]
  ignored_samples: []

scan:
  min_target_support: 1.0
  max_off_target_support: 0.0
  merge_distance: 1000
  min_qual: 30
  pass_only: true

gfa:
  min_segment_length: 1

compare:
  overlap_mode: contained
  min_reciprocal_overlap: 0.5
```

Use it like this:

```bash
privy --config configs/privy.yaml scan \
  --vcf variants.vcf.gz \
  --outdir results/
```

In that command, the target/off-target samples can come from
`configs/privy.yaml`. You can still override them at run time:

```bash
privy --config configs/privy.yaml scan \
  --vcf variants.vcf.gz \
  --targets T4 T5 \
  --off-targets O4 O5 O6 \
  --outdir results/
```

`--config` and `--cohort-file` are related but not identical:

- `--config` is a broader YAML settings file. It can include cohorts and many
  analysis parameters.
- `--cohort-file` only defines sample roles: target, off-target, and ignored.
  Use it when you want the same cohort file for `scan`, `pangenome`, and
  `landscape`.
- For `privy scan`, role-specific CLI inputs such as `--targets` or
  `--targets-file` override `--cohort-file`, and `--cohort-file` overrides the
  `cohorts:` section of `--config`.

## VCF Scan

Start with a multisample VCF when your discovery question is based on called
variants. Panex Privus streams through the VCF, tests each alternate allele
against the target/off-target cohort definition, and writes ranked candidate
loci.

First, make sure the VCF is bgzip-compressed and indexed:

```bash
bgzip -c variants.vcf > variants.vcf.gz
tabix -p vcf variants.vcf.gz
```

Then run the scan:

```bash
privy scan \
  --vcf variants.vcf.gz \
  --targets T1 T2 T3 \
  --off-targets O1 O2 O3 \
  --outdir results/
```

List all target samples after `--targets`, and all comparison samples after
`--off-targets`. You do not need to repeat the flags for `privy scan`.

VCF outputs are written under `results/vcf/`.

Scan outputs:

- `hits.tsv`: ranked candidate private alleles
- `regions.tsv`: merged candidate regions
- `sample_support.tsv`: per-sample genotype/evidence table
- `evidence.tsv`: per-locus evidence records
- `qc.tsv`: scan metrics
- `run.json`: run metadata and resolved configuration

Key scan options:

| Option | Description |
|--------|-------------|
| `--vcf PATH` | Indexed multisample VCF |
| `--targets TEXT [TEXT ...]` | Target sample names, for example `--targets T1 T2 T3` |
| `--targets-file PATH` | One target sample per line |
| `--off-targets TEXT [TEXT ...]` | Off-target sample names, for example `--off-targets O1 O2 O3` |
| `--off-targets-file PATH` | One off-target sample per line |
| `--ignore-samples TEXT [TEXT ...]` | Samples to exclude from both groups |
| `--ignore-samples-file PATH` | One ignored sample per line |
| `--cohort-file PATH` | YAML or TSV cohort definition |
| `--region TEXT` | Restrict to `contig:start-end` |
| `--contig TEXT` | Restrict to one contig |
| `--merge-distance INT` | Merge nearby hits into regions (default = 0) |

See [Figures and Tables](figures-and-tables.md#privy-scan) for example scan
tables and captions.

## GFA Scan

Use a GFA scan when your discovery question is based on a pangenome graph rather
than called variants. This asks a graph-specific version of the same biological
question: which coordinate-tagged graph segments are traversed by the target
samples and not traversed by off-target samples?

GFA hits are private graph-node calls, not VCF ALT-allele calls. In a graph, an
off-target sample may take a different path through the same coordinate interval
instead of traversing the target segment. Review `graph_segments.tsv` to see
segment length, length class, same-segment traversal counts, coordinate-coverage
counts, and the graph-specific interpretation for each GFA hit.

For large GFA graphs, build a reusable Privy GFA index before scanning:

```bash
privy index gfa --gfa pangenome.gfa.gz
```

This writes `pangenome.gfa.gz.privy.gfaidx` beside the GFA by default. The first
indexing run can take time because Privy must stream the full graph and all
sample walks. Later `privy scan --gfa pangenome.gfa.gz ...` runs auto-detect
that sidecar index, validate that it still matches the GFA file, and skip the
slow GFA walk-parsing step. Use `--gfa-index PATH` if the index is stored
somewhere else.

If the index was built with an older development version, or after pulling a
Privy update that changes GFA indexing, rebuild it once:

```bash
privy index gfa --gfa pangenome.gfa.gz --force
```

```bash
privy scan \
  --gfa pangenome.gfa.gz \
  --targets T1 T2 T3 \
  --off-targets O1 O2 O3 \
  --outdir results/
```

GFA outputs are written under `results/gfa/`. In addition to the common scan
files, GFA scans write `graph_segments.tsv`, a companion table for interpreting
private graph-node evidence without pretending it is the same thing as a VCF
alternate allele.

The `--gfa` input may be a plain-text `.gfa` file or a gzip-compressed
`.gfa.gz` file, matching the compressed output commonly written by
minigraph-cactus.

For large pangenome graphs, `privy scan --gfa` builds a single-pass,
scan-specific streaming index instead of retaining full segment sequences,
links, walks, and paths in memory. It logs progress while indexing and then
reports per-contig scan progress. You can point `--gfa` directly at the
minigraph-cactus `.gfa.gz` output; you do not need to decompress it first.

GFA segments must have coordinate tags such as `SN:Z:chr1`, `SO:i:1000`, and
`LN:i:500`. Minigraph-cactus output usually includes these tags. Without them,
Panex Privus cannot place graph segments back onto genomic coordinates for
region merging, comparison, or annotation.

Key GFA options:

| Option | Description |
|--------|-------------|
| `--gfa PATH` | GFA graph file, `.gfa` or `.gfa.gz` |
| `--gfa-index PATH` | Optional prebuilt Privy GFA index; auto-detected as `<GFA>.privy.gfaidx` when present |
| `--min-segment-length INT` | Minimum GFA segment length (default = 1) |
| `--region TEXT` | Restrict to `contig:start-end` |
| `--contig TEXT` | Restrict to one contig |

See [Figures and Tables](figures-and-tables.md#privy-scan) for example scan
tables and captions.

## Add BAM Evidence

BAM support is the next layer after VCF discovery. The VCF scan asks, "Which
called alleles match the target-private pattern?" BAM support asks, "Do the
reads at those candidate loci agree with that pattern?"

BAM support does not discover new loci. It revisits loci discovered from the VCF
and adds per-sample depth, allele fraction, and evidence classes.

Before using BAMs, make sure each file is coordinate-sorted and indexed:

```bash
samtools sort sample.bam -o sample.sorted.bam
samtools index sample.sorted.bam
```

For short reads, this usually means aligning FASTQ files to the reference genome
used for the VCF with an aligner such as BWA-MEM2 or minimap2, then sorting and
indexing the resulting BAM. For long reads, minimap2 is commonly used against a
reference assembly; graph-based read mapping can also be useful, but the BAM
support layer currently expects coordinate-sorted BAM files that can be queried
at the VCF hit coordinates.

For a small number of BAMs, pass repeated `--bam` flags:

```bash
privy scan \
  --vcf variants.vcf.gz \
  --targets T1 T2 \
  --off-targets O1 O2 \
  --bam T1.sorted.bam \
  --bam O1.sorted.bam \
  --outdir results/
```

For real projects, a manifest is usually clearer:

```text
bam_path	sample_id
T1.sorted.bam	T1
T2.sorted.bam	T2
O1.sorted.bam	O1
O2.sorted.bam	O2
```

```bash
privy scan \
  --vcf variants.vcf.gz \
  --targets T1 T2 \
  --off-targets O1 O2 \
  --bam-manifest manifest.tsv \
  --outdir results/
```

Panex Privus matches BAM files to samples by the BAM header `@RG SM` sample tag.
If no `SM` tag is present, it falls back to the filename stem.

Key BAM options:

| Option | Description |
|--------|-------------|
| `--bam PATH` | BAM support file; repeat for multiple files |
| `--bam-manifest PATH` | TSV mapping BAM files to sample IDs |
| `--bam-min-depth INT` | Minimum depth for informative evidence (default = 8) |
| `--bam-min-alt-count INT` | Minimum alternate-supporting read count (default = 2) |
| `--bam-min-alt-fraction FLOAT` | Minimum alternate allele fraction (default = 0.2) |

## Run VCF and GFA Together

Many users will get both a GFA and a genotyped VCF from the same pangenome
workflow. You can provide both inputs in one command and let Panex Privus create
`vcf/`, `gfa/`, and `compare/` result directories.

```bash
privy scan \
  --vcf variants.vcf.gz \
  --gfa pangenome.gfa.gz \
  --targets T1 T2 T3 \
  --off-targets O1 O2 O3 \
  --outdir results/
```

This writes VCF discovery outputs to `results/vcf/`, graph discovery outputs to
`results/gfa/`, and reconciliation outputs to `results/compare/`.

## Analyze The Pangenome

Use `privy pangenome` when you want to summarize the whole input pangenome and
then ask how the target and off-target sub-pangenomes differ. This is separate
from `privy scan`: a scan looks for target-private candidates, while pangenome
analysis describes the feature space those candidates come from.

GFA segments and VCF alternate alleles both use the same feature-matrix model,
so the same tables are available for graph and variant inputs. Plots are
generated afterwards with `privy plot --plot-set pangenome`, or immediately
with `--plots` when you want a single command to make both tables and figures.

```bash
privy pangenome \
  --gfa pangenome.gfa.gz \
  --targets T1 T2 T3 \
  --outdir results/pangenome/
```

For a VCF, each alternate allele is treated as one pangenome feature:

```bash
privy pangenome \
  --vcf variants.vcf.gz \
  --targets T1 T2 T3 \
  --outdir results/pangenome/
```

If you provide targets but omit off-targets, Panex Privus treats every other
sample in the GFA or VCF as off-target. You can also use list files or a cohort
file:

```bash
privy pangenome \
  --gfa pangenome.gfa.gz \
  --targets-file targets.txt \
  --off-targets-file off_targets.txt \
  --permutations 100 \
  --outdir results/pangenome/
```

If you provide both a GFA and VCF, Panex Privus writes separate source
directories under the chosen output directory:

```bash
privy pangenome \
  --gfa pangenome.gfa.gz \
  --vcf variants.vcf.gz \
  --cohort-file cohort.tsv \
  --outdir results/pangenome/
```

This writes `results/pangenome/gfa/` and `results/pangenome/vcf/`.

Pangenome outputs:

- `feature_summary.tsv`: one row per feature with full, target, and off-target
  presence counts
- `coverage_histogram.tsv`: number of features and bp present in 0, 1, 2, ...
  samples for each group
- `composition.tsv`: core, accessory, private, and absent feature counts
- `growth_curves.tsv`: permutation-based pangenome growth data
- `pangenome.json`: run metadata, resolved groups, and output list

Render pangenome figures from an existing pangenome output directory:

```bash
privy plot \
  --plot-set pangenome \
  --input-dir results/pangenome/ \
  --output-format pdf
```

Pangenome plot outputs:

- `pangenome_growth.pdf`: full, target, and off-target growth curves
- `pangenome_coverage.pdf`: feature coverage distribution
- `pangenome_composition.pdf`: stacked composition summary

See [Figures and Tables](figures-and-tables.md) for example output snippets,
figure captions, and guidance on using pangenome plots in research reports or
publications.

## Explore VCF Landscapes

Use `privy landscape` when you want genome-wide window context around VCF
signal. It is not the private-locus discovery step and it is not a replacement
for specialized population-genetic tools. It is a Panex-native view of how
missingness, non-reference burden, private/rare ALT burden, and local sample
similarity change along chromosomes under the same target/off-target cohort
definition used by `privy scan`.

The default is fixed-record windows. This keeps the number of variants per
window stable even when variant density changes along the chromosome:

```bash
privy landscape \
  --vcf variants.vcf.gz \
  --targets T1 T2 T3 \
  --off-targets O1 O2 O3 \
  --window-records 200 \
  --step-records 50 \
  --outdir results/landscape/
```

Use base-pair windows when physical coordinates are more important than a
stable number of variants:

```bash
privy landscape \
  --vcf variants.vcf.gz \
  --targets T1 T2 T3 \
  --off-targets O1 O2 O3 \
  --window-bp 1000000 \
  --step-bp 250000 \
  --outdir results/landscape-bp/
```

To reproduce a classic filtered SNP-density workflow inside `privy landscape`,
filter to biallelic SNP records and use fixed base-pair windows:

```bash
privy landscape \
  --vcf variants.vcf.gz \
  --targets T1 T2 T3 \
  --off-targets O1 O2 O3 \
  --variant-type snp \
  --biallelic-only \
  --max-site-missing-rate 0.2 \
  --min-alt-carriers 2 \
  --window-bp 1000000 \
  --step-bp 250000 \
  --plots \
  --outdir results/landscape-snp-density/
```

The resulting `windows.tsv` contains `density_variants_per_kb`, and the plot
set writes `variant_density_profile.*` figures. Because the records were
filtered to SNPs first, those density values are SNPs per kilobase.

If you provide targets but omit off-targets, every other non-ignored sample in
the VCF becomes off-target. You can also use `--targets-file`,
`--off-targets-file`, `--ignore-samples I1 I2`, and
`--ignore-samples-file`, or use `--cohort-file cohort.tsv`, matching the scan
and pangenome command style.

Landscape outputs:

- `sample_windows.tsv`: per-sample metrics for every emitted window
- `windows.tsv`: target/off-target summary metrics for every window
- `filter_summary.tsv`: audit counts for record-level VCF filters
- `background_blocks.tsv`: adjacent windows merged by nearest local background
- `candidate_introgression_blocks.tsv`: target windows merged when the nearest
  local background is an off-target sample
- `similarity.tsv`: pairwise genotype similarity. By default this writes every
  window-by-pair row so chromosome-level similarity plots can be rendered
  later; use `--similarity-output summary` for compact genome-wide pair means
  or `--similarity-output none` to skip the table.
- `local_pca.tsv`: optional PCA-like local similarity coordinates, written
  when `--local-pca` is used
- `landscape.json`: run metadata, resolved samples, parameters, and outputs

Render landscape figures from an existing landscape output directory:

```bash
privy plot \
  --plot-set landscape \
  --input-dir results/landscape/ \
  --output-format pdf
```

Landscape plot outputs:

- `plots/landscape_plot_index.tsv`: index of rendered landscape plots
- `plots/variant_density_profile.<contig>.pdf`: window-level variant density;
  when `--variant-type snp` was used, this is SNP density
- `plots/missingness_heatmap.<contig>.pdf`: sample-by-window missingness
- `plots/private_burden_heatmap.<contig>.pdf`: sample-by-window private ALT burden
- `plots/local_background_map.<contig>.pdf`: nearest-background assignment
- `plots/similarity_cluster_map.<contig>.pdf`: chromosome-level sample
  similarity when per-window similarity rows are available

Key landscape options:

| Option | Description |
|--------|-------------|
| `--vcf PATH` | Multisample VCF or BCF |
| `--targets TEXT [TEXT ...]` | Target sample names |
| `--targets-file PATH` | One target sample per line |
| `--off-targets TEXT [TEXT ...]` | Off-target sample names |
| `--off-targets-file PATH` | One off-target sample per line |
| `--ignore-samples TEXT [TEXT ...]` | Samples to exclude from both groups |
| `--ignore-samples-file PATH` | One ignored sample per line |
| `--cohort-file PATH` | YAML or TSV cohort definition |
| `--window-records INT` | Number of VCF records per fixed-record window (default = 200) |
| `--step-records INT` | Record step between windows (default = 50) |
| `--window-bp INT` | Use base-pair windows of this size |
| `--step-bp INT` | Base-pair step (default = `--window-bp`) |
| `--variant-type TEXT` | Variant class to include: `all`, `snp`, `indel`, or `sv` (default = all) |
| `--biallelic-only` | Restrict to records with exactly one ALT allele |
| `--max-site-missing-rate FLOAT` | Maximum missing genotype fraction across active samples before windowing |
| `--require-active-alt` | Keep only records where at least one target/off-target sample carries an ALT allele |
| `--min-alt-carriers INT` | Minimum number of active samples carrying any ALT allele |
| `--min-alt-carrier-freq FLOAT` | Minimum active-sample ALT carrier frequency |
| `--max-alt-carrier-freq FLOAT` | Maximum active-sample ALT carrier frequency |
| `--rare-max-count INT` | Carrier-count threshold for rare ALT burden (default = 1) |
| `--rare-max-freq FLOAT` | Carrier-frequency threshold for rare ALT burden (default = 0.05) |
| `--min-background-similarity FLOAT` | Minimum nearest-sample similarity for assigning background blocks (default = 0.65) |
| `--min-introgression-similarity FLOAT` | Minimum target-to-off-target similarity for candidate introgression blocks (default = `--min-background-similarity`) |
| `--min-introgression-delta FLOAT` | Minimum advantage over the nearest target sample (default = 0.05) |
| `--max-introgression-missing-rate FLOAT` | Maximum target missingness allowed in candidate introgression windows (default = 0.5) |
| `--min-introgression-windows INT` | Minimum adjacent windows needed to emit a candidate block (default = 10) |
| `--similarity-output TEXT` | Pairwise similarity table mode: `full`, `summary`, or `none` (default = full) |
| `--vcf-engine TEXT` | VCF parser: `auto`, `pysam`, or `cyvcf2` (default = auto) |
| `--local-pca` / `--no-local-pca` | Write or skip optional local PCA coordinates (default = no local PCA) |
| `--plot-format TEXT` | Plot format for immediate `--plots`: `png`, `svg`, or `pdf` (default = png) |
| `--plots` / `--no-plots` | Write or skip landscape figures during analysis (default = no plots) |

For large runs, it is often better to separate table generation from plotting:

```bash
privy landscape \
  --vcf variants.vcf.gz \
  --targets T1 T2 T3 \
  --off-targets O1 O2 O3 \
  --window-records 200 \
  --step-records 50 \
  --outdir results/landscape/

privy plot \
  --plot-set landscape \
  --input-dir results/landscape/ \
  --output-format pdf
```

`privy plot --plot-set landscape` renders one set of plots per contig by
default. Use `--contig Gm10` or `--contigs Gm01,Gm02,Gm03` to limit output,
`--plot-scope genome` for the older whole-genome overview, or
`--plot-scope both` to write chromosome-level plots plus whole-genome summaries.

Use `pdf` or `svg` when you need vector text and axes for publication. Dense
heatmap panels may still be embedded as raster image layers inside the vector
file, which keeps very large windowed plots usable.

Interpret local background blocks as exploratory shared-genomic-background
segments. They are useful for seeing which genomes are locally similar, but
they are not by themselves a formal recombination-rate map. For controlled
crosses, founder panels, MAGIC populations, or pedigrees, the landscape outputs
can help choose regions and samples for more formal recombination or
founder-haplotype analyses.

Interpret `candidate_introgression_blocks.tsv` as a prioritized donor-like
local-background table. A block is emitted when a target sample is locally
closest to an off-target sample and passes the configured thresholds. Shared
ancestry, low recombination, selection, structural variation, missingness, and
VCF representation can produce similar patterns, so these rows should be
treated as candidate intervals for follow-up rather than definitive
introgression calls.

See [Figures and Tables](figures-and-tables.md#privy-landscape) for example
landscape table snippets, figure titles, and caption language.

## Compare Existing Scan Outputs

If you ran VCF and GFA discovery separately, compare their hits tables directly:

```bash
privy compare \
  --hits-a results/vcf/hits.tsv \
  --hits-b results/gfa/hits.tsv \
  --outdir results/compare/
```

By default, `privy compare` normalizes minigraph-cactus GFA contig names such
as `Sample#0#Gm01` to `Gm01` before matching them to VCF contigs. It also uses
`contained` overlap mode, which works well when short GFA graph segments fall
inside longer VCF intervals. Use `--overlap-mode reciprocal` when you want a
stricter interval comparison.

Compare outputs:

- `compare.tsv`: per-locus-pair source agreement table
- `compare_summary.tsv`: match-class counts and summary statistics
- `compare.json`: compare run metadata and diagnostics, including raw vs
  normalized contig overlap and candidate-match counts

Key compare options:

| Option | Description |
|--------|-------------|
| `--hits-a PATH` | First `hits.tsv` |
| `--hits-b PATH` | Second `hits.tsv` |
| `--source-a TEXT` | Optional label for source A (default = inferred from `locus_id`) |
| `--source-b TEXT` | Optional label for source B (default = inferred from `locus_id`) |
| `--overlap-mode TEXT` | Match mode: `contained`, `reciprocal`, or `any` (default = `contained`) |
| `--min-reciprocal-overlap FLOAT` | Minimum overlap score for `contained` or `reciprocal` matching (default = 0.5) |
| `--breakpoint-tolerance-bp INT` | Gap tolerance for near misses (default = 200) |
| `--require-state-compatibility` | Require strictness compatibility (default = false) |
| `--normalize-contigs` / `--no-normalize-contigs` | Normalize minigraph-cactus contig names before comparing (default = true) |

See [Figures and Tables](figures-and-tables.md#privy-compare) for example
compare tables, a compare-summary figure, and captions.

## Generate a Report

Use `privy report` when you want a compact, shareable summary of a scan. The
report command reads the TSV files produced by `privy scan` and turns them into
ranked tables, strictness summaries, QC summaries, and Markdown or HTML.

```bash
privy report \
  --hits results/vcf/hits.tsv \
  --regions results/vcf/regions.tsv \
  --evidence results/vcf/evidence.tsv \
  --qc results/vcf/qc.tsv \
  --compare results/compare/compare.tsv \
  --format both \
  --outdir results/report/
```

Report outputs:

- `summary.tsv`: run-level summary metrics
- `ranked_hits.tsv`: top hits with an explicit `rank` column
- `strictness_summary.tsv`: counts and percentages by strictness class
- `support_summary.tsv`: evidence counts by source and evidence class
- `contradiction_summary.tsv`: contradiction metrics from QC and compare inputs
- `report.md`: human-readable Markdown report
- `report.html`: browser-friendly HTML report when `--format html` or `both`

See [Figures and Tables](figures-and-tables.md#privy-report) for example report
tables and captions.

## Plot Diagnostics

Use `privy plot` to make figures from existing output tables. The default
`--plot-set scan` makes quick diagnostic figures from scan and compare outputs:
ranking, score distributions, strictness classes, BAM evidence, and VCF/GFA
concordance at a glance. The `landscape` and `pangenome` plot sets render
figures from existing result directories after those analyses finish.

```bash
privy plot \
  --plot-set scan \
  --hits results/vcf/hits.tsv \
  --evidence results/vcf/evidence.tsv \
  --compare results/compare/compare.tsv \
  --plot-type all \
  --outdir results/plots/
```

Plot outputs:

- `locus_panel.png`: ranked view of the top loci by `final_score`
- `strictness_bar.png`: strictness-class distribution
- `score_distribution.png`: `final_score` distribution by strictness class
- `support_bar.png`: evidence-class counts by source, enabled by `--evidence`
- `compare_summary.png`: VCF/GFA match-class distribution, enabled by `--compare`

Plot types:

- `locus_panel`
- `strictness_bar`
- `score_distribution`
- `support_bar`
- `compare_summary`
- `all`

Use `--output-format svg` or `--output-format pdf` if you want vector graphics
for editing or publication layouts.

See [Figures and Tables](figures-and-tables.md#privy-plot) for example plot
titles and captions.

## Interactive Focus Regions

`privy interactive` builds static, self-contained HTML dashboards for one or
more focus regions. The option is named `--focus`, but the object you are
choosing is a genomic region such as `Gm15:1-4000000`.

For novice users, start with a focus region around 4 Mbp or smaller. Larger
regions can be tractable when variant density and annotation density are low,
but the HTML file embeds its own data and JavaScript, so very large regions can
become slow to open or awkward to browse. If you are unsure, split a broader
locus into adjacent focus regions and compare the generated dashboards.

By default, point `--vcf` at a multisample VCF/BCF and Privy will extract a
focus-region sites table before rendering the dashboard. The extraction writes
`focus_<contig>_<start>_<end>.sites.tsv` beside the HTML so the browser remains
auditable and can be rebuilt later.

```bash
privy interactive \
  --focus Gm15:1-4000000 \
  --vcf cohort.vcf.gz \
  --gff3 Wm82.gene_exons.gff3.gz \
  --functional-tsv Wm82.functional_annotations.tsv \
  --samples Harosoy Harosoy-sharp Kingawa \
  --track-gff RepeatMasker=Wm82.repeats.gff3.gz \
  --track-gff SSR=Wm82.ssr_markers.gff3 \
  --sample-abbrev HS=Harosoy-sharp \
  --keyword-group Trichome=trichome,epidermal,auxin,bhlh,microtubule,cell_wall \
  --outdir results/interactive/
```

`--samples` is interpreted as `OFFTARGET DERIVED DONOR` for the first focus
browser workflow. The dashboard uses those columns in the sites table to label
genotypes and to interpret target-private or donor-like patterns.

For large VCFs, use bgzip-compressed and indexed inputs (`.vcf.gz` plus `.tbi`
or `.csi`) so focus-region extraction can use indexed random access. Plain VCFs
are acceptable for small tests.

Repeat `--focus` for multiple regions:

```bash
privy interactive \
  --focus Gm15:1-4000000 \
  --focus Gm12:2340000-2440000 \
  --vcf cohort.vcf.gz \
  --gff3 Wm82.gene_exons.gff3.gz \
  --samples Harosoy Harosoy-sharp Kingawa \
  --outdir results/interactive/
```

Privy writes one HTML file per focus region. For multi-region runs, it also
writes an `index.html` with links to each region dashboard and feature table.
Use `--sites-tsv` instead of `--vcf` when you already have a precomputed focal
genotype table and want to rebuild the same dashboard without touching the VCF.

Focus dashboard outputs:

- `focus_<contig>_<start>_<end>.html`: shareable interactive region browser
- `focus_<contig>_<start>_<end>.features.tsv`: ranked variant-supported feature table
- `focus_<contig>_<start>_<end>.sites.tsv`: extracted focal genotype table when `--vcf` is used
- `focus_<contig>_<start>_<end>.json`: reproducibility metadata for that region
- `interactive.json`: run-level input and output index
- `index.html`: multi-region dashboard index, written only when more than one focus region is supplied

The focus browser can display gene models, exons, CDS, computed introns,
strand-aware promoter windows, target-private SNPs, INDEL/complex records,
size/symbol-based SV-like records, and optional generic GFF3 tracks such as
repeat annotations or SSR markers. Candidate feature lists are variant-supported:
functional annotation can help rank or group features, but features without
compatible focal variation are not elevated.

### Functional Annotation TSV

`--functional-tsv` is optional, but useful when you want the browser to display
gene functions, boost candidate rankings, and create phenotype keyword groups.
The table must be tab-delimited and include at least one join key:

- `gene`: short gene/locus name matching the GFF3 `Name` value or browser gene label
- `gene_id`: full gene ID matching the GFF3 `ID` value

Recommended columns:

| Column | Meaning |
|--------|---------|
| `gene` | Short gene name, such as `Glyma.15G001600` |
| `gene_id` | Full reference annotation gene ID |
| `representative_predicted_function` | Main function label to show in the browser |
| `functional_category_keywords` | Semicolon-separated category terms |
| `screening_priority` | Optional priority label such as `high_screening_interest` |
| `screening_note` | Short rationale or curation note |

Minimal example:

```text
gene	gene_id	representative_predicted_function	functional_category_keywords	screening_priority	screening_note
Glyma.15G001600	Glyma.15G001600.Wm82.a6.v1	K07300 - Ca2+:H+ antiporter	trichome_or_epidermal_development; signal_transduction	high_screening_interest	target-private variants in gene body; calcium transport candidate
```

You can add phenotype-oriented browser groups with `--keyword-group`. Each group
searches feature type, label, gene name, functional category, predicted
function, screening priority, screening note, and extra GFF3 metadata. The group
still only lists features that already have target-private variation support.
Terms are case-insensitive substring matches. Use specific words or stems that
fit the phenotype question, and prefer several biologically related terms over
one very broad word. The exact group names and terms are written into the focus
dashboard JSON metadata for reproducibility.

```bash
privy interactive \
  --focus Gm15:1-4000000 \
  --vcf cohort.vcf.gz \
  --gff3 Wm82.gene_exons.gff3.gz \
  --functional-tsv Wm82.functional_annotations.tsv \
  --samples Harosoy Harosoy-sharp Kingawa \
  --keyword-group Trichome=trichome,epidermal,auxin,bhlh,microtubule,cell_wall \
  --keyword-group Insect-resistance=defense,insect,jasmonate,protease,transporter,phenylpropanoid,detox \
  --outdir results/interactive/
```

### Interactive Scan Dashboards

Use `privy interactive --scan` to build a self-contained dashboard from
existing `privy scan` outputs. This is a review dashboard for ranked private
loci, candidate regions, strictness classes, score distributions, QC metrics,
and optional VCF/GFA comparison summaries.

```bash
privy interactive \
  --scan results/scan/ \
  --max-hits 5000 \
  --max-regions 1000 \
  --outdir results/interactive/
```

`--scan` accepts either:

- a direct scan source directory containing `hits.tsv`, such as
  `results/scan/vcf/` or `results/scan/gfa/`
- a combined scan run directory containing source subdirectories such as
  `vcf/`, `gfa/`, and optional `compare/`

The dashboard aggregates the full `hits.tsv` files for charts and summary
counts, but embeds only the first `--max-hits` hit rows and `--max-regions`
region rows per source. This keeps the HTML shareable for very large scans
while preserving run-level counts in the dashboard metadata. Increase those
limits when you need deeper table browsing, or keep them modest when the HTML
will be emailed or opened on a laptop.

Scan dashboard outputs:

- `scan_dashboard.html`: shareable interactive scan dashboard
- `scan_dashboard.json`: reproducibility metadata and source summaries

### Interactive Landscape Dashboards

Use `privy interactive --landscape` to build a self-contained dashboard from
existing `privy landscape` outputs. This dashboard is designed for local
background and window-metric review: it includes a contig selector,
sample-by-window heatmap, window profile, candidate introgression block table,
pairwise similarity summaries, filtering provenance, and run metadata.

```bash
privy interactive \
  --landscape results/landscape/ \
  --max-windows 20000 \
  --max-sample-windows 80000 \
  --max-blocks 5000 \
  --outdir results/interactive/
```

`--landscape` expects a directory containing `windows.tsv` and
`sample_windows.tsv`. It also reads these optional files when present:
`candidate_introgression_blocks.tsv`, `background_blocks.tsv`,
`filter_summary.tsv`, `similarity.tsv`, and `landscape.json`.

The dashboard preserves full row counts in `landscape_dashboard.json`, but
embeds bounded rows in the HTML so large genome-wide landscapes remain
shareable. Increase `--max-windows`, `--max-sample-windows`, or `--max-blocks`
when you need more complete in-browser review.

Landscape dashboard outputs:

- `landscape_dashboard.html`: shareable interactive landscape dashboard
- `landscape_dashboard.json`: reproducibility metadata and source summaries

### Interactive Pangenome Dashboards

Use `privy interactive --pangenome` to build a self-contained dashboard from
existing `privy pangenome` outputs. This dashboard is designed for inventory
and composition review: it includes source-aware feature counts, target-private
and off-target-private feature counts, composition bars, coverage histograms,
growth curves, feature-type summaries, top contigs, and a searchable feature
table.

```bash
privy interactive \
  --pangenome results/pangenome/ \
  --max-features 10000 \
  --max-private-features 5000 \
  --outdir results/interactive/
```

`--pangenome` accepts either:

- a direct pangenome source directory containing `feature_summary.tsv`
- a combined pangenome run directory containing source subdirectories such as
  `vcf/` and `gfa/`

It reads `composition.tsv`, `coverage_histogram.tsv`, `growth_curves.tsv`, and
`pangenome.json` when present. The dashboard preserves full feature counts in
`pangenome_dashboard.json`, but embeds bounded feature tables so large
pangenomes remain shareable.

Pangenome dashboard outputs:

- `pangenome_dashboard.html`: shareable interactive pangenome dashboard
- `pangenome_dashboard.json`: reproducibility metadata and source summaries

## Annotate Hits

Use `privy annotate` when you want to connect candidate private loci to gene
models or other GFF3 features. This is usually the first biological
interpretation step after ranking: it tells you whether each hit is coding,
UTR, exonic, intronic, or intergenic.

```bash
privy annotate \
  --hits results/vcf/hits.tsv \
  --gff annotation.gff3.gz \
  --outdir results/annotated/
```

Annotation outputs:

- `annotated_hits.tsv`: all hit columns plus annotation class and gene context
- `annotation_summary.tsv`: counts and percentages by annotation class
- `annotate.json`: annotation run metadata

If your hit contig names differ from the GFF3 contig names, use a two-column
contig alias file with `--contig-alias`.

See [Figures and Tables](figures-and-tables.md#privy-annotate) for annotation
table caption guidance.

## Export Intervals

Use `privy export` when you want to move candidates into genome browsers or
interval-based tools. BED is convenient for IGV and bedtools; GFF3 is useful
when you want feature-style records with attributes.

```bash
privy export \
  --hits results/vcf/hits.tsv \
  --regions results/vcf/regions.tsv \
  --format gff3 \
  --outdir results/exported/
```

Export outputs:

- `hits.bed` or `hits.gff3`: one interval or feature per candidate hit
- `regions.bed` or `regions.gff3`: merged candidate regions
- `export.json`: export run metadata and written file paths

Use `--format bed` for BED output. Use `--kind hits`, `--kind regions`, or
`--kind both` to control which interval sets are exported.

See [Figures and Tables](figures-and-tables.md#privy-export) for export track
caption guidance.
