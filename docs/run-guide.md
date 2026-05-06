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

## Command Shape

Global options come before the subcommand:

```bash
privy --config configs/privy.yaml --project-name soybean-protein scan ...
```

Common global options:

| Option | Description |
|--------|-------------|
| `--config PATH` | YAML configuration file |
| `--project-name TEXT` | Project name written into outputs |
| `--outdir PATH` | Default output directory |
| `--threads INT` | Reserved for supported parallel paths |
| `--log-level TEXT` | `debug`, `info`, `warning`, or `error` |
| `--quiet` | Reduce console output |
| `--version` | Show package version |

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
| `--off-targets TEXT [TEXT ...]` | Off-target sample names, for example `--off-targets O1 O2 O3` |
| `--cohort-file PATH` | YAML or TSV cohort definition |
| `--region TEXT` | Restrict to `contig:start-end` |
| `--contig TEXT` | Restrict to one contig |
| `--merge-distance INT` | Merge nearby hits into regions |

See [Figures and Tables](figures-and-tables.md#privy-scan) for example scan
tables and captions.

## GFA Scan

Use a GFA scan when your discovery question is based on a pangenome graph rather
than called variants. This asks the same biological question, but at the level of
graph segments: which segments are traversed by the target samples and absent
from off-target samples?

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

```bash
privy scan \
  --gfa pangenome.gfa.gz \
  --targets T1 T2 T3 \
  --off-targets O1 O2 O3 \
  --outdir results/
```

GFA outputs are written under `results/gfa/`.

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
| `--min-segment-length INT` | Minimum GFA segment length |
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
| `--bam-min-depth INT` | Minimum depth for informative evidence |
| `--bam-min-alt-count INT` | Minimum alternate-supporting read count |
| `--bam-min-alt-fraction FLOAT` | Minimum alternate allele fraction |

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
so the same tables and plots are available for graph and variant inputs.

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
sample in the GFA or VCF as off-target. You can also use list files:

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
  --targets-file targets.txt \
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
- `pangenome_growth.png`: full, target, and off-target growth curves
- `pangenome_coverage.png`: feature coverage distribution
- `pangenome_composition.png`: stacked composition summary
- `pangenome.json`: run metadata, resolved groups, and output list

See [Figures and Tables](figures-and-tables.md) for example output snippets,
figure captions, and guidance on using pangenome plots in research reports or
publications.

## Compare Existing Scan Outputs

If you ran VCF and GFA discovery separately, compare their hits tables directly:

```bash
privy compare \
  --hits-a results/vcf/hits.tsv \
  --hits-b results/gfa/hits.tsv \
  --outdir results/compare/
```

Compare outputs:

- `compare.tsv`: per-locus source agreement table
- `compare_summary.tsv`: match-class counts and summary statistics
- `compare.json`: compare run metadata

Key compare options:

| Option | Description |
|--------|-------------|
| `--hits-a PATH` | First `hits.tsv` |
| `--hits-b PATH` | Second `hits.tsv` |
| `--source-a TEXT` | Optional label for source A |
| `--source-b TEXT` | Optional label for source B |
| `--min-reciprocal-overlap FLOAT` | Minimum overlap for interval matching |
| `--breakpoint-tolerance-bp INT` | Gap tolerance for near misses |
| `--require-state-compatibility` | Require strictness compatibility |

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

Use `privy plot` to make quick diagnostic figures from scan and compare outputs.
These plots are meant to help you see ranking, score distributions, strictness
classes, BAM evidence, and VCF/GFA concordance at a glance.

```bash
privy plot \
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
