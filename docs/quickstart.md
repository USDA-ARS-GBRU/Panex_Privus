---
title: Quickstart
description: First commands for scanning, comparing, reporting, plotting, annotating, and exporting Panex Privus results.
---

# Quickstart

This walkthrough takes you from an input file to a ranked list of candidate
private regions. The story is the same whether you begin with a VCF or a graph:
define the group you care about, define the comparison group, ask Panex Privus
what is private to the target group, then add supporting evidence and summaries
around the candidates.

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

## VCF Scan

Start with a multisample VCF when your discovery question is based on called
variants. Panex Privus streams through the VCF, tests each alternate allele
against the target/off-target cohort definition, and writes ranked candidate
loci.

First, make sure the VCF is bgzip-compressed and indexed. The index lets Panex
Privus and downstream tools jump directly to genomic intervals.

```bash
bgzip -c variants.vcf > variants.vcf.gz
tabix -p vcf variants.vcf.gz
```

Next, choose the samples that define the biological contrast. In this example,
`T1`, `T2`, and `T3` are the group expected to share the signal, while `O1`,
`O2`, and `O3` are the samples where that signal should be absent.

```bash
privy scan \
  --vcf variants.vcf.gz \
  --targets T1 T2 T3 \
  --off-targets O1 O2 O3 \
  --outdir results/
```

The first file to inspect is `hits.tsv`. It is sorted by confidence, so the top
rows are the best candidates under the current filters and cohort definition.
The surrounding files explain how each hit was scored and how the run behaved.

- `hits.tsv`: ranked candidate private alleles
- `regions.tsv`: merged candidate regions
- `sample_support.tsv`: per-sample genotype/evidence table
- `evidence.tsv`: per-locus evidence records
- `qc.tsv`: scan metrics
- `run.json`: run metadata and resolved configuration

VCF scan outputs are written under `results/vcf/`.

## GFA Scan

Use a GFA scan when your discovery question is based on a pangenome graph rather
than called variants. This asks the same biological question, but at the level of
graph segments: which segments are traversed by the target samples and absent
from off-target samples?

Many users will get both a GFA and a genotyped VCF from the same pangenome
workflow. You can run them separately, or provide both inputs in one command and
let Panex Privus create `vcf/`, `gfa/`, and `compare/` result directories.

Run the graph scan with the same cohort logic:

```bash
privy scan \
  --gfa pangenome.gfa \
  --targets T1 T2 T3 \
  --off-targets O1 O2 O3 \
  --outdir results/
```

GFA segments must have coordinate tags such as `SN:Z:chr1`, `SO:i:1000`, and
`LN:i:500`. Minigraph-cactus output usually includes these tags. Without them,
Panex Privus cannot place graph segments back onto genomic coordinates for
region merging, comparison, or annotation.

GFA scan outputs are written under `results/gfa/`.

## Add BAM Evidence to a VCF Scan

BAM support is the next layer after VCF discovery. The VCF scan asks, "Which
called alleles match the target-private pattern?" BAM support asks, "Do the
reads at those candidate loci agree with that pattern?"

This is useful when you want read-level confidence around high-priority VCF
hits. BAM support does not discover new loci. It revisits loci discovered from
the VCF and adds per-sample depth, allele fraction, and evidence classes.

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

For a small number of BAMs, pass repeated `--bam` flags. Panex Privus matches
BAM files to samples by the BAM header `@RG SM` sample tag. If no `SM` tag is
present, it falls back to the filename stem. Repeated flags are most convenient
when your BAM headers already use the same sample names as the VCF.

```bash
privy scan \
  --vcf variants.vcf.gz \
  --targets T1 T2 \
  --off-targets O1 O2 \
  --bam T1.sorted.bam \
  --bam O1.sorted.bam \
  --outdir results/
```

For real projects, a manifest is usually clearer. It makes the sample mapping
explicit, avoids filename-stem surprises, and keeps long commands readable.

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

After the run, compare these files to the VCF-only output:

- `hits.tsv`: `support_score` and `final_score` reflect BAM evidence
- `evidence.tsv`: includes BAM evidence rows with `support`, `absence`,
  `contradiction`, `ambiguous`, or `uninformative` classes
- `sample_support.tsv`: includes depth and allele-fraction values where BAM
  pileup was informative
- `qc.tsv`: reports how many loci and BAM observations were evaluated

The evidence classes are intentionally conservative. Low-depth loci are marked
`uninformative` instead of being treated as absence, and off-target reads that
support the allele are marked as `contradiction`.

## Compare VCF and GFA Scans

If you ran both VCF and GFA discovery, compare the two result sets. Concordant
loci are strong follow-up candidates; source-specific loci are also useful,
because VCFs and graphs represent variation at different resolutions.

The simplest path is to provide both discovery inputs to `privy scan`:

```bash
privy scan \
  --vcf variants.vcf.gz \
  --gfa pangenome.gfa \
  --targets T1 T2 T3 \
  --off-targets O1 O2 O3 \
  --outdir results/
```

This writes VCF discovery outputs to `results/vcf/`, graph discovery outputs to
`results/gfa/`, and reconciliation outputs to `results/compare/`.

You can also compare scan directories explicitly:

```bash
privy compare \
  --hits-a results/vcf/hits.tsv \
  --hits-b results/gfa/hits.tsv \
  --outdir results/compare/
```

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

Use `--output-format svg` or `--output-format pdf` if you want vector graphics
for editing or publication layouts.

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
