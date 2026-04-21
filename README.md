# Panex Privus (`privy`)

**Find genomic variants shared within your target group and absent from everything else.**

Panex Privus is a comparative genomics toolkit for discovering *target-private alleles*: DNA variants that are present in a focal cohort of genomes and absent from a set of off-target genomes. It is designed for plant pangenome research and any study where you need to ask: *What is genetically unique to my group of interest?*

---

## Table of Contents

- [Background — what problem does this solve?](#background)
- [Requirements](#requirements)
- [Installation](#installation)
- [Your first scan — a step-by-step walkthrough](#your-first-scan)
- [Understanding the output files](#understanding-the-output-files)
- [Strictness classes — why missing data matters](#strictness-classes)
- [Configuration file](#configuration-file)
- [Command reference](#command-reference)
- [Current status](#current-status)
- [Troubleshooting](#troubleshooting)
- [Contributing](#contributing)
- [Citation](#citation)
- [License](#license)

---

## Background

### What is a target-private allele?

Imagine you are studying soybean breeding lines. You have five varieties that all share a desirable trait — say, high seed protein content — and ten other reference varieties that do not have that trait. A natural question is:

> **What DNA variants are present in all five high-protein varieties and absent from the ten reference varieties?**

Those variants are called *target-private alleles*. They may or may not *cause* the trait, but they are candidates worth investigating.

Panex Privus automates finding these variants from a VCF file or a pangenome graph (GFA). It handles:

- **Missingness explicitly** — some samples may have no genotype call at a position due to low coverage or quality filtering. Panex Privus never silently pretends missing data is an absence. Every result is labeled with a *strictness class* that tells you exactly what data was available.
- **Region merging** — nearby candidate variants are merged into genomic windows for downstream analysis.
- **Scoring and ranking** — results are scored by how confidently they fit the target-private pattern, so your highest-confidence candidates come first.

### What do I need to use this?

Privy uses a two-tier input model: a **primary discovery input** (required) and
optional **support layers** that add read-level validation on top.

#### Primary discovery inputs — choose one

VCF and GFA are two independent, standalone discovery backends. They answer the
same biological question from different representations of the data. Each backend
runs on its own; neither is a "support layer" for the other.

| Input | Role | Requirement |
|-------|------|-------------|
| **VCF** | Primary discovery — genotype-call-based allele detection | bgzip-compressed (`.vcf.gz`) + tabix index (`.tbi`) |
| **GFA** | Primary discovery — graph-traversal-based segment detection | Plain text file; no index needed |

Providing both in the same run does **not** combine them: when `--vcf` is given,
the VCF backend runs and the `--gfa` argument is ignored. To get independent
evidence from both representations, run them as separate `privy scan` calls and
compare the two result directories with `privy compare`.

**Recommended upstream tool**: [minigraph-cactus](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md)
(`cactus-pangenome`) produces both a `.gfa` and a genotype `.vcf.gz` in a single
workflow — exactly the files Panex Privus expects.

#### Support layers — optional

Support layers do not discover new loci. They query additional data at loci already
found by the primary backend and add evidence for or against each candidate.
The only current support layer is BAM.

| Layer | Flag | What it adds |
|-------|------|-------------|
| **BAM** | `--bam FILE` (repeat per BAM) or `--bam-manifest TSV` | Per-sample read depth and allele fraction at each VCF hit locus, classified as SUPPORT / ABSENCE / CONTRADICTION / AMBIGUOUS / UNINFORMATIVE |

When BAM files are provided, `privy scan` queries every hit locus via pysam pileup,
and updates `support_score` in `hits.tsv`, appends BAM-source rows to `evidence.tsv`,
and populates `depth` and `allele_fraction` in `sample_support.tsv`.

BAM files must be sorted by coordinate and indexed (`samtools index your.bam`).

---

## Requirements

### Python package

| Requirement | Version |
|-------------|---------|
| Python | 3.10 or higher |
| pysam | 0.22 or higher |
| Operating system | Linux or macOS (Windows via WSL) |

**New to bioinformatics?** If you are not sure whether you have Python 3.10+, open a terminal and run:

```bash
python --version
```

If the version shown is below 3.10, consider installing [Miniforge](https://github.com/conda-forge/miniforge) (a lightweight conda distribution) and creating a fresh environment.

### Upstream pangenome tool — minigraph-cactus

[minigraph-cactus](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md)
(`cactus-pangenome`) is the primary upstream tool for generating Panex Privus inputs.
Given a set of genome assemblies, it builds a pangenome graph and optionally calls
variants against a linear reference — producing both a GFA and a VCF in one workflow.

```bash
# Example: build a pangenome and output both GFA and VCF
cactus-pangenome ./js seqFile.txt \
    --outDir pangenome/ --outName my_pangenome \
    --reference RefSample \
    --vcf --giraffe --gfa
```

This gives you `pangenome/my_pangenome.gfa` and `pangenome/my_pangenome.vcf.gz` —
exactly the files Panex Privus expects.

Follow the [minigraph-cactus installation guide](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md#installation)
or install via conda:

```bash
conda install -c bioconda -c conda-forge cactus
```

---

### External command-line tools (required for VCF input preparation)

Panex Privus reads compressed, indexed VCF files. Preparing those files requires three standard bioinformatics tools that are **not** installed by `pip install .` and must be installed separately:

> **GFA users**: plain GFA text files do not need bgzip or tabix. If minigraph-cactus
> already produced a `.vcf.gz` with a `.tbi` index, you are ready to go. Skip this
> section if you are starting from a `.gfa` file only.

| Tool | What it does | Install via |
|------|-------------|-------------|
| **bgzip** | Block-compresses a VCF into the `.vcf.gz` format that tabix can index | Ships with `htslib` / `samtools` |
| **tabix** | Creates the `.tbi` position index that lets Panex Privus jump directly to any genomic region | Ships with `htslib` / `samtools` |
| **bcftools** | Swiss-army knife for VCF manipulation — used to list sample names, merge per-sample VCFs, filter records, etc. | Separate `bcftools` package |
| **samtools** | Sorts and indexes BAM files before the BAM support layer can query them (`samtools sort`, `samtools index`) | `samtools` package; ships alongside `htslib` |

The easiest way to get all three is through conda/bioconda:

```bash
# If you are using a conda environment (recommended)
conda install -c bioconda bcftools htslib

# Or install both in one shot with the samtools bundle
conda install -c bioconda samtools bcftools htslib
```

On Ubuntu/Debian Linux:
```bash
sudo apt-get install bcftools tabix
```

On macOS with Homebrew:
```bash
brew install bcftools htslib
```

Verify the tools are on your PATH:
```bash
bgzip  --version | head -1
tabix  --version | head -1
bcftools --version | head -1
```

---

## Installation

### Option 1: Standard install (recommended for users)

```bash
git clone https://github.com/USDA-ARS-GBRU/Panex_Privus.git
cd Panex_Privus
pip install .
```

### Option 2: Conda environment (recommended if pip gives errors)

```bash
conda create -n privy python=3.11
conda activate privy
conda install -c bioconda pysam   # installs pysam with compiled C extensions
pip install .
```

### Option 3: Developer install (for editing the source code)

```bash
git clone https://github.com/USDA-ARS-GBRU/Panex_Privus.git
cd Panex_Privus
pip install -e ".[dev]"
```

### Verify the installation

```bash
privy --version
privy --help
```

You should see version information and a list of subcommands.

> **Note on pysam**: `pysam` is a Python wrapper around the htslib C library. On most Linux and macOS systems `pip install .` compiles it automatically. If compilation fails, the safest fix is to install `pysam` via conda first (`conda install -c bioconda pysam`) and then install Panex Privus with `pip install .`.

---

## Your First Scan

This walkthrough takes you from an input file to a ranked list of candidate private regions.
Two paths are shown: **VCF** (variant-level) and **GFA** (graph-level).

Both inputs are most commonly produced by **minigraph-cactus** (`cactus-pangenome`).
A typical run produces `my_pangenome.gfa` and `my_pangenome.vcf.gz` in the same output
directory, giving you both paths immediately.

### VCF scan walkthrough

#### Step 1: Make sure your VCF is indexed

Panex Privus requires a bgzip-compressed, tabix-indexed VCF. If minigraph-cactus
already produced a `.vcf.gz` with a `.tbi` alongside it, skip to Step 2.
Otherwise, if your file is a plain `.vcf`, run:

```bash
bgzip -c your_variants.vcf > your_variants.vcf.gz
tabix -p vcf your_variants.vcf.gz
```

If you already have a `.vcf.gz` file but no `.tbi` file alongside it:

```bash
tabix -p vcf your_variants.vcf.gz
```

#### Step 2: List the sample names in your VCF

```bash
bcftools query -l your_variants.vcf.gz
```

You will see one sample name per line. These are the names you will pass to `--targets` and `--off-targets`.

#### Step 3: Run the scan

Suppose your VCF contains five target samples (`Benning`, `Harosoy`, `Clark`, `Williams`, `Essex`) and six off-target reference samples (`Jack`, `Lee`, `Minsoy`, `Richland`, `Dunfield`, `CNS`):

```bash
privy scan \
  --vcf your_variants.vcf.gz \
  --targets Benning Harosoy Clark Williams Essex \
  --off-targets Jack Lee Minsoy Richland Dunfield CNS \
  --outdir results/
```

This command will:

1. Stream through every variant in the VCF
2. For each variant, count how many of your target and off-target samples carry each allele
3. Classify the result using the strictness framework (see below)
4. Merge nearby passing variants into candidate regions
5. Score and rank all results
6. Write six output files to `results/`

For a typical plant genome (hundreds of millions of variants, tens of samples) this will finish in minutes.

#### Step 4: Explore your results

```bash
head -5 results/hits.tsv
```

The first rows are your highest-confidence candidates (rank 1 = best score).

---

### GFA scan walkthrough

If you have a pangenome graph from minigraph-cactus (`cactus-pangenome`), you can
run the discovery scan directly on the `.gfa` output file. No bgzip or tabix is needed.

The typical minigraph-cactus workflow produces a GFA where each sample's haplotype
appears as a W-line (walk) with the naming convention `SAMPLE#haplotype#contig` — for
example `Benning#1#chr1`. Panex Privus parses this automatically and uses the part
before the first `#` as the sample name.

#### Step 1: Check the sample names in your GFA

W-lines from minigraph-cactus follow the convention `SAMPLE#haplotype#contig`.
Panex Privus extracts the sample name automatically. For GFA1 P-lines, path names
are parsed the same way; plain names (no `#`) are also supported.

To preview the sample names Panex Privus will find, open the file and look at the
`W` or `P` lines:

```bash
grep "^W" your_graph.gfa | awk '{print $2}' | sort -u
grep "^P" your_graph.gfa | awk '{print $2}' | cut -d'#' -f1 | sort -u
```

#### Step 2: Run the scan

```bash
privy scan \
  --gfa your_graph.gfa \
  --targets Benning Harosoy Clark Williams Essex \
  --off-targets Jack Lee Minsoy Richland Dunfield CNS \
  --outdir results_gfa/
```

This command will:

1. Parse the entire GFA into memory and build traversal indices
2. For each graph segment that has reference-coordinate tags (SN/SO/LN), determine
   which samples traverse it and which have coverage at that position but traverse
   a different segment
3. Classify each segment using the same StrictnessClass framework as the VCF scan
4. Merge nearby passing segments into candidate regions
5. Score, rank, and write six output files to `results_gfa/`

#### Step 3: Explore your results

```bash
head -5 results_gfa/hits.tsv
```

GFA hits look the same as VCF hits: same columns, same strictness classes, same scoring.
The `variant_type` column will show `graph_region` and the `allele_key` column will
contain `contig:pos:SEG:segment_name` instead of `contig:pos:REF:ALT`.

#### Note on coordinate tags

Panex Privus can only place a segment on the output coordinate grid if that segment
carries `SN:Z:`, `SO:i:`, and `LN:i:` optional tags. **minigraph-cactus always writes
these tags** — if you built your GFA with `cactus-pangenome`, this will work
automatically. If your GFA came from a different tool and `hits.tsv` is empty, check
whether your segments have coordinate tags:

```bash
grep "^S" your_graph.gfa | head -3
```

You should see tags like `SN:Z:chr1  SO:i:1000  LN:i:500` after the sequence field.

---

### Adding BAM evidence to a VCF scan

BAM files provide a second, independent line of evidence: instead of trusting that
genotype calls in the VCF are correct, `privy scan` goes back to the raw reads and
asks *"Do we actually see the private allele in the reads?"*

#### Step 1: Make sure your BAMs are sorted and indexed

```bash
samtools sort -o T1.sorted.bam T1.bam
samtools index T1.sorted.bam
# Repeat for each sample BAM
```

If your BAMs are already coordinate-sorted (e.g., output from a standard aligner
pipeline), skip the sort step and just run `samtools index`.

#### Step 2: Add --bam flags to your scan

```bash
privy scan \
  --vcf your_variants.vcf.gz \
  --bam T1.sorted.bam \
  --bam T2.sorted.bam \
  --bam O1.sorted.bam \
  --bam O2.sorted.bam \
  --cohort-file cohort.yaml \
  --outdir results_with_bam/
```

For many samples, a manifest file is easier:

```bash
# manifest.tsv
# bam_path             sample_id
T1.sorted.bam          Benning
T2.sorted.bam          Harosoy
O1.sorted.bam          Jack
O2.sorted.bam          Lee
```

```bash
privy scan \
  --vcf your_variants.vcf.gz \
  --bam-manifest manifest.tsv \
  --cohort-file cohort.yaml \
  --outdir results_with_bam/
```

#### What changes in the outputs

- **`hits.tsv`**: `support_score` is now non-zero for loci where BAM reads confirm or
  contradict the VCF call.  `final_score` reflects this update.
- **`evidence.tsv`**: Gains one row per (locus, sample) pair for every BAM query, with
  `source_type=bam` and an `evidence_class` of `support`, `absence`, `contradiction`,
  `ambiguous`, or `uninformative`.
- **`sample_support.tsv`**: The `depth` and `allele_fraction` columns are now populated
  for BAM-covered samples (previously `NA`).

#### Evidence classes for BAM observations

| Class | When assigned |
|-------|--------------|
| `support` | Target sample: alt allele confirmed (depth ≥ min_depth, alt fraction ≥ threshold) |
| `ambiguous` | Target sample: adequate depth but alt allele not confirmed by BAM |
| `absence` | Off-target sample: alt allele absent with adequate depth |
| `contradiction` | Off-target sample: alt allele present — the private-allele model is undermined |
| `uninformative` | Any sample: depth below threshold, or locus is an indel (no per-allele counting) |

UNINFORMATIVE observations are excluded from the `support_score` mean so that
low-coverage samples do not penalise well-supported loci.

---

### Running a report

Once a scan is complete, generate a ranked summary and human-readable report:

```bash
privy report \
  --hits results/hits.tsv \
  --regions results/regions.tsv \
  --qc results/qc.tsv \
  --format both \
  --outdir report/
```

This writes:

| File | Contents |
|------|----------|
| `report/summary.tsv` | Run-level key/value summary |
| `report/ranked_hits.tsv` | Top 20 hits with explicit rank column |
| `report/strictness_summary.tsv` | Count and percentage per strictness class |
| `report/report.md` | Human-readable Markdown report |
| `report/report.html` | HTML version (requires `--format html` or `--format both`) |

Open `report/report.html` in a browser for a formatted view with tables and
navigation, or share `report/report.md` directly with collaborators.

To include evidence-level source support (requires `evidence.tsv`):

```bash
privy report \
  --hits results/hits.tsv \
  --regions results/regions.tsv \
  --evidence results/evidence.tsv \
  --qc results/qc.tsv \
  --format both \
  --title "Soybean protein scan" \
  --outdir report/
```

---

## Understanding the Output Files

After running `privy scan`, your output directory contains:

| File | What it contains |
|------|-----------------|
| `hits.tsv` | One row per candidate private allele, sorted by confidence score (rank 1 = highest). Start here. `support_score` is non-zero when BAM evidence was collected. |
| `regions.tsv` | Nearby hits merged into genomic intervals. Useful for downstream analysis. |
| `evidence.tsv` | One evidence record per hit per source type. With VCF only: one row per locus. With BAM added: one VCF row plus one row per (locus, sample) BAM observation, each with an `evidence_class`. |
| `sample_support.tsv` | Per-sample genotype at every hit locus. The `depth` and `allele_fraction` columns are populated when a BAM was provided for that sample; otherwise `NA`. |
| `qc.tsv` | Run-level quality control: how many records were evaluated, skipped, and passed. Review this to understand how your filters affected the analysis. |
| `run.json` | Complete record of every parameter used in this run. Keep this file to reproduce your results later. |

### Key columns in hits.tsv

| Column | Type | Description |
|--------|------|-------------|
| `locus_id` | text | Unique identifier for this locus (e.g., `PPX00000001`) |
| `contig` | text | Chromosome or scaffold name |
| `start` | integer | Start position (0-based, half-open — same convention as BED files) |
| `end` | integer | End position (0-based, half-open) |
| `allele_key` | text | VCF: `contig:pos:REF:ALT` (1-based). GFA: `contig:pos:SEG:segment_name` |
| `variant_type` | text | VCF: `snp`, `indel`, or `sv`. GFA: `graph_region` |
| `target_support_n` | integer | Number of target samples carrying this allele |
| `target_total_n` | integer | Total number of target samples |
| `offtarget_support_n` | integer | Number of off-target samples carrying this allele (should be 0) |
| `offtarget_total_n` | integer | Total number of off-target samples |
| `target_missing_n` | integer | Target samples with no genotype call at this position |
| `offtarget_missing_n` | integer | Off-target samples with no genotype call at this position |
| `strictness_class` | text | Confidence classification (see next section) |
| `discovery_score` | float | Score derived from the VCF/GFA cohort pattern alone |
| `support_score` | float | Score from secondary evidence (BAM); 0.0 when no BAM was provided |
| `penalty_score` | float | Deduction for missingness or contradiction |
| `final_score` | float | `discovery_score + support_score − penalty_score`. Higher is better. |

> **Coordinate note**: `start` and `end` use 0-based half-open coordinates (the BED/pysam convention). To convert to the 1-based position shown in your VCF, add 1 to `start`. For a SNP at VCF POS 12345, `start=12344`, `end=12345`.

---

## Strictness Classes

One of Panex Privus's most important features is that it never silently ignores missing genotype data. Every hit carries a `strictness_class` label that tells you exactly how complete the supporting data is.

| Class | Meaning | Emitted as a hit? |
|-------|---------|:-----------------:|
| `strict_complete` | All target samples carry the allele. All off-target samples are confidently absent. No missing data anywhere. **The gold standard.** | Yes |
| `strict_target_missing` | Pattern is consistent, but some target samples had no genotype call at this position. The called targets all support it. | Yes |
| `strict_offtarget_missing` | Pattern is consistent, but some off-target samples had no genotype call. The called off-targets are all absent. | Yes |
| `strict_both_missing` | Both target and off-target groups have some missing calls. The samples that were called fit the private-allele pattern. | Yes |
| `relaxed_threshold` | The allele passes your support thresholds, but missingness exceeded a tolerance you configured. Use with care. | Yes |
| `contradicted` | At least one off-target sample carries this allele — the private-allele model fails. | **No** |

### How to interpret this in practice

A `strict_complete` hit is the most trustworthy. A `strict_target_missing` hit is still worth examining, but you should check `sample_support.tsv` to understand which samples were missing and why before including it in a publication.

The `contradicted` class is never emitted as a hit, but it is counted in `qc.tsv` (`alleles_contradicted`). A high number of contradicted alleles is informative: it means there is substantial shared variation between your target and off-target groups in general.

---

## Configuration File

For reproducible analyses, define your parameters in a YAML file instead of typing them on the command line every time.

Create a file called `privy.yaml`:

```yaml
project_name: soybean_protein_scan

cohorts:
  targets:    [Benning, Harosoy, Clark, Williams, Essex]
  off_targets: [Jack, Lee, Minsoy, Richland, Dunfield, CNS]

scan:
  min_target_support:    1.0    # all called targets must carry the allele
  max_off_target_support: 0.0  # no off-targets allowed to carry the allele
  merge_distance:        1000  # merge hits within 1000 bp into a region
  pass_only:             true  # skip VCF records where FILTER != PASS
  min_qual:              30    # skip VCF records with QUAL < 30
  allow_multiallelic:    true  # evaluate multi-allelic records (one allele at a time)

scoring:
  discovery_weight: 1.0
  support_weight:   0.7
  penalty_weight:   0.8
```

Run with the config:

```bash
privy --config privy.yaml scan --vcf your_variants.vcf.gz --outdir results/
```

Any CLI flag you provide overrides the corresponding config value. Config values override the built-in defaults. This three-tier priority (defaults → config file → CLI flags) lets you store your standard settings in the file and override individual parameters on the command line when needed.

A fully annotated example config is in [`configs/privy.yaml`](configs/privy.yaml).

### Optional cohort file

If you prefer, you can store the cohort definition in a separate file and pass
it with `--cohort-file`.

YAML example:

```yaml
targets: [Benning, Harosoy, Clark, Williams, Essex]
off_targets: [Jack, Lee, Minsoy, Richland, Dunfield, CNS]
ignored_samples: [LowCoverageControl]
```

Run with:

```bash
privy scan --gfa your_graph.gfa --cohort-file cohort.yaml --outdir results_gfa/
```

TSV example:

```tsv
sample_id	cohort_role
Benning	target
Harosoy	target
Jack	off_target
Lee	off_target
LowCoverageControl	ignored
```

Run with:

```bash
privy scan --vcf your_variants.vcf.gz --cohort-file cohort.tsv --outdir results/
```

If you also pass `--targets`, `--off-targets`, or `--ignore-samples` on the
command line, those explicit CLI flags take precedence over the cohort file.

---

## Command Reference

### Global options

```
privy [--config PATH] [--outdir PATH] [--threads N] [--log-level LEVEL] [--quiet] COMMAND
```

| Option | Description |
|--------|-------------|
| `--config PATH` | Path to a YAML configuration file |
| `--outdir PATH` | Output directory (default: current directory) |
| `--threads N` | Number of threads (default: 1; multi-threading coming in a later release) |
| `--log-level` | Logging verbosity: `debug`, `info`, `warning`, `error` (default: `info`) |
| `--quiet` | Suppress all log output |

### privy scan

The primary discovery engine. Accepts a VCF **or** a GFA as the input — provide
whichever you have. When both are given, VCF is used as the primary source.

```bash
# VCF scan
privy scan --vcf PATH --targets SAMPLE [SAMPLE ...] --off-targets SAMPLE [SAMPLE ...] [OPTIONS]

# GFA scan
privy scan --gfa PATH --targets SAMPLE [SAMPLE ...] --off-targets SAMPLE [SAMPLE ...] [OPTIONS]
```

**Primary input options:**

| Option | Default | Description |
|--------|---------|-------------|
| `--vcf PATH` | — | Indexed multisample VCF (.vcf.gz + .tbi) |
| `--gfa PATH` | — | Pangenome graph file (GFA1 or GFA1.1, plain text) |
| `--cohort-file PATH` | none | YAML or TSV cohort definition file |
| `--targets TEXT` | required | Target sample names (repeat for each sample) |
| `--off-targets TEXT` | required | Off-target sample names (repeat for each sample) |
| `--outdir PATH` | `.` | Output directory |
| `--min-target-support FLOAT` | `1.0` | Minimum fraction of called targets that must carry the allele |
| `--max-off-target-support FLOAT` | `0.0` | Maximum fraction of called off-targets allowed to carry the allele |
| `--merge-distance INT` | `0` | Merge loci within this many bp into regions (0 = no merging) |
| `--pass-only / --no-pass-only` | `true` | (VCF) Require `FILTER=PASS` |
| `--min-qual FLOAT` | none | (VCF) Minimum VCF QUAL score |
| `--allow-multiallelic` | `true` | (VCF) Evaluate multiallelic records |
| `--min-segment-length INT` | `1` | (GFA) Skip segments shorter than this many bp |
| `--region TEXT` | none | Restrict scan to a region (format: `chr1:1000-2000`) |
| `--contig TEXT` | none | Restrict scan to a single contig |

**BAM support layer options** (all optional; BAM support activates when `--bam` or `--bam-manifest` is provided):

| Option | Default | Description |
|--------|---------|-------------|
| `--bam PATH` | none | Sorted, indexed BAM file. Repeat flag for multiple files. |
| `--bam-manifest PATH` | none | TSV with `bam_path` and `sample_id` columns; alternative to repeating `--bam` |
| `--bam-min-depth INT` | `8` | Minimum read depth to call evidence at a locus (below → UNINFORMATIVE) |
| `--bam-min-alt-count INT` | `2` | Minimum alt-supporting reads to call SUPPORT or CONTRADICTION |
| `--bam-min-alt-fraction FLOAT` | `0.2` | Minimum alt allele fraction to call SUPPORT or CONTRADICTION |
| `--bam-min-mapq INT` | `20` | Minimum mapping quality for BAM reads |
| `--bam-min-baseq INT` | `20` | Minimum base quality at the pileup position (SNPs only) |

Run `privy scan --help` for the full option list.

### privy report

Generate ranked summaries and a human-readable report from a previous scan.

```bash
privy report \
  --hits results/hits.tsv \
  --regions results/regions.tsv \
  --qc results/qc.tsv \
  --format both \
  --top-n 50 \
  --outdir report/
```

Key options:

| Option | Default | Description |
|--------|---------|-------------|
| `--hits PATH` | required | `hits.tsv` from `privy scan` |
| `--regions PATH` | none | `regions.tsv` from `privy scan` |
| `--evidence PATH` | none | `evidence.tsv` from `privy scan` |
| `--qc PATH` | none | `qc.tsv` from `privy scan` |
| `--run-json PATH` | none | `run.json` from `privy scan` |
| `--format TEXT` | `markdown` | `markdown`, `html`, or `both` |
| `--top-n INTEGER` | `20` | Number of top loci to include |
| `--title TEXT` | project name | Report title |
| `--outdir PATH` | `.` | Output directory |
| `--no-include-qc` | — | Omit the QC section |
| `--no-include-strictness` | — | Omit the strictness distribution |
| `--no-include-regions` | — | Omit the regions section |

Output files written to `--outdir`:

| File | Contents |
|------|----------|
| `summary.tsv` | Run-level key/value summary |
| `ranked_hits.tsv` | Top-N hits with explicit rank column |
| `strictness_summary.tsv` | Count and percentage per strictness class |
| `support_summary.tsv` | Evidence breakdown by source type (if `--evidence` provided) |
| `contradiction_summary.tsv` | Contradiction metrics from QC/compare |
| `report.md` | Human-readable Markdown report |
| `report.html` | HTML version (when `--format html` or `--format both`) |

### privy compare

`privy compare` reconciles two `privy scan` result sets — typically a VCF scan and a GFA
scan run on the same cohort — by matching loci on coordinate overlap and classifying each
pair according to how well the two evidence sources agree.

**Typical workflow:**

```bash
# Step 1: run VCF scan
privy scan --vcf cohort.vcf.gz \
    --targets S1 S2 --off-targets S3 S4 \
    --outdir results/vcf/

# Step 2: run GFA scan
privy scan --gfa pangenome.gfa \
    --targets S1 S2 --off-targets S3 S4 \
    --outdir results/gfa/

# Step 3: compare the two result sets
privy compare \
    --hits-a results/vcf/hits.tsv \
    --hits-b results/gfa/hits.tsv \
    --outdir results/compare/
```

**Match classes:**

| Class | Meaning |
|-------|---------|
| `supported` | Both sources agree; locus is target-private in both |
| `partially_supported` | Overlap present but strictness classes differ (one strict, one relaxed) |
| `contradicted` | One source reports a contradicted locus |
| `source_specific` | Locus found in only one source — not necessarily a failure |
| `uninformative` | Match found but evidence is too weak to classify |
| `missing_data` | One or both sources lack data at this position |

**Key options:**

| Option | Default | Description |
|--------|---------|-------------|
| `--hits-a PATH` | required | hits.tsv from the first scan run (source A) |
| `--hits-b PATH` | required | hits.tsv from the second scan run (source B) |
| `--source-a TEXT` | auto | Display label for source A (inferred from locus_id prefix) |
| `--source-b TEXT` | auto | Display label for source B |
| `--min-reciprocal-overlap FLOAT` | `0.5` | Minimum reciprocal overlap for a coordinate match |
| `--breakpoint-tolerance-bp INT` | `200` | Gap tolerance for near-miss breakpoint matching |
| `--require-state-compatibility` | off | Also require strictness-class compatibility |

**Outputs:**

| File | Contents |
|------|----------|
| `compare.tsv` | One row per locus pair: coordinates, match class, overlap, scores |
| `compare_summary.tsv` | Match-class breakdown with counts, percentages, mean scores |
| `compare.json` | Run metadata and resolved configuration |

> **Source-specific loci are expected.** GFA operates at the segment/bubble level while VCF
> operates at single-nucleotide variants — not every VCF hit will have a GFA counterpart
> and vice versa. `source_specific` rows are informative, not errors.

### privy plot

`privy plot` is planned for v0.6. It will provide visualization of hits, regions, and
cross-source comparison results. See [Current Status](#current-status) for the roadmap.

---

## Current Status

Panex Privus is under active development. Version 0.5.0-dev.

### What works now

- **`privy scan` with VCF input** — fully operational
  - Streaming, indexed VCF reading via pysam
  - All six strictness classes
  - Filter/QUAL/multiallelic controls
  - Region merging (`--merge-distance`)
  - Hit scoring and ranking
  - All six output files
- **`privy scan` with GFA input** — fully operational
  - GFA1 and GFA1.1 (W-line walks) parsing
  - Private-segment discovery using the same StrictnessClass framework
  - Coordinate-based missingness detection (samples absent from a locus vs. traversing
    an alternative bubble arm)
  - Same six output files — directly comparable with VCF output via `privy compare`
- **BAM support layer** — fully operational
  - Pass `--bam sample.bam` (repeat for multiple) or `--bam-manifest manifest.tsv`
  - Depth and allele-fraction queries at each discovered VCF hit locus
  - Target samples: SUPPORT / AMBIGUOUS evidence; off-target: ABSENCE / CONTRADICTION
  - Low-depth and indel loci reported as UNINFORMATIVE (never silently ignored)
  - BAM evidence rows appended to `evidence.tsv`; `depth` and `allele_fraction` columns
    populated in `sample_support.tsv`; `support_score` updated in `hits.tsv`
- **`privy report`** — fully operational
  - Reads `hits.tsv` plus optional `regions.tsv`, `evidence.tsv`, `qc.tsv`, `run.json`
  - Writes `summary.tsv`, `ranked_hits.tsv`, `strictness_summary.tsv`,
    `support_summary.tsv` (when evidence.tsv provided), `contradiction_summary.tsv`
  - Renders `report.md` (default) and/or `report.html` (`--format html|both`)
  - Strictness class distribution table
  - QC section from scan metrics
  - Scientific caveats section
- **`privy compare`** — fully operational
  - Reconciles two `hits.tsv` files (VCF scan vs. GFA scan, or any two scans)
  - Coordinate-overlap matching with reciprocal-overlap and breakpoint-tolerance controls
  - Six match classes: supported / partially_supported / contradicted / source_specific /
    uninformative / missing_data
  - Writes `compare.tsv` (per-locus-pair), `compare_summary.tsv`, `compare.json`
  - Source labels auto-inferred from locus_id prefix (PPX → vcf, GPX → gfa)
- 500+ unit and integration tests passing
- YAML configuration with three-tier priority

### Roadmap

| Version | Focus |
|---------|-------|
| v0.1 | VCF scan, strictness classification, scoring, all outputs |
| v0.2 | GFA scan — standalone pangenome graph discovery |
| v0.3 | `privy report` — ranked summaries and QC reports |
| v0.4 | BAM support layer — read-level evidence at discovered loci |
| v0.5 (current) | `privy compare` — cross-evidence reconciliation between VCF and GFA result sets |
| v0.6 | `privy plot` — visualization of hits, regions, and comparison results |
| v1.0 | Polished docs, example datasets, manuscript-ready outputs, GitHub release |

---

## Troubleshooting

### "VCF index not found"

Your VCF must have a `.tbi` or `.csi` index file in the same directory. Create one with:

```bash
tabix -p vcf your_variants.vcf.gz
```

### "No target samples from the cohort definition were found in the VCF header"

The sample names you passed to `--targets` do not match any names in the VCF. List the VCF's sample names to check spelling:

```bash
bcftools query -l your_variants.vcf.gz
```

Sample names are case-sensitive.

### hits.tsv is empty (0 rows)

Check `qc.tsv` to see how many records were evaluated vs. skipped. Common causes:

- **All records skipped by filter**: try `--no-pass-only` to include non-PASS records
- **All records fail target support threshold**: if `--min-target-support 1.0`, all called targets must carry the allele. Try lowering to `0.8`
- **All alleles contradicted**: off-target samples carry the same alleles. Check whether your cohort definition is correct

### GFA scan produces no hits

Check that your GFA segments have coordinate tags. Open the GFA and look at a few S-lines:

```bash
grep "^S" your_graph.gfa | head -3
```

Each segment needs `SN:Z:chr1  SO:i:1000  LN:i:500`-style tags. Segments without
these tags cannot be placed on the coordinate grid and are skipped.

If your graph uses W-lines (GFA1.1), presence/absence at a locus is detected from the
walk coordinates. Make sure your sample names in W-lines match the names you pass to
`--targets` and `--off-targets` exactly (they are case-sensitive).

### "No target samples from the cohort definition were found in the GFA"

Check how sample names are encoded in your GFA. For W-lines:

```bash
grep "^W" your_graph.gfa | awk '{print $2}' | sort -u
```

For P-lines (the sample name is the part before the first `#`):

```bash
grep "^P" your_graph.gfa | awk '{print $2}' | cut -d'#' -f1 | sort -u
```

Pass those exact names to `--targets` and `--off-targets`.

### pysam installation errors

Install `pysam` via bioconda before running `pip install .`:

```bash
conda install -c bioconda pysam
pip install .
```

---

## Contributing

Contributions are welcome. Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

Areas where help is especially welcome:

- Additional test fixtures and regression cases (especially real-world GFA and BAM files)
- Documentation improvements and worked examples
- `privy plot` (v0.6) — visualization of hits, regions, and comparison results
- Edge-case handling for unusually structured GFA path-name conventions

Please open an issue before making large architectural changes.

---

## Citation

If you use Panex Privus in published work, please cite the software and associated manuscript once available.

```bibtex
@software{panex_privus,
  title  = {Panex Privus: a comparative genomics toolkit for discovering target-private genomic signal},
  author = {Panex Privus Contributors},
  year   = {2026},
  url    = {https://github.com/USDA-ARS-GBRU/Panex_Privus}
}
```

A `CITATION.cff` file is included in this repository for use with GitHub's "Cite this repository" feature.

---

## License

MIT License — see [LICENSE](LICENSE) for details.

---

## Name

**Panex Privus** is the official project name. **`privy`** is the command-line tool.

*Panex* evokes breadth — across the whole pangenome. *Privus* means "set apart" or "belonging to one alone." Together: what belongs to your group, found across the whole genome.
