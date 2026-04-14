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

Panex Privus automates finding these variants from a standard multi-sample VCF file. It handles:

- **Missingness explicitly** — some samples may have no genotype call at a position due to low coverage or quality filtering. Panex Privus never silently pretends missing data is an absence. Every result is labeled with a *strictness class* that tells you exactly what data was available.
- **Region merging** — nearby candidate variants are merged into genomic windows for downstream analysis.
- **Scoring and ranking** — results are scored by how confidently they fit the target-private pattern, so your highest-confidence candidates come first.

### What do I need to use this?

- A **multi-sample VCF file** containing genotype calls for all your target and off-target samples in a single file.
- The VCF must be **compressed with bgzip** and **indexed with tabix** (standard practice in bioinformatics).
- Sample names in the VCF header that you can pass to `privy scan`.

If your samples are in separate VCF files, you can merge them first using `bcftools merge`.

---

## Requirements

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

This walkthrough takes you from a VCF file to a ranked list of candidate private alleles.

### Step 1: Make sure your VCF is indexed

Panex Privus requires a bgzip-compressed, tabix-indexed VCF. If your file is a plain `.vcf`, run:

```bash
bgzip -c your_variants.vcf > your_variants.vcf.gz
tabix -p vcf your_variants.vcf.gz
```

If you already have a `.vcf.gz` file but no `.tbi` file alongside it:

```bash
tabix -p vcf your_variants.vcf.gz
```

### Step 2: List the sample names in your VCF

```bash
bcftools query -l your_variants.vcf.gz
```

You will see one sample name per line. These are the names you will pass to `--targets` and `--off-targets`.

### Step 3: Run the scan

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

### Step 4: Explore your results

```bash
head -5 results/hits.tsv
```

The first rows are your highest-confidence candidates (rank 1 = best score).

---

## Understanding the Output Files

After running `privy scan`, your output directory contains:

| File | What it contains |
|------|-----------------|
| `hits.tsv` | One row per candidate private allele, sorted by confidence score (rank 1 = highest). Start here. |
| `regions.tsv` | Nearby hits merged into genomic intervals. Useful for downstream analysis. |
| `evidence.tsv` | One evidence record per hit, showing the source of each finding (VCF in Phase 1). |
| `sample_support.tsv` | Per-sample genotype at every hit locus. Tells you exactly which samples carried the allele and which were missing. |
| `qc.tsv` | Run-level quality control: how many records were evaluated, skipped, and passed. Review this to understand how your filters affected the analysis. |
| `run.json` | Complete record of every parameter used in this run. Keep this file to reproduce your results later. |

### Key columns in hits.tsv

| Column | Type | Description |
|--------|------|-------------|
| `locus_id` | text | Unique identifier for this locus (e.g., `PPX00000001`) |
| `contig` | text | Chromosome or scaffold name |
| `start` | integer | Start position (0-based, half-open — same convention as BED files) |
| `end` | integer | End position (0-based, half-open) |
| `allele_key` | text | Human-readable variant description: `contig:pos:REF:ALT` (1-based VCF position) |
| `variant_type` | text | `snp`, `indel`, or `sv` (structural variant) |
| `target_support_n` | integer | Number of target samples carrying this allele |
| `target_total_n` | integer | Total number of target samples |
| `offtarget_support_n` | integer | Number of off-target samples carrying this allele (should be 0) |
| `offtarget_total_n` | integer | Total number of off-target samples |
| `target_missing_n` | integer | Target samples with no genotype call at this position |
| `offtarget_missing_n` | integer | Off-target samples with no genotype call at this position |
| `strictness_class` | text | Confidence classification (see next section) |
| `final_score` | float | Composite confidence score. Higher is better. |

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

The primary discovery engine. Finds target-private alleles in a VCF.

```bash
privy scan --vcf PATH --targets SAMPLE [SAMPLE ...] --off-targets SAMPLE [SAMPLE ...] [OPTIONS]
```

Key options:

| Option | Default | Description |
|--------|---------|-------------|
| `--vcf PATH` | required | Indexed multisample VCF (.vcf.gz + .tbi) |
| `--targets TEXT` | required | Target sample names (repeat for each sample) |
| `--off-targets TEXT` | required | Off-target sample names (repeat for each sample) |
| `--outdir PATH` | `.` | Output directory |
| `--min-target-support FLOAT` | `1.0` | Minimum fraction of called targets that must carry the allele |
| `--max-off-target-support FLOAT` | `0.0` | Maximum fraction of called off-targets allowed to carry the allele |
| `--merge-distance INT` | `0` | Merge loci within this many bp into regions (0 = no merging) |
| `--pass-only / --no-pass-only` | `true` | Require `FILTER=PASS` in the VCF |
| `--min-qual FLOAT` | none | Minimum VCF QUAL score |
| `--allow-multiallelic` | `true` | Evaluate multiallelic records (one alt allele at a time) |
| `--region TEXT` | none | Restrict scan to a region (format: `chr1:1000-2000`) |
| `--contig TEXT` | none | Restrict scan to a single contig |

Run `privy scan --help` for the full option list.

### privy compare, privy report, privy plot

These subcommands are under active development and will be available in upcoming releases. See [Current Status](#current-status) for the roadmap.

---

## Current Status

Panex Privus is under active development. Version 0.1.0-dev.

### What works now

- **`privy scan` with VCF input** — fully operational
  - Streaming, indexed VCF reading via pysam
  - All six strictness classes
  - Filter/QUAL/multiallelic controls
  - Region merging (`--merge-distance`)
  - Hit scoring and ranking
  - All six output files
- 183 unit and integration tests passing
- YAML configuration with three-tier priority

### Roadmap

| Version | Focus |
|---------|-------|
| v0.1 (current) | VCF scan, strictness classification, scoring, all outputs |
| v0.2 | `privy report` — ranked summaries and QC reports |
| v0.3 | BAM support layer — read-level evidence at discovered loci |
| v0.4 | GFA support layer — graph-context annotation |
| v0.5 | XMFA support, `privy compare` — cross-evidence reconciliation |
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

- Additional test fixtures and regression cases
- Documentation improvements and worked examples
- BAM support layer (Phase 3)
- Report and plot commands

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
