---
title: Installation
description: Install Panex Privus for use or development.
---

# Installation

## Requirements

| Requirement | Version |
|-------------|---------|
| Python | 3.10 or higher |
| pysam | 0.22 or higher |
| Operating system | Linux or macOS (Windows via WSL) |

Check your Python version:

```bash
python --version
```

If your version is below Python 3.10, install a newer Python environment first.
Miniforge or another conda distribution is a good choice for bioinformatics
workflows.

## Standard Install

```bash
git clone https://github.com/USDA-ARS-GBRU/Panex_Privus.git
cd Panex_Privus
pip install .
```

## Conda-Friendly Install

This route is useful when compiled Python packages are difficult to install with
plain `pip`.

```bash
conda create -n privy python=3.11
conda activate privy
conda install -c bioconda pysam
pip install .
```

## Developer Install

```bash
git clone https://github.com/USDA-ARS-GBRU/Panex_Privus.git
cd Panex_Privus
pip install -e ".[dev]"
```

Verify the install:

```bash
privy --version
privy --help
```

## External Tools

VCF scans require bgzip-compressed, indexed VCF files. The following command-line
tools are often needed before running Privy:

| Tool | What it does |
|------|--------------|
| `bgzip` | Block-compresses VCF files |
| `tabix` | Creates `.tbi` indexes for VCF files |
| `bcftools` | Lists samples, filters records, and manipulates VCFs |
| `samtools` | Sorts and indexes BAM files |

Install with conda:

```bash
conda install -c bioconda samtools bcftools htslib
```

On macOS with Homebrew:

```bash
brew install bcftools htslib samtools
```

On Ubuntu/Debian:

```bash
sudo apt-get install bcftools tabix samtools
```

## Upstream Pangenome Inputs

Panex Privus works well with outputs from
[minigraph-cactus](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md).
A typical run can produce both a `.gfa` and a `.vcf.gz`:

```bash
cactus-pangenome ./js seqFile.txt \
    --outDir pangenome/ --outName my_pangenome \
    --reference RefSample \
    --vcf --giraffe --gfa
```
