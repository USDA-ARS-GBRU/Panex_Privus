---
title: Installation
description: Install Panex Privus for use or development.
---

# Installation

Panex Privus is a Python command-line tool, but it lives in a bioinformatics
ecosystem where compiled packages and external command-line tools matter. The
smoothest path is to create a dedicated environment, install the compiled
bioinformatics dependencies there, then install Panex Privus from this source
tree.

## Requirements

| Requirement | Version |
|-------------|---------|
| Python | 3.10 or higher |
| pysam | 0.22 or higher |
| Operating system | Linux or macOS; Windows users should use WSL |

Check your Python version:

```bash
python --version
```

If your Python is older than 3.10, create a new environment with mamba, pixi,
conda, or another environment manager before installing Panex Privus.

## Option 1: Mamba Install

This is the recommended route for most bioinformatics users. Mamba resolves
conda environments quickly, and Bioconda provides reliable builds of `pysam`,
`samtools`, `bcftools`, and `htslib`.

```bash
mamba create -n privy -c conda-forge -c bioconda \
  python=3.11 pysam samtools bcftools htslib
conda activate privy

git clone https://github.com/USDA-ARS-GBRU/Panex_Privus.git
cd Panex_Privus
python -m pip install .
```

Use this option if you want a clean named environment that works well with other
genomics tools.

## Option 2: Pixi Install

Pixi is a modern conda-compatible project environment manager. This route is
useful if you prefer reproducible, per-project environments with a lockfile.

```bash
git clone https://github.com/USDA-ARS-GBRU/Panex_Privus.git
cd Panex_Privus

pixi init --format pixi --channel conda-forge --channel bioconda
pixi add python=3.11 pip pysam samtools bcftools htslib
pixi run python -m pip install .
```

Then run Privy through pixi:

```bash
pixi run privy --help
```

This creates `pixi.toml` and, after solving, `pixi.lock` in the repository
checkout. Keep them if you want a local locked environment; ignore or remove
them if you only needed a temporary install.

## Option 3: Standard Pip Install

Use plain `pip` if you already have a working Python environment and compiled
dependencies install cleanly on your system.

```bash
git clone https://github.com/USDA-ARS-GBRU/Panex_Privus.git
cd Panex_Privus
python -m pip install .
```

If `pysam` fails to build, switch to the mamba or conda-friendly path so `pysam`
comes from Bioconda instead of being compiled locally by pip.

## Option 4: Conda-Friendly Install

This is the same idea as the mamba route, using conda directly.

```bash
conda create -n privy -c conda-forge -c bioconda \
  python=3.11 pysam samtools bcftools htslib
conda activate privy

git clone https://github.com/USDA-ARS-GBRU/Panex_Privus.git
cd Panex_Privus
python -m pip install .
```

## Option 5: Developer Install

Use an editable install when you plan to modify the code and run tests.

```bash
git clone https://github.com/USDA-ARS-GBRU/Panex_Privus.git
cd Panex_Privus
python -m pip install -e ".[dev]"
```

If you use mamba or conda for development, create and activate the environment
first, then run the editable install command inside it.

## Update to the Most Recent Version

Users following active development can update an existing source checkout
without recreating the environment:

```bash
conda activate privy
cd Panex_Privus
git pull origin main
python -m pip install -U .
privy --version
```

For an editable developer install, run:

```bash
conda activate privy
cd Panex_Privus
git pull origin main
python -m pip install -U -e ".[dev]"
pytest
```

On shared clusters, first make sure `python`, `pip`, and `privy` resolve inside
the same environment:

```bash
which python
which pip
which privy
python -c "import sys; print(sys.executable)"
```

If `pip` says it is defaulting to a user installation, or if paths point to a
cluster module or `~/.local` instead of your conda environment, deactivate or
purge conflicting modules and reactivate the `privy` environment before
installing.

## Verify the Install

After any install path, check that the command is available:

```bash
privy --version
privy --help
```

If you installed with pixi, use:

```bash
pixi run privy --version
pixi run privy --help
```

## External Tools

Panex Privus can read GFA files directly, but VCF and BAM workflows usually
depend on a few standard command-line tools before Privy starts.

| Tool | What it does |
|------|--------------|
| `bgzip` | Block-compresses VCF files into `.vcf.gz` |
| `tabix` | Creates `.tbi` or `.csi` indexes for VCF files |
| `bcftools` | Lists samples, filters records, and manipulates VCFs |
| `samtools` | Sorts and indexes BAM files |

If you used the mamba, pixi, or conda-friendly install above, these tools should
already be in the environment. Otherwise install them separately.

With mamba or conda:

```bash
mamba install -c conda-forge -c bioconda samtools bcftools htslib
```

On macOS with Homebrew:

```bash
brew install bcftools htslib samtools
```

On Ubuntu/Debian:

```bash
sudo apt-get install bcftools tabix samtools
```

## Prepare Inputs

VCF scans require a bgzip-compressed, indexed VCF:

```bash
bgzip -c variants.vcf > variants.vcf.gz
tabix -p vcf variants.vcf.gz
```

BAM support requires coordinate-sorted, indexed BAM files:

```bash
samtools sort sample.bam -o sample.sorted.bam
samtools index sample.sorted.bam
```

GFA scans can read plain-text `.gfa` or gzip-compressed `.gfa.gz` files and do
not require bgzip, tabix, or a separate index. Graph segments need coordinate
tags such as `SN:Z:chr1`, `SO:i:1000`, and `LN:i:500` so Privy can place graph
hits back on genomic coordinates.

## Upstream Pangenome Inputs

Panex Privus works well with outputs from
[minigraph-cactus](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md).
A typical run can produce both a `.gfa.gz` and a `.vcf.gz`:

```bash
cactus-pangenome ./js seqFile.txt \
    --outDir pangenome/ --outName my_pangenome \
    --reference RefSample \
    --vcf --giraffe --gfa
```

Those outputs can be scanned separately or together:

```bash
privy scan \
  --vcf pangenome/my_pangenome.vcf.gz \
  --gfa pangenome/my_pangenome.gfa.gz \
  --targets T1 T2 T3 \
  --off-targets O1 O2 O3 \
  --outdir results/
```
