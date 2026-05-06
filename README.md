# Panex Privus (`privy`)

> [!NOTE]
> Panex Privus is under active development. Suggestions, corrections, teaching
> examples, and issue reports are welcome through GitHub Issues or pull requests.

**Find genomic variants shared within your target group and absent from everything else.**

Panex Privus is a comparative genomics toolkit for discovering *target-private*
alleles and graph segments: genomic signal present in a focal cohort and absent
from off-target genomes. It is designed for plant pangenome work and other
comparative genomics projects where the core question is:

> What is genetically unique to my group of interest?

## Background

### What is a target-private allele?

Imagine you are studying soybean breeding lines. You have several varieties that
all share a desirable trait, such as high seed protein content, and a broader
reference set that does not have that trait. A natural question is:

> What DNA variants are present in the high-protein varieties and absent from
> the reference varieties?

Those variants are called *target-private alleles*. They may or may not cause
the trait, but they are candidates worth investigating.

Panex Privus automates finding this signal from a VCF file or a pangenome graph
(GFA). It handles missing data explicitly, merges nearby candidate variants into
regions, and scores results so the highest-confidence candidates surface first.

## What It Does

- Scans indexed multisample **VCF** files for target-private alleles.
- Scans **GFA** pangenome graphs for target-private graph segments.
- Summarizes full, target, and off-target pangenomes from GFA segments or VCF alleles.
- Keeps missing data explicit with `strictness_class` labels.
- Adds optional **BAM** read-level support at VCF hit loci.
- Compares VCF and GFA result sets with coordinate-aware matching.
- Generates reports, plots, GFF3 annotations, and BED/GFF3 exports.

## Quick Start

### Install

The recommended install path is a dedicated mamba environment:

```bash
mamba create -n privy -c conda-forge -c bioconda \
  python=3.11 pysam samtools bcftools htslib
conda activate privy
```

Then install Panex Privus from source:

```bash
git clone https://github.com/USDA-ARS-GBRU/Panex_Privus.git
cd Panex_Privus
python -m pip install .
```

See the [installation guide](https://usda-ars-gbru.github.io/Panex_Privus/installation/)
for pixi, conda, pip, and developer install options.

### Update to the Most Recent Version

If you are following active development from a source checkout, update with:

```bash
conda activate privy
cd Panex_Privus
git pull origin main
python -m pip install -U .
privy --version
```

For an editable developer install, use `python -m pip install -U -e ".[dev]"`
after pulling updates.

### Run a Scan

Run a VCF scan:

```bash
privy scan \
  --vcf cohort.vcf.gz \
  --targets Benning Harosoy Clark \
  --off-targets Jack Lee Minsoy \
  --outdir results/
```

Run a GFA scan:

```bash
privy scan \
  --gfa pangenome.gfa.gz \
  --targets Benning Harosoy Clark \
  --off-targets Jack Lee Minsoy \
  --outdir results/
```

Run both discovery backends and compare them:

```bash
privy scan \
  --vcf cohort.vcf.gz \
  --gfa pangenome.gfa.gz \
  --targets Benning Harosoy Clark \
  --off-targets Jack Lee Minsoy \
  --outdir results/
```

Analyze the graph pangenome and target/off-target sub-pangenomes:

```bash
privy pangenome \
  --gfa pangenome.gfa.gz \
  --targets Benning Harosoy Clark \
  --outdir results/pangenome/
```

Run the same pangenome summaries from VCF alleles:

```bash
privy pangenome \
  --vcf cohort.vcf.gz \
  --targets Benning Harosoy Clark \
  --outdir results/pangenome/
```

Compare two scan outputs:

```bash
privy compare \
  --hits-a results/vcf/hits.tsv \
  --hits-b results/gfa/hits.tsv \
  --outdir results/compare/
```

Generate a report:

```bash
privy report \
  --hits results/vcf/hits.tsv \
  --regions results/vcf/regions.tsv \
  --qc results/vcf/qc.tsv \
  --format both \
  --outdir results/report/
```

See the [run guide](https://usda-ars-gbru.github.io/Panex_Privus/run-guide/)
for the full workflow, command options, BAM support, VCF/GFA comparison,
reports, plots, annotation, and export.

## Documentation

The public documentation site is:

- [Panex Privus documentation](https://usda-ars-gbru.github.io/Panex_Privus/)
- [Installation](https://usda-ars-gbru.github.io/Panex_Privus/installation/)
- [Run guide](https://usda-ars-gbru.github.io/Panex_Privus/run-guide/)
- [Figures and tables](https://usda-ars-gbru.github.io/Panex_Privus/figures-and-tables/)
- [Core concepts](https://usda-ars-gbru.github.io/Panex_Privus/concepts/)
- [Output files](https://usda-ars-gbru.github.io/Panex_Privus/outputs/)
- [Configuration](https://usda-ars-gbru.github.io/Panex_Privus/configuration/)
- [Current status and roadmap](https://usda-ars-gbru.github.io/Panex_Privus/status/)
- [Troubleshooting](https://usda-ars-gbru.github.io/Panex_Privus/troubleshooting/)
- [Architecture](https://usda-ars-gbru.github.io/Panex_Privus/architecture/)
- [Contributing](CONTRIBUTING.md)

## Current Status

Panex Privus is currently `0.8.0-dev`.

Operational commands:

- `privy scan`
- `privy compare`
- `privy report`
- `privy plot`
- `privy annotate`
- `privy export`

The current test suite has 632 passing unit and integration tests.

## Contact

Roth E Conrad - roth.conrad@uga.edu

## Citation

If you use Panex Privus in a publication, please cite the repository and see
[CITATION.cff](CITATION.cff) for citation metadata.

## License

Panex Privus is released under the MIT License. See [LICENSE](LICENSE).

## Acknowledgments

Panex Privus is developed for comparative genomics and plant pangenome research.

### Funding Support
This is a project supported by the U.S. Department of Agriculture - Agricultural Research Service (USDA-ARS) - Genomics and Bioinformatics Research Unit (GBRU) through CRIS Project No. 6066-21310-006-000-D.
Additional project support was through <and any additional agreements or grants>.
