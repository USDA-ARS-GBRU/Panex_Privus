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
- Builds **VCF landscape** summaries across sliding windows, including
  missingness, private-allele burden, local sample similarity, and local
  background blocks.
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

All cohort-aware commands accept the same cohort inputs: grouped sample flags
such as `--targets Benning Harosoy Clark`, role-specific files such as
`--targets-file targets.txt`, or a single YAML/TSV cohort file with
`--cohort-file cohort.tsv`.

For large GFA graphs, build the reusable Privy GFA index first. This can take
some time, but later scans auto-detect the sidecar index and skip the expensive
GFA walk-parsing step:

```bash
privy index gfa --gfa pangenome.gfa.gz
```

If you built the index with an older development version, or after pulling a
Privy update that changes GFA indexing, rebuild it once with
`privy index gfa --gfa pangenome.gfa.gz --force`.

Then run a GFA scan:

```bash
privy scan \
  --gfa pangenome.gfa.gz \
  --targets Benning Harosoy Clark \
  --off-targets Jack Lee Minsoy \
  --outdir results/
```

Run both discovery backends and compare them. If
`pangenome.gfa.gz.privy.gfaidx` exists beside the GFA, Privy uses it
automatically:

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

Explore genome-wide VCF landscapes with sliding windows:

```bash
privy landscape \
  --vcf cohort.vcf.gz \
  --targets Benning Harosoy Clark \
  --off-targets Jack Lee Minsoy \
  --window-records 200 \
  --step-records 50 \
  --outdir results/landscape/
```

`privy landscape` complements discovery. It writes per-sample window metrics,
per-window target/off-target summaries, local background blocks, candidate
donor-like or introgressed blocks, per-window pairwise similarity, and run
metadata. The default `--similarity-output full` keeps chromosome-level
similarity plots ready for `privy plot --plot-set landscape`; use
`--similarity-output summary` or `none` for leaner large runs. Use
`--local-pca` for optional local-structure coordinates. Use `--window-bp` and
`--step-bp` when physical base-pair windows are easier to interpret than
fixed-record windows. For filtered SNP-density landscapes, add filters such as
`--variant-type snp`, `--biallelic-only`, `--max-site-missing-rate`, and
`--min-alt-carriers`; `windows.tsv` reports density and the landscape plot set
writes a `variant_density_profile` figure.

Render figures from existing pangenome or landscape results:

```bash
privy plot \
  --plot-set landscape \
  --input-dir results/landscape/ \
  --output-format pdf

privy plot \
  --plot-set pangenome \
  --input-dir results/pangenome/ \
  --output-format pdf
```

Compare two scan outputs:

```bash
privy compare \
  --hits-a results/vcf/hits.tsv \
  --hits-b results/gfa/hits.tsv \
  --outdir results/compare/
```

For VCF/GFA comparisons, `privy compare` normalizes minigraph-cactus contig
names like `Sample#0#Gm01` to `Gm01` and uses contained-overlap matching by
default. Use `--overlap-mode reciprocal` for stricter interval matching.

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
- `privy index`
- `privy landscape`
- `privy pangenome`

The current test suite has 684 passing unit and integration tests.

## Contact

Roth E Conrad - roth.conrad@uga.edu

## Citation

If you use Panex Privus in a publication, please cite the repository and see
[CITATION.cff](CITATION.cff) for citation metadata.

## License

Panex Privus is released under the MIT License. See [LICENSE](LICENSE).

## Acknowledgments

Panex Privus is developed for comparative genomics and plant pangenome research.

The `privy landscape` module is inspired by established VCF and population
genomics tools and methods. Panex Privus focuses on target/off-target-aware
interpretation rather than replacing these projects.

Landscape inspiration and citation sources include:

- [VCFtools](https://vcftools.github.io/man_latest.html) and the VCF format paper
- [scikit-allel](https://scikit-allel.readthedocs.io/)
- [cyvcf2](https://github.com/brentp/cyvcf2)
- [pixy](https://pixy.readthedocs.io/)
- [VCF-kit](https://vcf-kit.readthedocs.io/en/latest/tajima/)
- [PopGenome](https://www.rdocumentation.org/packages/PopGenome/versions/2.7.5/topics/sliding.window.transform-methods)
- [SnpSift Private](https://pcingola.github.io/SnpEff/snpsift/private/)
- [Haploview](https://broadinstitute.org/haploview/blocks-and-haplotypes)
- [LDBlockShow](https://academic.oup.com/bib/article/22/4/bbaa227/5939575)
- [PLINK](https://www.cog-genomics.org/plink/1.9/ibd)
- [BCFtools/RoH](https://samtools.github.io/bcftools/howtos/roh-calling.html)
- [R/qtl2](https://pmc.ncbi.nlm.nih.gov/articles/PMC6366910/)
- [RFMix](https://www.sciencedirect.com/science/article/pii/S0002929713002899)
- [Loter](https://academic.oup.com/mbe/article/35/9/2318/5040668)
- [local PCA/lostruct](https://academic.oup.com/genetics/article/211/1/289/5931130)
- [WinPCA](https://academic.oup.com/bioinformatics/article/41/10/btaf529/8261369)

### Funding Support
This is a project supported by the U.S. Department of Agriculture - Agricultural Research Service (USDA-ARS) - Genomics and Bioinformatics Research Unit (GBRU) through CRIS Project No. 6066-21310-006-000-D.
Additional project support was through <and any additional agreements or grants>.
