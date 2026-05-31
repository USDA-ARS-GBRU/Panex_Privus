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
- Builds self-contained interactive HTML dashboards for focus regions, scan
  outputs, landscape outputs, and pangenome outputs.

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

### Cohort Inputs

All cohort-aware commands accept the same cohort inputs: grouped sample flags
such as `--targets Benning Harosoy Clark`, role-specific files such as
`--targets-file targets.txt`, or a single YAML/TSV cohort file with
`--cohort-file cohort.tsv`.

See the [run guide](https://usda-ars-gbru.github.io/Panex_Privus/run-guide/)
for the full workflow, command options, BAM support, VCF/GFA comparison,
reports, plots, annotation, and export.

### `privy scan`

Use `privy scan` to discover target-private alleles from VCF calls, private
graph segments from GFA pangenomes, or both sources in one run.

Run a VCF scan:

```bash
privy scan \
  --vcf cohort.vcf.gz \
  --targets Benning Harosoy Clark \
  --off-targets Jack Lee Minsoy \
  --outdir results/
```

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

More: [Run Guide](https://usda-ars-gbru.github.io/Panex_Privus/run-guide/#vcf-scan) | [Figures](https://usda-ars-gbru.github.io/Panex_Privus/figures-and-tables/#privy-scan) | [Outputs](https://usda-ars-gbru.github.io/Panex_Privus/outputs/#key-hitstsv-columns) | [Architecture](https://usda-ars-gbru.github.io/Panex_Privus/architecture/#discovery-algorithms)

### `privy index gfa`

For large GFA graphs, build the reusable Privy GFA index first. This can take
some time, but later scans auto-detect the sidecar index and skip the expensive
GFA walk-parsing step.

```bash
privy index gfa --gfa pangenome.gfa.gz
```

If you built the index with an older development version, or after pulling a
Privy update that changes GFA indexing, rebuild it once with
`privy index gfa --gfa pangenome.gfa.gz --force`.

More: [Run Guide](https://usda-ars-gbru.github.io/Panex_Privus/run-guide/#gfa-scan) | [Figures](https://usda-ars-gbru.github.io/Panex_Privus/figures-and-tables/#privy-scan) | [Outputs](https://usda-ars-gbru.github.io/Panex_Privus/outputs/#gfa-graph-segment-outputs) | [Architecture](https://usda-ars-gbru.github.io/Panex_Privus/architecture/#implementation-architecture)

### `privy pangenome`

Use `privy pangenome` to summarize full, target, and off-target pangenomes from
GFA graph segments or VCF alleles.

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

More: [Run Guide](https://usda-ars-gbru.github.io/Panex_Privus/run-guide/#analyze-the-pangenome) | [Figures](https://usda-ars-gbru.github.io/Panex_Privus/figures-and-tables/#privy-pangenome) | [Outputs](https://usda-ars-gbru.github.io/Panex_Privus/outputs/#pangenome-outputs) | [Architecture](https://usda-ars-gbru.github.io/Panex_Privus/architecture/#pangenome-feature-architecture)

### `privy landscape`

Use `privy landscape` to explore genome-wide VCF behavior with sliding windows,
including missingness, private ALT burden, local similarity, local background
blocks, and candidate donor-like intervals.

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

More: [Run Guide](https://usda-ars-gbru.github.io/Panex_Privus/run-guide/#explore-vcf-landscapes) | [Figures](https://usda-ars-gbru.github.io/Panex_Privus/figures-and-tables/#privy-landscape) | [Outputs](https://usda-ars-gbru.github.io/Panex_Privus/outputs/#landscape-outputs) | [Architecture](https://usda-ars-gbru.github.io/Panex_Privus/architecture/#landscape-algorithms)

### `privy plot`

Use `privy plot` to render publication-review figures from existing scan,
pangenome, or landscape result directories.

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

More: [Run Guide](https://usda-ars-gbru.github.io/Panex_Privus/run-guide/#plot-diagnostics) | [Figures](https://usda-ars-gbru.github.io/Panex_Privus/figures-and-tables/#privy-plot) | [Outputs](https://usda-ars-gbru.github.io/Panex_Privus/outputs/#pangenome-outputs) | [Architecture](https://usda-ars-gbru.github.io/Panex_Privus/architecture/#implementation-architecture)

### `privy interactive`

Use `privy interactive` to build self-contained HTML dashboards for candidate
focus regions and for existing scan, landscape, or pangenome result directories.
For collaborator sharing, a trusted HTTPS page is best. If you share through
OneDrive, SharePoint, Outlook, or Teams, send a zip and ask recipients to
download and open the `.html` file directly in Edge, Chrome, or Firefox rather
than using the cloud preview pane. Privy writes `SHARING_NOTES.txt` beside each
dashboard with this guidance for recipients and IT support.

Build a shareable interactive dashboard for a focus region. In the current
development implementation, `--focus` names a genomic region to render. Start
with a region around 4 Mbp or smaller unless you already know the variant and
annotation density is low enough for a larger static HTML file.

```bash
privy interactive \
  --focus Gm15:1-4000000 \
  --vcf cohort.vcf.gz \
  --gff3 Wm82.gene_exons.gff3.gz \
  --samples Harosoy Harosoy-sharp Kingawa \
  --track-gff RepeatMasker=Wm82.repeats.gff3.gz \
  --outdir results/interactive/
```

Repeat `--focus` to build one HTML file per focus region. When multiple focus
regions are supplied, Privy also writes an `index.html` linking to each region
dashboard. Use `--sites-tsv` instead of `--vcf` when you want to rebuild a
dashboard from a previously extracted focus-region genotype table.
Add `--functional-tsv` to join gene-level annotation and `--keyword-group`
to create phenotype-oriented feature lists; the run guide defines the TSV
schema and keyword matching behavior.

Build an interactive dashboard from existing `privy scan` outputs:

```bash
privy interactive \
  --scan results/scan/ \
  --max-hits 5000 \
  --max-regions 1000 \
  --outdir results/interactive/
```

`--scan` accepts either a direct scan source directory containing `hits.tsv`,
or a combined scan directory with `vcf/`, `gfa/`, and optional `compare/`
children.

Build an interactive dashboard from existing `privy landscape` outputs:

```bash
privy interactive \
  --landscape results/landscape/ \
  --max-windows 20000 \
  --max-sample-windows 80000 \
  --outdir results/interactive/
```

The landscape dashboard summarizes windows, sample-by-window metrics, local
background similarity, and candidate introgression blocks from an existing
landscape output directory.

Build an interactive dashboard from existing `privy pangenome` outputs:

```bash
privy interactive \
  --pangenome results/pangenome/ \
  --max-features 10000 \
  --max-private-features 5000 \
  --outdir results/interactive/
```

`--pangenome` accepts either a direct pangenome directory containing
`feature_summary.tsv`, or a combined pangenome directory with `vcf/` and `gfa/`
children.

Example interactive dashboards generated from the checked-in example outputs:

- [Scan dashboard](https://usda-ars-gbru.github.io/Panex_Privus/assets/examples/interactive/scan/scan_dashboard.html)
- [Landscape dashboard](https://usda-ars-gbru.github.io/Panex_Privus/assets/examples/interactive/landscape/landscape_dashboard.html)
- [Pangenome dashboard](https://usda-ars-gbru.github.io/Panex_Privus/assets/examples/interactive/pangenome-gfa/pangenome_dashboard.html)
- [Focus-region browser](https://usda-ars-gbru.github.io/Panex_Privus/assets/examples/interactive/focus/focus_chr1_1_1000.html)

More: [Run Guide](https://usda-ars-gbru.github.io/Panex_Privus/run-guide/#interactive-focus-regions) | [Figures](https://usda-ars-gbru.github.io/Panex_Privus/figures-and-tables/#privy-interactive) | [Outputs](https://usda-ars-gbru.github.io/Panex_Privus/outputs/#interactive-outputs) | [Architecture](https://usda-ars-gbru.github.io/Panex_Privus/architecture/#implementation-architecture)

### `privy compare`

Use `privy compare` to reconcile two scan outputs, such as VCF and GFA
candidate sets, with coordinate-aware interval matching.

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

More: [Run Guide](https://usda-ars-gbru.github.io/Panex_Privus/run-guide/#compare-existing-scan-outputs) | [Figures](https://usda-ars-gbru.github.io/Panex_Privus/figures-and-tables/#privy-compare) | [Outputs](https://usda-ars-gbru.github.io/Panex_Privus/outputs/#compare-outputs) | [Architecture](https://usda-ars-gbru.github.io/Panex_Privus/architecture/#cross-source-comparison)

### `privy report`

Use `privy report` to package scan outputs into collaborator-friendly Markdown
or HTML summaries.

Generate a report:

```bash
privy report \
  --hits results/vcf/hits.tsv \
  --regions results/vcf/regions.tsv \
  --qc results/vcf/qc.tsv \
  --format both \
  --outdir results/report/
```

More: [Run Guide](https://usda-ars-gbru.github.io/Panex_Privus/run-guide/#generate-a-report) | [Figures](https://usda-ars-gbru.github.io/Panex_Privus/figures-and-tables/#privy-report) | [Outputs](https://usda-ars-gbru.github.io/Panex_Privus/outputs/#report-outputs) | [Architecture](https://usda-ars-gbru.github.io/Panex_Privus/architecture/#implementation-architecture)

### `privy annotate`

Use `privy annotate` to intersect candidate private loci with GFF3 gene models
or other feature annotations.

```bash
privy annotate \
  --hits results/vcf/hits.tsv \
  --gff annotation.gff3.gz \
  --outdir results/annotated/
```

More: [Run Guide](https://usda-ars-gbru.github.io/Panex_Privus/run-guide/#annotate-hits) | [Figures](https://usda-ars-gbru.github.io/Panex_Privus/figures-and-tables/#privy-annotate) | [Outputs](https://usda-ars-gbru.github.io/Panex_Privus/outputs/#annotate-outputs) | [Architecture](https://usda-ars-gbru.github.io/Panex_Privus/architecture/#implementation-architecture)

### `privy export`

Use `privy export` to write BED or GFF3 intervals for genome browsers and
downstream interval tools.

```bash
privy export \
  --hits results/vcf/hits.tsv \
  --regions results/vcf/regions.tsv \
  --format gff3 \
  --outdir results/exported/
```

More: [Run Guide](https://usda-ars-gbru.github.io/Panex_Privus/run-guide/#export-intervals) | [Figures](https://usda-ars-gbru.github.io/Panex_Privus/figures-and-tables/#privy-export) | [Outputs](https://usda-ars-gbru.github.io/Panex_Privus/outputs/#export-outputs) | [Architecture](https://usda-ars-gbru.github.io/Panex_Privus/architecture/#implementation-architecture)

## Comparative Pangenome Layer (new — under active development)

> [!NOTE]
> The comparative-pangenome commands below are new and currently
> **synthetic-data-validated**: they pass an extensive unit + integration test
> suite on built-in synthetic graphs, but have **not yet been validated on real
> crop pangenomes or benchmarked against established tools** (that work is planned
> on the UGA Sapelo2 cluster — see the [validation plan](docs/benchmarking.md)).
> Treat results as candidate signal pending that validation.

This layer consumes an existing pangenome graph (GFA from minigraph-cactus or
PGGB) and adds a comparative + breeder-facing layer on top of the target-private
core: derive syntenic regions, project any coordinate to any reference, resolve
multi-allelic microhaplotypes, compute breeder population genetics, and render
interactive visualizations. Full guide:
[Comparative pangenome workflows](docs/comparative-pangenome.md).

### `privy synteny`

Graph-native (or PAF-anchored) synteny with typed rearrangements
(collinear / inversion / translocation / duplication) and **target-private
structural regions**.

```bash
privy synteny --gfa pangenome.gfa --reference Wm82#0#Gm01 \
  --targets Benning,Harosoy --off-targets Jack,Lee --outdir results/synteny/
# also: --paf untangle.paf   (anchors from odgi untangle / minimap2 / wfmash)
```

### `privy project`

Project a region — or a raw set of graph segments — onto **any** reference
genome in the graph (and lift annotation tracks across genomes).

```bash
privy project --gfa pangenome.gfa --region Wm82#0#Gm01:1000000-1500000 --outdir results/project/
privy project --gfa pangenome.gfa --node-set s120,s121,s122 --outdir results/project/
```

### `privy microhap`

Detect local **multi-allelic microhaplotypes** and flag **target-private
alleles**; writes per-locus and allele-matrix tables.

```bash
privy microhap --gfa pangenome.gfa --reference Wm82#0#Gm01 \
  --targets Benning,Harosoy --off-targets Jack,Lee --outdir results/microhap/
```

### `privy popgen`

Breeder population genetics: allelic diversity, target-vs-off-target
differentiation (F_ST / Jost's D, **diagnostic markers**), PCA, F_IS, private
allelic richness, plus **GP-ready dosage matrix + VanRaden GRM**.

```bash
privy popgen --gfa pangenome.gfa --reference Wm82#0#Gm01 \
  --targets Benning,Harosoy --off-targets Jack,Lee --outdir results/popgen/
```

### `privy dashboard`

Build a self-contained, offline, interactive HTML dashboard (linked
riparian ↔ dotplot, target-private highlighting, SVG export, allele panel) from
synteny (and optional microhap) outputs.

```bash
privy dashboard --synteny results/synteny/ --microhap results/microhap/
```

More: [Comparative pangenome guide](docs/comparative-pangenome.md) |
[Team tutorial](docs/team-guide.md) |
[Testing guide](docs/testing-guide.md) |
[Sapelo2 live runs](docs/sapelo2-runs.md) |
[Validation plan](docs/benchmarking.md)

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

Panex Privus is currently `0.9.1`.

Operational commands (released line):

- `privy scan`
- `privy compare`
- `privy report`
- `privy plot`
- `privy annotate`
- `privy export`
- `privy index gfa`
- `privy landscape`
- `privy pangenome`
- `privy interactive`

Comparative-pangenome layer (new; synthetic-validated, on the
`feature/comparative-pangenome-p0` branch, pending real-data validation):

- `privy synteny`
- `privy project`
- `privy microhap`
- `privy popgen`
- `privy dashboard`

The released test suite has 701 passing tests; with the comparative-pangenome
layer the branch suite has **944 passing** unit + integration tests (1 skipped
when the optional `scikit-learn` extra is absent).

Compact version history and roadmap:

| Version | Focus | Status |
|---------|-------|--------|
| `v0.1` | VCF private-allele scan, missingness-aware scoring, candidate regions, reports | Complete |
| `v0.2` | GFA private-segment scan and coordinate-aware VCF/GFA comparison | Complete |
| `v0.3` | Plotting, annotation, BED/GFF3 export, and documentation examples | Complete |
| `v0.4` | Reproducible configuration, cohort files, validation, and packaging hardening | Complete |
| `v0.5` | Pangenome feature summaries, composition tables, coverage histograms, and growth curves | Complete |
| `v0.6` | VCF landscape summaries and exploratory local-background blocks | Complete |
| `v0.7` | Publication-oriented figure, report, and output documentation | Complete |
| `v0.8` | GFA indexing, larger example workflows, and reviewer-facing architecture docs | Complete |
| `v0.9` | VCF landscape summaries, local background blocks, and candidate donor-like intervals | Complete |
| `v0.9.1` | Interactive HTML dashboards for scan, landscape, pangenome, and focus-region review | Current (released line) |
| `v1.0` | Manuscript-ready release hardening, archive-ready examples, and expanded benchmark validation | Planned |

Comparative-pangenome layer (branch `feature/comparative-pangenome-p0`; ✅ Built =
implemented + synthetic-data-validated in CI; real-data validation pending):

| Version | Focus | Status |
|---------|-------|--------|
| `v1.1` | PAF/GFA path-coordinate model; dependency tiers; synthetic fixtures | ✅ Built |
| `v1.2` | Coordinate projection engine (`privy project`) | ✅ Built |
| `v1.3` | Graph-native typed synteny (`privy synteny`, GFA + PAF) | ✅ Built |
| `v1.4` | Alignment- and gene-anchored synteny chainer | ✅ Built |
| `v1.4.5` | Microhaplotype / allele-space layer (`privy microhap`); MADC + PHG hVCF I/O | ✅ Built |
| `v1.5` | Static riparian / dotplot / density figures (`privy plot --plot-set synteny`) | ✅ Built |
| `v1.6` | Interactive comparative dashboard (`privy dashboard`) | ✅ Built (GenomeSpy track-browser deferred) |
| `v1.7` | Polyploid dosage + GP matrices; chromosome-structure binning | ✅ Built |
| `v1.7.5` | Breeder population genetics (`privy popgen`) | ✅ Built |
| `v2.0` | vg-deconstruct ingest; Sapelo2 real-data validation + benchmarks + manuscript | ◐ In progress |

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
