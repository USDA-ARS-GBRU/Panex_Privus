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

## What It Does

- Scans indexed multisample **VCF** files for target-private alleles.
- Scans **GFA** pangenome graphs for target-private graph segments.
- Keeps missing data explicit with `strictness_class` labels.
- Adds optional **BAM** read-level support at VCF hit loci.
- Compares VCF and GFA result sets with coordinate-aware matching.
- Generates reports, plots, GFF3 annotations, and BED/GFF3 exports.

## Quick Start

Install from source:

```bash
git clone https://github.com/USDA-ARS-GBRU/Panex_Privus.git
cd Panex_Privus
pip install .
```

Run a VCF scan:

```bash
privy scan \
  --vcf cohort.vcf.gz \
  --targets Benning Harosoy Clark \
  --off-targets Jack Lee Minsoy \
  --outdir results/vcf/
```

Run a GFA scan:

```bash
privy scan \
  --gfa pangenome.gfa \
  --targets Benning Harosoy Clark \
  --off-targets Jack Lee Minsoy \
  --outdir results/gfa/
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

## Documentation

The detailed user guide lives in `docs/` and is set up for GitHub Pages:

- **Docs site**: <https://usda-ars-gbru.github.io/Panex_Privus/>

- [Documentation home](docs/index.md)
- [Installation](docs/installation.md)
- [Quickstart](docs/quickstart.md)
- [Core concepts](docs/concepts.md)
- [Command reference](docs/commands.md)
- [Output files](docs/outputs.md)
- [Configuration](docs/configuration.md)
- [Current status and roadmap](docs/status.md)
- [Troubleshooting](docs/troubleshooting.md)
- [Architecture](docs/architecture.md)
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

The current test suite has 603 passing unit and integration tests.

## Citation

If you use Panex Privus in a publication, please cite the repository and see
[CITATION.cff](CITATION.cff) for citation metadata.

## License

Panex Privus is released under the MIT License. See [LICENSE](LICENSE).

## Acknowledgments

Panex Privus is developed for comparative genomics and plant pangenome research,
with support from the USDA-ARS Genomics and Bioinformatics Research Unit.
