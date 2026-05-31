---
title: Documentation
description: User guide for Panex Privus.
---

# Panex Privus Documentation

Panex Privus (`privy`) is a comparative genomics toolkit for discovering
target-private genomic signal: alleles and graph segments shared within a
target cohort and absent from off-target genomes.

The package also includes pangenome summaries, VCF landscape analysis for
sliding-window context, and self-contained interactive HTML dashboards for
sharing focus-region browsers and run-level summaries with collaborators.

This documentation is organized as a user guide rather than one very long
README.

## Start Here

- [Installation]({{ '/installation/' | relative_url }})
- [Run Guide]({{ '/run-guide/' | relative_url }})
- [Core concepts]({{ '/concepts/' | relative_url }})

## Running Analyses

- [Configuration]({{ '/configuration/' | relative_url }})
- [Output files]({{ '/outputs/' | relative_url }})
- [Figures and tables]({{ '/figures-and-tables/' | relative_url }})
- [Troubleshooting]({{ '/troubleshooting/' | relative_url }})

## Comparative Pangenome (new — under active development)

A new layer that consumes a pangenome graph (GFA) and adds synteny, coordinate
projection, microhaplotypes, breeder population genetics, and interactive
visualization on top of the target-private core. Synthetic-data-validated in CI;
real-data validation pending (see the validation plan).

- [Comparative pangenome workflows]({{ '/comparative-pangenome/' | relative_url }}) — the command guide
- [Team tutorial]({{ '/team-guide/' | relative_url }}) — hands-on, concepts + worked examples
- [Testing guide]({{ '/testing-guide/' | relative_url }}) — how to run and read the tests
- [Sapelo2 live runs]({{ '/sapelo2-runs/' | relative_url }}) — per-crop real-data recipes
- [Validation & benchmarking plan]({{ '/benchmarking/' | relative_url }})

## Project Information

- [Current status and roadmap]({{ '/status/' | relative_url }})
- [Architecture]({{ '/architecture/' | relative_url }})
- [Development log]({{ '/dev-log/' | relative_url }})
- [Archived longform README]({{ '/archive/longform-readme/' | relative_url }})
- [Contributing](https://github.com/USDA-ARS-GBRU/Panex_Privus/blob/main/CONTRIBUTING.md)

## Main Workflow

1. Define target and off-target samples.
2. For large GFA inputs, build the reusable sidecar index with `privy index gfa`.
3. Run `privy scan` on VCF or GFA input.
4. Optionally add BAM support to VCF scans.
5. Compare VCF and GFA scan outputs with `privy compare`.
6. Run `privy pangenome` when you want feature composition and growth context.
7. Run `privy landscape` when you want genome-wide VCF window context.
8. Build `privy interactive` dashboards for browsable review and sharing.
9. Generate reports and plots.
10. Annotate or export candidate loci for downstream interpretation.
