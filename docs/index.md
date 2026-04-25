---
title: Documentation
description: User guide for Panex Privus.
---

# Panex Privus Documentation

Panex Privus (`privy`) is a comparative genomics toolkit for discovering
target-private genomic signal: alleles and graph segments shared within a
target cohort and absent from off-target genomes.

This documentation is organized as a user guide rather than one very long
README.

## Start Here

- [Installation]({{ '/installation/' | relative_url }})
- [Quickstart]({{ '/quickstart/' | relative_url }})
- [Core concepts]({{ '/concepts/' | relative_url }})

## Running Analyses

- [Command reference]({{ '/commands/' | relative_url }})
- [Configuration]({{ '/configuration/' | relative_url }})
- [Output files]({{ '/outputs/' | relative_url }})
- [Troubleshooting]({{ '/troubleshooting/' | relative_url }})

## Project Information

- [Current status and roadmap]({{ '/status/' | relative_url }})
- [Architecture]({{ '/architecture/' | relative_url }})
- [Development log]({{ '/dev-log/' | relative_url }})
- [Contributing](https://github.com/USDA-ARS-GBRU/Panex_Privus/blob/main/CONTRIBUTING.md)

## Main Workflow

1. Define target and off-target samples.
2. Run `privy scan` on VCF or GFA input.
3. Optionally add BAM support to VCF scans.
4. Compare VCF and GFA scan outputs with `privy compare`.
5. Generate reports and plots.
6. Annotate or export candidate loci for downstream interpretation.
