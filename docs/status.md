---
title: Current Status and Roadmap
description: Current Panex Privus command status, roadmap, and contribution areas.
---

# Current Status and Roadmap

Panex Privus is under active development. Current version: `0.8.1-dev`.

## What Works Now

- `privy scan` with VCF input
- `privy scan` with GFA input
- BAM support for VCF hits
- `privy compare`
- `privy report`
- `privy plot`
- `privy annotate`
- `privy export` to BED and GFF3
- `privy index gfa`
- `privy pangenome` with GFA or VCF input
- `privy landscape` with VCF input

The current test suite has 699 passing unit and integration tests.

## Roadmap

| Version | Focus |
|---------|-------|
| v0.1 | VCF scan, strictness classification, scoring, and scan outputs |
| v0.2 | GFA scan as standalone graph discovery |
| v0.3 | `privy report` |
| v0.4 | BAM support layer |
| v0.5 | `privy compare` |
| v0.6 | `privy plot` |
| v0.7 | `privy annotate` |
| v0.8 | `privy export` to BED/GFF3 |
| v0.9 | `privy landscape` VCF sliding-window summaries, local background maps, and candidate donor-like blocks |
| v1.0 | Polished docs, example datasets, manuscript-ready outputs, release hardening |

## Helpful Contributions

- Real-world regression fixtures, especially for GFA and BAM edge cases
- Worked examples for new users
- Multi-cohort workflows
- Annotated VCF-style export
- Landscape examples from real VCFs, especially pangenome-derived VCFs and
  controlled-cross/MAGIC populations
- Documentation improvements
