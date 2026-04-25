---
title: Configuration
description: YAML and cohort configuration patterns for Panex Privus.
---

# Configuration

Panex Privus supports YAML configuration for reproducible runs.

Config priority is:

1. package defaults
2. YAML config file
3. explicit CLI flags

## Example

```yaml
project_name: soybean_privy_scan

cohorts:
  targets: [Benning, Harosoy, Clark]
  off_targets: [Jack, Lee, Minsoy]

scan:
  min_target_support: 1.0
  max_off_target_support: 0.0
  merge_distance: 1000
  strictness_report: true

filters:
  min_qual: 30
  pass_only: true
  allow_multiallelic: true

bam:
  enabled: true
  min_depth: 8
  min_alt_count: 2
  allele_fraction_min: 0.2
  min_mapq: 20
  min_baseq: 20

gfa:
  min_segment_length: 1

compare:
  min_reciprocal_overlap: 0.5
  breakpoint_tolerance_bp: 200

scoring:
  discovery_weight: 1.0
  support_weight: 0.7
  penalty_weight: 0.8
```

Run with:

```bash
privy --config configs/privy.yaml scan --vcf variants.vcf.gz --outdir results/
```

## Cohort Files

You can define cohorts separately from the main config.

YAML:

```yaml
targets:
  - T1
  - T2
off_targets:
  - O1
  - O2
ignored_samples:
  - LowQualitySample
```

TSV:

```text
sample_id	cohort_role
T1	target
T2	target
O1	off_target
O2	off_target
LowQualitySample	ignored
```

Use either with:

```bash
privy scan --vcf variants.vcf.gz --cohort-file cohort.tsv --outdir results/
```

CLI sample flags override cohort-file entries when both are supplied.
