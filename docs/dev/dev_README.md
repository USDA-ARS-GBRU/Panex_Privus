Below is a launch-ready first pass you can hand directly to Claude Code. I’ve written it in the voice and structure of a serious open-source bioinformatics package, not as notes.

⸻

README.md

# Panex Privus

**Panex Privus** is a comparative genomics toolkit for discovering **target-private genomic signal**: alleles and regions shared within a focal cohort and absent from off-target genomes.

Its command-line interface is:

privy

Panex Privus is designed for VCF-first discovery, with optional support and comparison layers from:
	•	VCF — primary discovery backend
	•	BAM — read-level support and contradiction checks
	•	GFA — graph-context annotation
	•	XMFA — alignment-based corroboration

The package is built for a single central question:

What genomic signal is genuinely shared within my target group and not shared outside it?

⸻

Why Panex Privus exists

Many existing genomics tools can call variants, build graphs, align assemblies, or extract regions. Fewer tools are designed around a stricter biological inference problem:
	•	define a target cohort
	•	define an off-target cohort
	•	identify private alleles or private regions
	•	evaluate missingness separately from biological contradiction
	•	compare support across multiple evidence layers
	•	produce outputs that are useful for both computation and publication

Panex Privus is built for that problem.

⸻

Core concepts

Private allele

A candidate allele that is present in the target cohort and absent from the off-target cohort.

Strictness class

A reported classification that separates biological support from technical missingness.

For example, a locus may be:
	•	strictly supported with complete data
	•	consistent with the target-private model but missing one or more target calls
	•	consistent with the target-private model but missing one or more off-target calls
	•	relaxed-threshold support only
	•	contradicted by off-target evidence

Cross-evidence comparison

A locus can be assessed across multiple data types. A private allele identified in VCF may be:
	•	supported by BAM read evidence
	•	located near a GFA branch or junction
	•	corroborated by XMFA alignment structure
	•	contradicted by off-target read support
	•	source-specific or uninformative

⸻

Features

privy scan

Primary discovery engine.
	•	VCF-first target-private allele discovery
	•	target/off-target cohort logic
	•	explicit missingness classification
	•	interval/region merging
	•	optional BAM, GFA, and XMFA support overlays
	•	machine-friendly and publication-friendly outputs

privy compare

Cross-evidence reconciliation engine.
	•	compare loci or regions across sources
	•	classify support, contradiction, or source-specific evidence
	•	quantify overlap and compatibility

privy report

Interpretation engine.
	•	ranked hit summaries
	•	QC summaries
	•	strictness class distributions
	•	support and contradiction summaries
	•	Markdown and optional HTML output

privy plot

Focused evidence visualization.
	•	locus-level evidence panels
	•	region summaries
	•	genotype heatmaps
	•	support and contradiction plots
	•	BAM depth and allele fraction panels
	•	GFA and XMFA context overlays where available

⸻

Supported inputs

Primary
	•	bgzip/tabix indexed multisample VCF
	•	one or more cohort definitions: targets and off-targets

Secondary support layers
	•	BAM files mapped to the same reference
	•	GFA files for graph-context annotation
	•	XMFA files for alignment-based corroboration

Current design priorities

Panex Privus is optimized for:
	•	plant pangenome scale
	•	multi-sample comparative genomics
	•	interval-centric outputs
	•	streaming and chunked processing
	•	auditable filtering and scoring

⸻

Installation

From source

git clone https://github.com/<your-org-or-user>/panex-privus.git
cd panex-privus
pip install .

Development install

git clone https://github.com/<your-org-or-user>/panex-privus.git
cd panex-privus
pip install -e ".[dev]"


⸻

Quick start

Minimal VCF-first scan

privy scan \
  --vcf cohort.vcf.gz \
  --targets Benning Harosoy Clark \
  --off-targets Jack Lee Minsoy \
  --mode private_allele \
  --outdir results/

Scan with BAM and GFA support

privy scan \
  --vcf cohort.vcf.gz \
  --targets Benning Harosoy Clark \
  --off-targets Jack Lee Minsoy \
  --mode private_allele \
  --bam bam_manifest.tsv \
  --gfa graph.gfa.gz \
  --merge-distance 1000 \
  --outdir results/

Compare evidence across layers

privy compare \
  --hits results/hits.tsv \
  --vcf cohort.vcf.gz \
  --bam bam_manifest.tsv \
  --gfa graph.gfa.gz \
  --mode multi_evidence \
  --outdir compare/

Generate a report

privy report \
  --hits results/hits.tsv \
  --regions results/regions.tsv \
  --evidence results/evidence.tsv \
  --compare compare/compare.tsv \
  --qc results/qc.tsv \
  --outdir report/

Plot a locus

privy plot \
  --hits results/hits.tsv \
  --locus-id PPX000123 \
  --vcf cohort.vcf.gz \
  --bam bam_manifest.tsv \
  --gfa graph.gfa.gz \
  --outdir plots/


⸻

YAML configuration

Panex Privus supports YAML configuration files for reproducible runs.

Example privy.yaml

project_name: soybean_privy_scan
mode: private_allele

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
  allele_fraction_min: 0.2
  summarize_softclips: true

gfa:
  enabled: true
  junction_window_bp: 1000
  report_path_membership: true

xmfa:
  enabled: false
  gap_aware: true

compare:
  overlap_mode: reciprocal
  min_reciprocal_overlap: 0.5
  breakpoint_tolerance_bp: 200

scoring:
  discovery_weight: 1.0
  support_weight: 0.7
  penalty_weight: 0.8

CLI flags override config values.

⸻

Output files

hits.tsv

One row per passing locus.

Typical columns:
	•	locus_id
	•	contig
	•	start
	•	end
	•	variant_type
	•	allele_key
	•	target_support_n
	•	target_total_n
	•	offtarget_support_n
	•	offtarget_total_n
	•	target_missing_n
	•	offtarget_missing_n
	•	strictness_class
	•	discovery_score
	•	support_score
	•	penalty_score
	•	final_score

regions.tsv

Merged candidate regions.

Typical columns:
	•	region_id
	•	contig
	•	start
	•	end
	•	n_loci
	•	variant_types
	•	dominant_strictness_class
	•	target_consistency
	•	offtarget_exclusion
	•	final_score

evidence.tsv

Normalized evidence records from VCF, BAM, GFA, and XMFA.

sample_support.tsv

Per-sample support summaries.

qc.tsv

Run-level QC and filtering metrics.

run.json

Full run provenance, parameterization, and scoring weights.

compare.tsv

Cross-source comparison records.

⸻

Strictness classes

Missingness is reported explicitly rather than silently hidden in pass/fail logic.

Expected strictness classes include:
	•	strict_complete
	•	strict_target_missing
	•	strict_offtarget_missing
	•	strict_both_missing
	•	relaxed_threshold
	•	contradicted

This makes the software more auditable and easier to interpret in real datasets.

⸻

Comparison classes

Cross-source comparisons classify evidence into:
	•	supported
	•	partially_supported
	•	contradicted
	•	source_specific
	•	uninformative
	•	missing_data

These classes are designed to distinguish:
	•	real biological agreement
	•	technical incompleteness
	•	evidence-layer disagreement
	•	source-specific discoveries

⸻

Design principles

Panex Privus is built around a few explicit principles:
	1.	VCF-first discovery
VCF is the primary discovery backend in v1.
	2.	Format-aware, logic-centered
File formats do not define truth. Each source contributes evidence for or against target-private status.
	3.	Missingness is not contradiction
Missing data should be surfaced explicitly, not buried.
	4.	Plant pangenome scale matters
The package is designed for indexed streaming, chunked processing, and interval-based work.
	5.	Cross-evidence comparison is a first-class feature
BAM, GFA, and XMFA are not decorative add-ons. They inform support, contradiction, and uncertainty.

⸻

Repository layout

panex-privus/
├── README.md
├── LICENSE
├── pyproject.toml
├── CHANGELOG.md
├── CITATION.cff
├── CONTRIBUTING.md
├── docs/
│   ├── overview.md
│   ├── installation.md
│   ├── cli.md
│   ├── config.md
│   ├── formats.md
│   ├── scoring.md
│   ├── examples.md
│   └── architecture.md
├── src/
│   └── privy/
├── tests/
│   ├── unit/
│   ├── integration/
│   ├── regression/
│   └── data/
└── .github/
    └── workflows/


⸻

Development roadmap

v0.1
	•	package scaffold
	•	YAML config
	•	cohort model
	•	VCF private-allele scan
	•	strictness classification
	•	hits/qc/run outputs

v0.2
	•	region merging
	•	scoring
	•	report command
	•	basic plots

v0.3
	•	BAM support layer
	•	compare v1: VCF vs BAM

v0.4
	•	GFA support layer
	•	compare v2: VCF vs GFA

v0.5
	•	XMFA support layer
	•	multi-evidence compare modes

v1.0
	•	polished docs
	•	example workflows
	•	benchmark datasets
	•	GitHub release
	•	citation metadata
	•	manuscript-ready outputs

⸻

Status

Panex Privus is under active development.

The initial release is focused on:
	•	VCF-first target-private discovery
	•	explicit strictness classification
	•	interval-aware outputs
	•	BAM/GFA/XMFA support architecture
	•	reproducible and testable CLI workflows

⸻

Contributing

Contributions are welcome, especially in:
	•	test fixtures and regression cases
	•	BAM support logic
	•	GFA context annotation
	•	XMFA parsing and corroboration
	•	report and plotting polish
	•	documentation and examples

Please open an issue before making large architectural changes.

⸻

Citation

If you use Panex Privus in published work, please cite the software release and associated manuscript once available.

A CITATION.cff file will be included in the repository.

⸻

License

See LICENSE.

⸻

Name

Panex Privus is the official project name.
privy is the command-line tool.

The project is named to reflect its central function: finding signal that belongs to the focal set and not to the rest.

---