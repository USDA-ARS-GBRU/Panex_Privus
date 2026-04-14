# Architecture

## Overview

Panex Privus is a comparative genomics toolkit centered on one inference problem:

> identify genomic signal shared within a target cohort and absent from off-target genomes

The software is designed to be:

- VCF-first
- interval-aware
- multi-evidence
- auditable
- scalable to plant pangenome workflows

This document describes the conceptual and software architecture of the package.

---

## Architectural priorities

### 1. Logic-centered rather than format-centered
Panex Privus supports multiple file formats, but no format is treated as the sole definition of truth.

Instead, all sources contribute evidence toward a common internal question:

- does this locus support target-private status?
- does it contradict target-private status?
- is the evidence ambiguous or missing?

This is the architectural core of the project.

### 2. VCF-first discovery
In v1, VCF is the primary discovery backend.

BAM, GFA, and XMFA are support and comparison layers. They provide corroboration, contradiction, or contextual annotation.

### 3. Missingness is explicit
Missing calls must not be silently folded into a generic pass/fail decision.

Panex Privus reports strictness classes so that biological support and technical incompleteness remain separable.

### 4. Streaming and chunked processing
The package must scale to plant pangenome datasets. Whole-file in-memory assumptions are not acceptable in the core analysis path.

---

## Command architecture

Panex Privus exposes four top-level commands:

privy scan
privy compare
privy report
privy plot

privy scan

Primary discovery engine.

Responsibilities:
	•	parse cohort definitions
	•	stream through primary input data
	•	identify target-private alleles or intervals
	•	classify strictness
	•	merge loci into candidate regions
	•	optionally annotate with BAM, GFA, and XMFA evidence
	•	emit machine-readable outputs and run metadata

privy compare

Cross-evidence reconciliation engine.

Responsibilities:
	•	compare loci or regions across sources
	•	quantify overlap and compatibility
	•	classify support, contradiction, or source specificity
	•	emit comparison tables and summaries

privy report

Interpretation engine.

Responsibilities:
	•	summarize discovery results
	•	rank loci and regions
	•	summarize strictness classes
	•	summarize evidence support and contradictions
	•	render Markdown and optional HTML reports

privy plot

Visualization engine.

Responsibilities:
	•	generate focused, publication-quality plots
	•	explain loci and regions
	•	visualize genotype patterns and support layers
	•	avoid becoming a general-purpose genome browser

⸻

Domain model

The internal architecture is built around a small set of domain objects.

CohortDefinition

Represents the biological grouping.

Fields:
	•	targets
	•	off_targets
	•	ignored_samples
	•	metadata

Responsibilities:
	•	validate cohort membership
	•	prevent overlap between target and off-target sets
	•	provide target/off-target lookup utilities

Locus

Represents a genomic site or interval under evaluation.

Fields:
	•	locus_id
	•	contig
	•	start
	•	end
	•	locus_type
	•	primary_source
	•	source_ids
	•	metadata

Responsibilities:
	•	coordinate normalization
	•	overlap logic
	•	merge behavior
	•	provenance tracking

AllelePattern

Represents cohort-level VCF logic for a candidate allele.

Fields:
	•	allele_key
	•	target_support_n
	•	target_total_n
	•	offtarget_support_n
	•	offtarget_total_n
	•	target_missing_n
	•	offtarget_missing_n
	•	strictness_class
	•	pattern_pass
	•	pattern_reason

Responsibilities:
	•	encode private-allele logic
	•	preserve missingness
	•	support downstream scoring

EvidenceRecord

Represents one normalized statement about a locus from any source.

Fields:
	•	locus_id
	•	source_type
	•	sample_id or group_id
	•	evidence_class
	•	metric_name
	•	metric_value
	•	qualifiers
	•	provenance

Responsibilities:
	•	normalize heterogeneous evidence
	•	preserve source traceability
	•	support compare and report steps

ComparisonRecord

Represents a cross-source comparison outcome.

Fields:
	•	locus_id
	•	source_a
	•	source_b
	•	match_class
	•	coordinate_overlap
	•	state_compatibility
	•	support_summary
	•	contradiction_summary
	•	comparison_score

ScoredHit

Represents a ranked final result.

Fields:
	•	locus_id
	•	discovery_score
	•	support_score
	•	penalty_score
	•	final_score
	•	rank
	•	strictness_class
	•	summary_label

⸻

Evidence model

The package uses a unified evidence model with four main evidence classes:
	•	support
	•	absence
	•	ambiguous
	•	contradiction
	•	uninformative

Different sources contribute different evidence types:

VCF

Primary discovery evidence.

VCF determines:
	•	target support
	•	off-target exclusion
	•	missingness
	•	candidate private alleles
	•	initial loci and intervals

BAM

Read-level support evidence.

BAM contributes:
	•	depth
	•	allele counts
	•	allele fraction
	•	optional soft-clip and split-read summaries
	•	absence-confidence context in off-target samples

GFA

Graph-context evidence.

GFA contributes:
	•	path membership
	•	local graph complexity
	•	branch and junction proximity
	•	structural-context annotation

XMFA

Alignment-corroboration evidence.

XMFA contributes:
	•	target-specific aligned patterns
	•	regional consistency
	•	gap-aware corroboration

⸻

Discovery architecture

Primary scan backend

The initial and primary backend is VCF-based.

Scan workflow
	1.	validate config and cohort definitions
	2.	open indexed VCF input
	3.	iterate by contig or chunk
	4.	evaluate alternate alleles for target-private patterns
	5.	classify strictness
	6.	emit locus hits
	7.	merge nearby loci into regions
	8.	annotate candidate loci or regions with optional BAM/GFA/XMFA evidence
	9.	score loci and regions
	10.	write outputs

Private allele logic

Default mode in v1 is private_allele.

A candidate allele passes when:
	•	target support meets threshold
	•	off-target support is zero
	•	filter criteria are met

Missingness is reported as a strictness class, not hidden.

⸻

Strictness classes

Strictness classes separate technical missingness from biological contradiction.

Expected classes:
	•	strict_complete
	•	strict_target_missing
	•	strict_offtarget_missing
	•	strict_both_missing
	•	relaxed_threshold
	•	contradicted

This distinction is central to the package and should remain visible in all outputs.

⸻

Region model

Single variants are often insufficient as biological objects. The package therefore supports interval construction.

Region construction

Nearby passing loci may be merged into candidate regions using configurable rules such as:
	•	merge distance
	•	variant-type compatibility
	•	allele/state consistency
	•	maximum gap without signal

Region summaries

Each region should summarize:
	•	number of constituent loci
	•	variant-type composition
	•	dominant strictness class
	•	target consistency
	•	off-target exclusion
	•	aggregate support score

⸻

Comparison architecture

privy compare is designed to compare sources using:
	•	interval overlap
	•	source-aware tolerances
	•	evidence compatibility

Match classes

Comparison outputs classify loci as:
	•	supported
	•	partially_supported
	•	contradicted
	•	source_specific
	•	uninformative
	•	missing_data

Compatibility dimensions

Comparisons may evaluate:
	•	coordinate overlap
	•	allele/state compatibility
	•	target support agreement
	•	off-target exclusion agreement
	•	boundary tolerance
	•	evidence sufficiency

This allows meaningful VCF vs BAM, VCF vs GFA, VCF vs XMFA, and multi-evidence comparisons.

⸻

Scoring architecture

Scoring is designed to be transparent and configurable.

Discovery score

Derived from the VCF pattern:
	•	target support fraction
	•	off-target exclusion fraction
	•	variant quality
	•	private-allele specificity

Support score

Derived from evidence overlays:
	•	BAM support
	•	GFA context support
	•	XMFA corroboration

Penalty score

Derived from:
	•	target missingness
	•	off-target missingness
	•	contradictory evidence
	•	low-confidence support
	•	highly ambiguous contexts

Final score

The final rank score is additive and user-configurable:

final_score = discovery_score + support_score - penalty_score

Scoring weights must be written into run.json.

⸻

Performance model

Panex Privus is designed for plant pangenome-scale workflows.

Required design constraints
	•	indexed file access
	•	streaming over primary records
	•	chunked contig/window processing
	•	localized evidence queries
	•	bounded memory usage
	•	support for manifest-based BAM inputs

Recommended execution phases

The package should separate:
	1.	discovery
	2.	support annotation
	3.	comparison
	4.	reporting

This improves:
	•	scalability
	•	testability
	•	debuggability
	•	reproducibility

⸻

Package layout

src/privy/
├── cli/
│   ├── main.py
│   ├── scan.py
│   ├── compare.py
│   ├── report.py
│   └── plot.py
├── core/
│   ├── cohort.py
│   ├── locus.py
│   ├── patterns.py
│   ├── evidence.py
│   ├── scoring.py
│   ├── intervals.py
│   └── config.py
├── io/
│   ├── vcf.py
│   ├── bam.py
│   ├── gfa.py
│   ├── xmfa.py
│   ├── bed.py
│   ├── tsv.py
│   └── jsonio.py
├── backends/
│   ├── vcf_scan.py
│   ├── bam_support.py
│   ├── gfa_support.py
│   └── xmfa_support.py
├── compare/
│   ├── engine.py
│   ├── overlap.py
│   ├── compatibility.py
│   └── classifiers.py
├── report/
│   ├── summary.py
│   ├── markdown.py
│   └── html.py
├── plot/
│   ├── loci.py
│   ├── regions.py
│   ├── summaries.py
│   └── themes.py
├── utils/
│   ├── logging.py
│   ├── validation.py
│   ├── parallel.py
│   ├── metrics.py
│   └── misc.py
└── __init__.py


⸻

Testing architecture

The project should include:

Unit tests
	•	cohort validation
	•	private-allele logic
	•	strictness classification
	•	interval merging
	•	overlap logic
	•	compatibility logic
	•	scoring

Integration tests
	•	VCF scan
	•	VCF + BAM support
	•	VCF + GFA annotation
	•	VCF + XMFA corroboration
	•	compare workflows

Regression tests

Curated fixtures for:
	•	missing data cases
	•	multiallelic records
	•	symbolic SVs
	•	contradictory off-target evidence
	•	region merge behavior

⸻

Configuration model

YAML config support is required for reproducible workflows.

Priority order:
	1.	package defaults
	2.	YAML config
	3.	CLI overrides

The resolved configuration must be written to run.json.

⸻

Non-goals for v1

Panex Privus is not intended to be:
	•	a variant caller
	•	an assembler
	•	a graph constructor
	•	a sequence extraction utility
	•	a genome browser replacement
	•	a generic workflow orchestrator

These boundaries protect clarity and keep the project coherent.

⸻

Summary

Panex Privus is designed as a logic-centered, VCF-first, multi-evidence comparative genomics framework for discovering and validating target-private genomic signal.

Its architectural center is simple:
	•	define the focal group
	•	define the background group
	•	detect what belongs only to the focal group
	•	compare evidence honestly
	•	preserve uncertainty instead of hiding it

---

# `privy --help` CLI contract

Panex Privus
CLI: privy

A comparative genomics toolkit for discovering target-private alleles and regions
shared within a focal cohort and absent from off-target genomes.

USAGE:
  privy [OPTIONS] COMMAND [ARGS]...

COMMANDS:
  scan      Discover target-private alleles or regions from VCF-first workflows
  compare   Compare loci or regions across VCF, BAM, GFA, and XMFA evidence
  report    Generate ranked summaries, QC tables, and Markdown/HTML reports
  plot      Create focused locus and region visualizations

GLOBAL OPTIONS:
  --config PATH              Path to YAML configuration file
  --project-name TEXT        Optional project name written into outputs
  --outdir PATH              Output directory
  --threads INTEGER          Number of worker threads to use where supported
  --log-level TEXT           Logging level: debug, info, warning, error
  --quiet                    Reduce console output
  --version                  Show version and exit
  -h, --help                 Show this message and exit


SCAN
  Discover target-private alleles and candidate private regions.

USAGE:
  privy scan [OPTIONS]

PRIMARY INPUT OPTIONS:
  --vcf PATH                 Indexed multisample VCF (.vcf.gz + .tbi); primary v1 backend
  --xmfa PATH                XMFA alignment file; optional secondary or alternate input
  --gfa PATH                 GFA graph file; optional graph-context support layer
  --bam PATH [PATH ...]      One or more BAM files mapped to the same reference
  --bam-manifest PATH        TSV manifest mapping BAM files to sample names/groups

COHORT OPTIONS:
  --targets TEXT [TEXT ...]      Target sample names
  --off-targets TEXT [TEXT ...]  Off-target sample names
  --ignore-samples TEXT [TEXT ...]
                                 Samples to ignore during discovery
  --cohort-file PATH             Optional cohort definition file

DISCOVERY OPTIONS:
  --mode TEXT                Discovery mode. Supported:
                               private_allele
                               private_genotype
                               private_sv_state
                             Default: private_allele
  --min-target-support FLOAT Minimum fraction of target samples supporting allele
  --max-off-target-support FLOAT
                             Maximum fraction of off-target samples supporting allele
  --allow-multiallelic / --no-allow-multiallelic
                             Whether to evaluate multiallelic records
  --pass-only / --no-pass-only
                             Require VCF FILTER=PASS
  --min-qual FLOAT           Minimum VCF QUAL
  --region TEXT              Restrict scan to region: contig:start-end
  --contig TEXT              Restrict scan to a contig
  --chunk-size INTEGER       Chunk size for streaming large contigs
  --merge-distance INTEGER   Merge nearby passing loci into candidate regions
  --same-variant-class-only / --no-same-variant-class-only
                             Only merge loci of the same variant class

STRICTNESS OPTIONS:
  --strictness-report / --no-strictness-report
                             Report strictness classes explicitly
  --relaxed-target-missing FLOAT
                             Optional threshold for tolerated target missingness
  --relaxed-offtarget-missing FLOAT
                             Optional threshold for tolerated off-target missingness

BAM SUPPORT OPTIONS:
  --bam-min-depth INTEGER    Minimum depth for BAM evidence evaluation
  --bam-min-alt-count INTEGER
                             Minimum alternate-supporting reads
  --bam-min-alt-fraction FLOAT
                             Minimum alternate allele fraction
  --summarize-softclips / --no-summarize-softclips
                             Summarize soft-clipped reads near candidate loci
  --summarize-splitreads / --no-summarize-splitreads
                             Summarize split-read support near candidate loci

GFA SUPPORT OPTIONS:
  --junction-window-bp INTEGER
                             Window around locus/region for branch-junction annotation
  --report-path-membership / --no-report-path-membership
                             Report GFA path membership where available
  --report-graph-complexity / --no-report-graph-complexity
                             Summarize local graph complexity

XMFA SUPPORT OPTIONS:
  --gap-aware / --no-gap-aware
                             Use gap-aware alignment corroboration
  --xmfa-window-bp INTEGER   Window for local XMFA corroboration

SCORING OPTIONS:
  --discovery-weight FLOAT   Weight for discovery score
  --support-weight FLOAT     Weight for support score
  --penalty-weight FLOAT     Weight for penalty score

OUTPUT OPTIONS:
  --write-hits / --no-write-hits
                             Write hits.tsv
  --write-regions / --no-write-regions
                             Write regions.tsv
  --write-evidence / --no-write-evidence
                             Write evidence.tsv
  --write-sample-support / --no-write-sample-support
                             Write sample_support.tsv
  --write-qc / --no-write-qc
                             Write qc.tsv
  --write-run-json / --no-write-run-json
                             Write run.json

SCAN OUTPUTS:
  hits.tsv
  regions.tsv
  evidence.tsv
  sample_support.tsv
  qc.tsv
  run.json

SCAN NOTES:
  - VCF is the primary discovery backend in v1
  - Missingness is reported via strictness_class
  - BAM, GFA, and XMFA provide support, contradiction, or context
  - Designed for indexed streaming and plant pangenome-scale workflows


COMPARE
  Compare loci or regions across evidence sources.

USAGE:
  privy compare [OPTIONS]

INPUT OPTIONS:
  --hits PATH                hits.tsv from privy scan
  --regions PATH             regions.tsv from privy scan
  --vcf PATH                 Indexed VCF for comparison
  --bam PATH [PATH ...]      BAM files for comparison
  --bam-manifest PATH        BAM manifest TSV
  --gfa PATH                 GFA graph file
  --xmfa PATH                XMFA alignment file
  --a PATH                   First comparison input
  --b PATH                   Second comparison input

COMPARE OPTIONS:
  --mode TEXT                Comparison mode. Supported:
                               vcf_vs_bam
                               vcf_vs_gfa
                               vcf_vs_xmfa
                               scan_vs_scan
                               multi_evidence
                             Default: multi_evidence
  --overlap-mode TEXT        Overlap mode:
                               any
                               reciprocal
                               contained
                             Default: reciprocal
  --min-reciprocal-overlap FLOAT
                             Minimum reciprocal overlap for interval matching
  --breakpoint-tolerance-bp INTEGER
                             Tolerance for breakpoint-aware comparisons
  --require-state-compatibility / --no-require-state-compatibility
                             Require allele/state compatibility in addition to overlap

OUTPUT OPTIONS:
  --write-compare-tsv / --no-write-compare-tsv
                             Write compare.tsv
  --write-summary-tsv / --no-write-summary-tsv
                             Write compare_summary.tsv
  --write-json / --no-write-json
                             Write compare.json

COMPARE CLASSES:
  supported
  partially_supported
  contradicted
  source_specific
  uninformative
  missing_data

COMPARE NOTES:
  Comparison is based on interval overlap plus evidence compatibility.
  Different evidence layers may support, contradict, or fail to inform a locus.


REPORT
  Generate ranked summaries and human-readable reports.

USAGE:
  privy report [OPTIONS]

INPUT OPTIONS:
  --hits PATH                hits.tsv
  --regions PATH             regions.tsv
  --evidence PATH            evidence.tsv
  --compare PATH             compare.tsv
  --qc PATH                  qc.tsv
  --run-json PATH            run.json

REPORT OPTIONS:
  --format TEXT              Report format:
                               markdown
                               html
                               both
                             Default: markdown
  --top-n INTEGER            Number of top loci/regions to summarize
  --include-qc / --no-include-qc
                             Include QC section
  --include-strictness / --no-include-strictness
                             Include strictness class summary
  --include-compare / --no-include-compare
                             Include compare summary
  --include-regions / --no-include-regions
                             Include candidate region summary
  --title TEXT               Optional report title

OUTPUTS:
  summary.tsv
  ranked_hits.tsv
  strictness_summary.tsv
  support_summary.tsv
  contradiction_summary.tsv
  report.md
  report.html (optional)

REPORT NOTES:
  The report command is designed to convert raw outputs into collaborator-ready summaries.


PLOT
  Generate focused plots for loci, regions, and summary diagnostics.

USAGE:
  privy plot [OPTIONS]

INPUT OPTIONS:
  --hits PATH                hits.tsv
  --regions PATH             regions.tsv
  --evidence PATH            evidence.tsv
  --vcf PATH                 Indexed VCF
  --bam PATH [PATH ...]      BAM files
  --bam-manifest PATH        BAM manifest TSV
  --gfa PATH                 GFA graph file
  --xmfa PATH                XMFA alignment file

SELECTION OPTIONS:
  --locus-id TEXT            Plot a specific locus
  --region-id TEXT           Plot a specific region
  --top-n INTEGER            Plot top N loci or regions by score
  --contig TEXT              Restrict plots to a contig
  --region TEXT              Restrict plots to contig:start-end

PLOT TYPES:
  --plot-type TEXT           Supported:
                               locus_panel
                               region_summary
                               genotype_heatmap
                               strictness_bar
                               support_bar
                               depth_panel
                             Default: locus_panel

PLOT OPTIONS:
  --width FLOAT              Figure width
  --height FLOAT             Figure height
  --dpi INTEGER              Figure DPI
  --output-format TEXT       png, pdf, svg
  --show-labels / --no-show-labels
                             Show sample or locus labels where applicable

PLOT NOTES:
  - Focused explanatory plots, not a genome browser
  - Designed for diagnostics and publication-ready figure generation


EXAMPLES

  Minimal scan:
    privy scan \
      --vcf cohort.vcf.gz \
      --targets Benning Harosoy Clark \
      --off-targets Jack Lee Minsoy \
      --mode private_allele \
      --outdir results/

  Scan with BAM and GFA support:
    privy scan \
      --vcf cohort.vcf.gz \
      --targets Benning Harosoy Clark \
      --off-targets Jack Lee Minsoy \
      --bam bam_manifest.tsv \
      --gfa graph.gfa.gz \
      --merge-distance 1000 \
      --outdir results/

  Compare evidence:
    privy compare \
      --hits results/hits.tsv \
      --vcf cohort.vcf.gz \
      --bam bam_manifest.tsv \
      --gfa graph.gfa.gz \
      --mode multi_evidence \
      --outdir compare/

  Generate report:
    privy report \
      --hits results/hits.tsv \
      --regions results/regions.tsv \
      --evidence results/evidence.tsv \
      --compare compare/compare.tsv \
      --qc results/qc.tsv \
      --outdir report/

  Plot top loci:
    privy plot \
      --hits results/hits.tsv \
      --top-n 10 \
      --vcf cohort.vcf.gz \
      --bam bam_manifest.tsv \
      --outdir plots/


⸻
