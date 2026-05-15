---
title: Architecture
description: Panex Privus architecture and CLI surface.
---

# Architecture

## Overview

Panex Privus is a comparative genomics toolkit centered on one inference problem:

> identify genomic signal shared within a target cohort and absent from off-target genomes

The software is designed to be:

- VCF/GFA primary-discovery capable
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

### 2. Explicit primary inputs
VCF and GFA are both primary discovery backends.

VCF discovers target-private alleles from genotype calls. GFA discovers
target-private graph segments from path/walk traversal. BAM is a support
layer for VCF hits: it adds read-level depth and allele-fraction evidence at
already discovered loci. XMFA is not part of the active v0.8 roadmap.

### 3. Missingness is explicit
Missing calls must not be silently folded into a generic pass/fail decision.

Panex Privus reports strictness classes so that biological support and technical incompleteness remain separable.

### 4. Streaming and chunked processing
The package must scale to plant pangenome datasets. Whole-file in-memory assumptions are not acceptable in the core analysis path.

---

## Command architecture

Panex Privus exposes nine top-level commands:

privy scan
privy compare
privy pangenome
privy landscape
privy report
privy plot
privy annotate
privy export
privy index

privy scan

Primary discovery engine.

Responsibilities:
	•	parse cohort definitions
	•	stream through primary input data
	•	identify target-private alleles or intervals
	•	classify strictness
	•	merge loci into candidate regions
	•	optionally annotate VCF hits with BAM evidence
	•	emit machine-readable outputs and run metadata

privy compare

Cross-evidence reconciliation engine.

Responsibilities:
	•	compare loci or regions across sources
	•	quantify overlap and compatibility
	•	classify support, contradiction, or source specificity
	•	emit comparison tables and summaries

privy pangenome

Whole-feature summary engine.

Responsibilities:
	•	turn GFA segments or VCF alternate alleles into a shared feature matrix
	•	summarize full, target, and off-target feature presence
	•	report core/accessory/private/absent composition
	•	build feature coverage histograms and pangenome growth curves
	•	emit pangenome tables and run metadata

privy landscape

Windowed VCF context engine.

Responsibilities:
	•	stream VCF records into fixed-record or base-pair windows
	•	report per-sample missingness, heterozygosity, non-reference burden,
	  rare/private ALT burden, and median genotype-class frequency
	•	summarize target/off-target window-level context
	•	compute pairwise local genotype similarity
	•	merge adjacent nearest-background assignments into local background blocks
	•	emit window tables, block tables, similarity summaries, and run metadata

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
	•	generate focused, publication-quality plots from existing output tables
	•	explain loci and regions from scan and compare outputs
	•	render landscape and pangenome plots after data generation
	•	visualize genotype patterns, support layers, and run-level summaries
	•	avoid becoming a general-purpose genome browser

privy annotate

Gene annotation engine.

Responsibilities:
	•	intersect hits with GFF3 annotations
	•	classify loci as CDS, UTR, exonic, intronic, or intergenic
	•	handle contig aliases between discovery and annotation references
	•	emit annotated hit tables and summary counts

privy export

Downstream-format export engine.

Responsibilities:
	•	convert scan TSV outputs into genome-tool-friendly interval files
	•	write BED and GFF3 files for hits and merged regions
	•	preserve strictness, variant class, score, and provenance details
	•	emit export metadata

privy index

Reusable index builder.

Responsibilities:
	•	build reusable sidecar indexes for expensive input parsing paths
	•	currently supports `privy index gfa`
	•	allow repeated GFA scans without repeating the full graph walk parse

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

Primary graph discovery evidence.

GFA contributes:
	•	path membership
	•	walk/path traversal by sample
	•	graph segment coordinates from SN/SO/LN tags
	•	private graph-node and graph-region candidates
	•	GFA-specific segment length and coordinate-coverage summaries
	•	missing-vs-absent classification at graph loci

Landscape

Windowed context evidence.

Landscape contributes:
	•	per-sample missingness and genotype-burden tracks
	•	target/off-target private ALT burden by window
	•	local sample-similarity matrices
	•	local background blocks based on nearest-neighbor similarity

Landscape outputs are exploratory context. They do not replace formal
population-genetic tools, QTL/genetic-map software, or local ancestry models.

⸻

Discovery architecture

Primary scan backends

The active primary backends are VCF and GFA. VCF scans evaluate alternate
alleles. GFA scans evaluate coordinate-tagged graph segments. GFA calls are
private graph-node evidence, not VCF-style ALT-allele calls.

Scan workflow
	1.	validate config and cohort definitions
	2.	open indexed VCF input or parse GFA input
	3.	iterate by contig, region, record, or graph segment
	4.	evaluate alternate alleles or graph segments for target-private patterns
	5.	classify strictness
	6.	emit locus hits
	7.	merge nearby loci into regions
	8.	annotate VCF candidate loci with optional BAM evidence
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

Landscape architecture

`privy landscape` is separate from discovery. It explains genomic context around
the VCF callset and around scan candidates, but it does not decide which loci
are target-private candidates.

Window modes
	•	record windows: a fixed number of VCF records per window
	•	base-pair windows: a fixed physical span per window

Default record windows keep the number of variants per window stable across
uneven variant density. Base-pair windows are easier to interpret on
chromosome-scale coordinate plots.

Core landscape outputs
	•	sample_windows.tsv
	•	windows.tsv
	•	background_blocks.tsv
	•	candidate_introgression_blocks.tsv
	•	similarity.tsv
	•	landscape.json

Core landscape figures
	•	missingness_heatmap
	•	private_burden_heatmap
	•	local_background_map
	•	similarity_cluster_map

Local background blocks

A local background block is a run of adjacent windows where a sample's nearest
genotypic neighbor stays the same and passes a similarity threshold. These
blocks are best interpreted as shared genomic background segments. A true
recombination map usually requires a formal cross or pedigree design and a
genetic-map model.

Candidate introgression blocks are derived from target-sample windows whose
nearest local background is an off-target sample. They are exploratory
donor-like intervals, not formal local ancestry calls.

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

This allows meaningful VCF-vs-GFA and scan-vs-scan comparisons.

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
│   ├── pangenome.py
│   ├── landscape.py
│   ├── report.py
│   ├── plot.py
│   ├── annotate.py
│   ├── export.py
│   └── index.py
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
│   ├── pangenome.py
│   └── landscape.py
├── pangenome/
│   ├── analysis.py
│   ├── gfa.py
│   ├── model.py
│   └── vcf.py
├── landscape/
│   └── vcf.py
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
│   ├── pangenome.py
│   ├── landscape.py
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
	•	GFA scan
	•	pangenome summaries
	•	landscape window summaries
	•	compare workflows
	•	report, plot, and annotate workflows

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

Panex Privus is designed as a logic-centered, VCF/GFA primary-discovery,
multi-evidence comparative genomics framework for discovering and validating
target-private genomic signal.

Its architectural center is simple:
	•	define the focal group
	•	define the background group
	•	detect what belongs only to the focal group
	•	compare evidence honestly
	•	preserve uncertainty instead of hiding it

---

# Current CLI Surface

This section is a compact architecture contract for the active command surface.
For exhaustive flags and examples, prefer `README.md` and the live `privy
--help` output.

Panex Privus
CLI: privy

A comparative genomics toolkit for discovering target-private alleles and regions
shared within a focal cohort and absent from off-target genomes.

USAGE:
  privy [OPTIONS] COMMAND [ARGS]...

COMMANDS:
  scan      Discover target-private alleles or graph segments from VCF or GFA
  compare   Reconcile two privy scan hits.tsv files
  pangenome Summarize full, target, and off-target pangenomes
  landscape Create VCF sliding-window landscapes and local background maps
  report    Generate ranked summaries, QC tables, and Markdown/HTML reports
  plot      Create plots from existing scan, landscape, or pangenome outputs
  annotate  Intersect private loci with GFF3 gene annotations
  export    Export scan hits and regions to downstream genome-tool formats
  index     Build reusable indexes for supported inputs

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
  --vcf PATH                 Indexed multisample VCF (.vcf.gz + .tbi/.csi)
  --gfa PATH                 GFA graph file (.gfa or .gfa.gz); primary backend when used without --vcf
  --bam PATH [PATH ...]      One or more BAM files mapped to the same reference
  --bam-manifest PATH        TSV manifest mapping BAM files to sample names/groups

COHORT OPTIONS:
  --targets TEXT [TEXT ...]      Target sample names
  --targets-file PATH            Text file with one target sample per line
  --off-targets TEXT [TEXT ...]  Off-target sample names
  --off-targets-file PATH        Text file with one off-target sample per line
  --ignore-samples TEXT [TEXT ...]
                                 Samples to ignore during discovery
  --ignore-samples-file PATH     Text file with one sample name to ignore per line
  --cohort-file PATH             Optional cohort definition file (TSV or YAML)

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
  --min-segment-length INT   Minimum GFA segment length to evaluate

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
  graph_segments.tsv  (GFA scans only)
  qc.tsv
  run.json

SCAN NOTES:
  - VCF and GFA are primary discovery backends
  - Missingness is reported via strictness_class
  - BAM provides read-level support/contradiction evidence for VCF hits
  - Designed for indexed streaming and plant pangenome-scale workflows


COMPARE
  Compare two privy scan result sets by coordinate overlap and state compatibility.

USAGE:
  privy compare [OPTIONS]

INPUT OPTIONS:
  --hits-a PATH              hits.tsv from the first scan run
  --hits-b PATH              hits.tsv from the second scan run
  --source-a TEXT            Optional display label for source A
  --source-b TEXT            Optional display label for source B

COMPARE OPTIONS:
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
  Generate plots from existing scan, landscape, or pangenome outputs.

USAGE:
  privy plot [OPTIONS]

INPUT OPTIONS:
  --plot-set TEXT            scan, landscape, or pangenome; default scan
  --input-dir PATH           Existing landscape or pangenome result directory
  --hits PATH                hits.tsv
  --regions PATH             regions.tsv
  --evidence PATH            evidence.tsv
  --compare PATH             compare.tsv

SELECTION OPTIONS:
  --top-n INTEGER            Plot top N loci or regions by score

PLOT TYPES:
  --plot-type TEXT           Supported:
                               locus_panel
                               strictness_bar
                               score_distribution
                               support_bar
                               compare_summary
                               all
                             Default: all

PLOT OPTIONS:
  --width FLOAT              Figure width
  --height FLOAT             Figure height
  --dpi INTEGER              Figure DPI
  --output-format TEXT       png, pdf, svg
  --plot-scope TEXT          Landscape only: chromosome, genome, or both;
                             default chromosome
  --contig TEXT              Landscape only: render one contig/chromosome
  --contigs TEXT             Landscape only: comma-separated contigs/chromosomes
  --show-labels / --no-show-labels
                             Show sample or locus labels where applicable

PLOT NOTES:
  - Focused explanatory plots, not a genome browser
  - scan plots use --hits and optional scan/compare tables
  - landscape and pangenome plots use --input-dir
  - Designed for diagnostics and publication-ready figure generation

LANDSCAPE
  Create target/off-target-aware VCF sliding-window landscapes.

USAGE:
  privy landscape [OPTIONS]

INPUT OPTIONS:
  --vcf PATH                 Multisample VCF or BCF

COHORT OPTIONS:
  --targets TEXT [TEXT ...]  Target sample names
  --targets-file PATH        Text file with one target sample per line
  --off-targets TEXT [TEXT ...]
                             Off-target sample names
  --off-targets-file PATH    Text file with one off-target sample per line
  --ignore-samples TEXT [TEXT ...]
                             Samples to exclude
  --ignore-samples-file PATH Text file with one sample name to ignore per line
  --cohort-file PATH         Optional cohort definition file (TSV or YAML)

WINDOW OPTIONS:
  --window-records INT       Records per fixed-record window; default 200
  --step-records INT         Records advanced per fixed-record step; default 50
  --window-bp INT            Use physical base-pair windows
  --step-bp INT              Physical base-pair step

LANDSCAPE OPTIONS:
  --rare-max-count INT       Carrier-count threshold for rare ALT burden
  --rare-max-freq FLOAT      Carrier-frequency threshold for rare ALT burden
  --min-background-similarity FLOAT
                             Minimum similarity for local background assignment
  --min-introgression-similarity FLOAT
                             Minimum similarity for candidate introgression
  --min-introgression-delta FLOAT
                             Minimum off-target advantage over nearest target
  --max-introgression-missing-rate FLOAT
                             Maximum missingness in candidate windows
  --min-introgression-windows INT
                             Minimum adjacent windows per candidate block
  --similarity-output TEXT   Pairwise similarity table mode: full, summary, none;
                             default full
  --vcf-engine TEXT          VCF parser: auto, pysam, cyvcf2
  --local-pca / --no-local-pca
                             Write or skip local PCA coordinate table
  --plot-format TEXT         png, svg, or pdf for immediate --plots
  --plots / --no-plots       Write or skip landscape figures during analysis;
                             default no plots

LANDSCAPE OUTPUTS:
  sample_windows.tsv
  windows.tsv
  background_blocks.tsv
  candidate_introgression_blocks.tsv
  similarity.tsv
  local_pca.tsv
  landscape.json
  plots/landscape_plot_index.tsv (from privy plot --plot-set landscape)
  plots/missingness_heatmap.<contig>.png
  plots/private_burden_heatmap.<contig>.png
  plots/local_background_map.<contig>.png
  plots/similarity_cluster_map.<contig>.png (when full similarity rows exist)

LANDSCAPE NOTES:
  - Complements discovery; does not replace privy scan
  - Local background blocks are exploratory shared-background segments
  - Formal recombination maps require cross/pedigree-aware modeling


EXAMPLES

  Minimal scan:
    privy scan \
      --vcf cohort.vcf.gz \
      --targets Benning Harosoy Clark \
      --off-targets Jack Lee Minsoy \
      --mode private_allele \
      --outdir results/

  Combined VCF and GFA scan with BAM support:
    privy scan \
      --vcf cohort.vcf.gz \
      --gfa graph.gfa.gz \
      --targets Benning Harosoy Clark \
      --off-targets Jack Lee Minsoy \
      --bam-manifest bam_manifest.tsv \
      --merge-distance 1000 \
      --outdir results/

  Compare evidence:
    privy compare \
      --hits-a results/vcf/hits.tsv \
      --hits-b results/gfa/hits.tsv \
      --outdir results/compare/

  Generate report:
    privy report \
      --hits results/vcf/hits.tsv \
      --regions results/vcf/regions.tsv \
      --evidence results/vcf/evidence.tsv \
      --compare results/compare/compare.tsv \
      --qc results/vcf/qc.tsv \
      --outdir report/

  Plot top loci:
    privy plot \
      --plot-set scan \
      --hits results/vcf/hits.tsv \
      --top-n 10 \
      --outdir plots/

  Plot an existing VCF landscape:
    privy plot \
      --plot-set landscape \
      --input-dir results/landscape/ \
      --output-format pdf

  Build a VCF landscape:
    privy landscape \
      --vcf cohort.vcf.gz \
      --targets Benning Harosoy Clark \
      --off-targets Jack Lee Minsoy \
      --window-records 200 \
      --step-records 50 \
      --outdir results/landscape/

  Export intervals:
    privy export \
      --hits results/vcf/hits.tsv \
      --regions results/vcf/regions.tsv \
      --outdir exported/


⸻
