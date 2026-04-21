# Changelog

All notable changes to Panex Privus are documented here.

Format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).
Versioning follows [Semantic Versioning](https://semver.org/).

---

## [Unreleased] — v0.6.0-dev

### Added

**Phase 1 — Repository scaffold and core domain (2026-04-10)**

- Full `src/privy/` package layout with `src`-layout and `hatchling` build backend
- Core domain objects (fully implemented and tested):
  - `CohortDefinition` — immutable, validated target/off-target sample assignments
  - `Locus` — genomic interval with overlap, distance, and merge logic (0-based half-open)
  - `AllelePattern` + `StrictnessClass` + `classify_strictness()` — the central
    private-allele classification kernel
  - `EvidenceRecord`, `ComparisonRecord` — format-agnostic normalised evidence objects
  - `ScoredHit` + scoring functions — transparent additive scoring
  - `merge_loci_to_regions()`, `reciprocal_overlap()` — interval merging logic
  - `PrivyConfig` + all section models — Pydantic config with YAML loading
- CLI scaffold: `privy scan`, `privy compare`, `privy report`, `privy plot` — all
  options wired; unimplemented backends raise `NotImplementedError`
- TSV and JSON output writer utilities (`TsvWriter`, `write_run_json`)
- Unit tests for all six core modules (113 tests)
- Example configs: `configs/privy.yaml`, `configs/privy_minimal.yaml`
- GitHub Actions CI workflow (Python 3.10/3.11/3.12)

**Phase 2 — VCF scan backend (2026-04-10 to 2026-04-14)**

- `src/privy/io/vcf.py` — complete VCF reader:
  - `validate_vcf_index()`, `get_vcf_samples()`, `get_vcf_contigs()`
  - `stream_vcf_records()` — indexed streaming via pysam with correct
    generator/file-handle lifecycle
  - `extract_cohort_counts()` — per-record cohort counting
  - `classify_variant_type()`, `format_allele_key()`, genotype helpers
- `src/privy/utils/metrics.py` — `ScanStats` accumulator for `qc.tsv`
- `src/privy/backends/vcf_scan.py` — complete scan orchestrator:
  - `run_vcf_scan()` — 11-step pipeline: validate → stream → filter →
    classify → accumulate → score → rank → merge → write
  - `_scan_contig()` — streaming per-record loop (FILTER, QUAL, multiallelic,
    per-alt allele evaluation)
  - All output writers: `hits.tsv`, `regions.tsv`, `evidence.tsv`,
    `sample_support.tsv`, `qc.tsv`, `run.json`
- `tests/conftest.py` — shared fixtures including `indexed_vcf` (pysam-based
  bgzip+tabix, no PATH dependency)
- `tests/data/small_cohort.vcf` — 9-record synthetic test VCF covering all
  strictness classes and filter paths
- `tests/unit/test_vcf_io.py` — 43 tests for VCF helper functions
- `tests/integration/test_vcf_scan.py` — 27 end-to-end scan tests

**Publication prep (2026-04-14)**

- New user-facing `README.md` for beginner bioinformatics students
- `CONTRIBUTING.md`, `CITATION.cff`, `LICENSE` (MIT), `.gitignore`
- Developer notes moved to `docs/dev/`

**Maintenance pass — senior review refinement (2026-04-14)**

- Added `tests/integration/test_scan_cli.py` with command-level coverage for
  successful `privy scan` execution and incomplete-cohort failure handling
- Tightened the pysam typing boundary in `src/privy/io/vcf.py` with explicit
  structural types for records, samples, and genotypes
- Updated `src/privy/backends/vcf_scan.py` and `src/privy/cli/scan.py` to use
  clearer modern annotations and a cleaner typed scan path
- Tightened config validation tests to assert `pydantic.ValidationError`
  instead of broad `Exception`

### Fixed

- `extract_cohort_counts()`: replaced `sample not in samples` guard with
  `try/except KeyError` because pysam's `VariantRecordSamples.__contains__`
  raises `KeyError` for sample names absent from the VCF header instead of
  returning `False`
- Fixed CLI/backend drift in `privy scan` by aligning the `write_run_json`
  parameter name between `src/privy/cli/scan.py` and
  `src/privy/backends/vcf_scan.py`
- Fixed several `privy scan` flags that were accepted by the CLI but not
  previously applied to the resolved configuration

**Phase 3 — GFA standalone scan (2026-04-14)**

- `src/privy/io/gfa.py` — complete GFA1/1.1 parser:
  - Parses S, L, P, W, H lines into typed dataclasses
  - Builds four inverted indices: `segment_to_paths`, `segment_to_walks`,
    `sample_to_paths`, `sample_to_walks`
  - Position index from SN/SO/LN segment tags for overlap queries
  - `get_samples_traversing_segment()` — who traverses a given segment
  - `get_samples_present_at_locus()` — who has any coverage at a position
    (distinguishes *absent* from *missing*)
  - `extract_cohort_segment_counts()` — same six-tuple interface as the VCF
    function; feeds `classify_strictness()` directly
  - Path-name parsing: `SAMPLE#HAP#CONTIG` (minigraph-cactus/PGGB) and plain names
- `src/privy/backends/gfa_scan.py` — new standalone discovery engine:
  - `run_gfa_scan()` mirrors `run_vcf_scan()` in structure and output schema
  - Locus IDs use `GPX` prefix; `variant_type = "graph_region"`
  - `sample_traversal` map replaces genotype: `"traverses"` / `"absent"` / `"missing"`
  - Same six output files as the VCF backend — directly comparable
- `src/privy/core/config.py` — `GfaConfig` gains `min_segment_length` and
  `path_name_format` fields
- `src/privy/cli/scan.py` — `privy scan --gfa` now routes to `run_gfa_scan()`
  when no `--vcf` is provided; error message updated to list GFA as a valid input
- `tests/data/small_cohort.gfa` — GFA1.1 W-line fixture with 5 samples, 7 segments,
  and two bubbles covering `strict_complete` and `strict_target_missing`
- `tests/unit/test_gfa_io.py` — 53 unit tests across 9 test classes
- `tests/integration/test_gfa_scan.py` — 37 end-to-end scan tests

**Maintenance pass — scan CLI truthfulness and GFA ergonomics (2026-04-14)**

- `src/privy/cli/scan.py` now applies only explicitly provided CLI options to
  the resolved config, preserving the intended priority order of defaults →
  config file → command-line flags
- Implemented `--cohort-file` loading for both YAML and TSV cohort definitions
- Added `--min-segment-length` to the `privy scan` CLI and wired it to
  `cfg.gfa.min_segment_length`
- Added command-level tests for GFA scans, boolean override behaviour, and
  YAML/TSV cohort-file loading

**Phase 4 — privy report (2026-04-20)**

- `src/privy/report/summary.py` — complete `run_report()` implementation:
  - Reads `hits.tsv` (required) + optional `regions.tsv`, `evidence.tsv`,
    `qc.tsv`, `run.json`
  - Writes `summary.tsv`, `ranked_hits.tsv`, `strictness_summary.tsv`,
    `support_summary.tsv` (when evidence provided), `contradiction_summary.tsv`
  - Three-tier format selection: `markdown`, `html`, `both`
- `src/privy/report/markdown.py` — `render_markdown_report()`:
  - Run summary, filtering/QC, top-N hits, strictness distribution, regions,
    source support, contradiction summary, caveats sections
  - Pipe-table formatting with `|`-safe cell escaping
- `src/privy/report/html.py` — `render_html_report()`:
  - Converts `report.md` to `report.html` using the `Markdown` library
    (`tables` + `fenced_code` extensions)
  - Self-contained HTML with minimal inline CSS
- `src/privy/io/tsv.py` — three new column schemas for report outputs:
  - `RANKED_HITS_COLUMNS` — `["rank", *HITS_COLUMNS]`
  - `STRICTNESS_SUMMARY_COLUMNS`
  - `SUPPORT_SUMMARY_COLUMNS`
- `pyproject.toml` — added `Markdown>=3.4` runtime dep,
  `types-Markdown` dev dep; bumped version to `0.3.0.dev0`
- `tests/unit/test_report_summary.py` — 35 unit tests across 8 classes
- `tests/integration/test_report.py` — 43 integration tests across 8 classes
  including CLI command-level tests

**Phase 4 — privy report (2026-04-20)**

- `src/privy/report/summary.py` — complete `run_report()` implementation:
  - Reads `hits.tsv` (required) + optional `regions.tsv`, `evidence.tsv`,
    `qc.tsv`, `run.json`
  - Writes `summary.tsv`, `ranked_hits.tsv`, `strictness_summary.tsv`,
    `support_summary.tsv` (when evidence provided), `contradiction_summary.tsv`
  - Three-tier format selection: `markdown`, `html`, `both`
- `src/privy/report/markdown.py` — `render_markdown_report()`:
  - Run summary, filtering/QC, top-N hits, strictness distribution, regions,
    source support, contradiction summary, caveats sections
  - Pipe-table formatting with `|`-safe cell escaping
- `src/privy/report/html.py` — `render_html_report()`:
  - Converts `report.md` to `report.html` using the `Markdown` library
    (`tables` + `fenced_code` extensions)
  - Self-contained HTML with minimal inline CSS
- `src/privy/io/tsv.py` — three new column schemas for report outputs
- `pyproject.toml` — added `Markdown>=3.4` runtime dep,
  `types-Markdown` dev dep; bumped version to `0.3.0.dev0`
- `tests/unit/test_report_summary.py` — 35 unit tests across 8 classes
- `tests/integration/test_report.py` — 43 integration tests across 8 classes

**Phase 5 — BAM support layer (2026-04-20)**

- `src/privy/core/config.py` — added `min_mapq` and `min_baseq` to
  `BamConfig` (both default 20)
- `src/privy/io/bam.py` — complete BAM I/O implementation:
  - `validate_bam_index()` — checks for `.bai` / `.csi` index
  - `get_bam_sample_name()` — reads SM tag from @RG header
  - `load_bam_manifest()` — parses TSV with `bam_path`, `sample_id` columns;
    skips `#` comment lines
  - `query_position_depth()` — per-position depth via `count_coverage` with
    mapping-quality filter and unmapped/secondary/supplementary exclusion
  - `query_allele_counts_at_locus()` — SNP pileup returning
    `(ref_count, alt_count, other_count)`; depth-only for indels
- `src/privy/backends/bam_support.py` — complete annotation engine:
  - `HitLocusInfo` — lightweight locus descriptor (decouples from VCF internals)
  - `BamAnnotationResult` — bundles evidence records, per-locus support scores,
    per-(locus, sample) depth/allele-fraction metrics
  - `resolve_bam_sample_pairs()` — manifest > explicit paths; SM-tag fallback
    to filename stem
  - `annotate_loci_with_bam()` — per-(locus, sample) pileup and classification;
    UNINFORMATIVE values excluded from support-score mean
  - `_classify_bam_evidence()` — target: SUPPORT / AMBIGUOUS;
    off-target: ABSENCE / CONTRADICTION; low depth: UNINFORMATIVE; indel: UNINFORMATIVE
- `src/privy/backends/vcf_scan.py` — BAM integration:
  - `HitRecord` gains `ref_allele` and `alt_allele` fields
  - BAM support block runs between scan accumulation and scoring when
    `--bam` or `--bam-manifest` is provided
  - `_score_hit_records()` accepts optional `support_scores` dict
  - `_write_evidence_tsv()` appends BAM `EvidenceRecord` rows after VCF rows
  - `_write_sample_support_tsv()` populates `depth` and `allele_fraction`
    columns from BAM metrics (previously always `"NA"`)
- `src/privy/cli/scan.py` — added `--bam-min-mapq` and `--bam-min-baseq`
  options, wired to `BamConfig.min_mapq` / `min_baseq` via `_apply_cli_overrides`
- `tests/conftest.py` — four new BAM fixtures built with pysam (no samtools
  required for fixture creation): `bam_target_t1`, `bam_offtarget_o1`,
  `bam_offtarget_o1_with_alt`, `bam_low_depth_t2`
- `tests/unit/test_bam_io.py` — 31 unit tests across 5 classes
- `tests/integration/test_bam_support.py` — 34 integration tests across 6
  classes including full VCF + BAM end-to-end tests

**Phase 6 — privy compare (2026-04-21)**

- `src/privy/backends/compare.py` — complete comparison engine:
  - `HitsRow` — typed dataclass for hits.tsv rows
  - `load_hits_tsv()` — parse hits.tsv with full type coercion and error messages
  - `infer_source_label()` — auto-infers "vcf" / "gfa" from locus_id prefix (PPX / GPX)
  - `reciprocal_overlap_rows()` — intersection / union overlap for HitsRow pairs
  - `is_state_compatible()` — strictness-class compatibility check; respects
    `require_state_compatibility` for strict_* vs relaxed_threshold discrimination
  - `classify_match()` — SUPPORTED / PARTIALLY_SUPPORTED / CONTRADICTED /
    SOURCE_SPECIFIC with clear precedence rules
  - `compute_comparison_score()` — overlap-scaled numeric score in [0, 1]
  - `find_best_match()` — contig-indexed search with reciprocal-overlap primary and
    breakpoint-tolerance fallback
  - `run_compare()` — full pipeline: load → index → match → classify → write outputs
  - Writes `compare.tsv` (18 columns), `compare_summary.tsv`, `compare.json`
- `src/privy/cli/compare.py` — complete rewrite:
  - Focused on scan-vs-scan (two hits.tsv files via `--hits-a`/`--hits-b`)
  - Removed XMFA option and multi-mode switching (XMFA support dropped)
  - Fixed import: `from privy.backends.compare import run_compare`
  - `--source-a`/`--source-b` explicit label overrides
  - `--min-reciprocal-overlap`, `--breakpoint-tolerance-bp`,
    `--require-state-compatibility` wired to `CompareConfig`
- `src/privy/io/tsv.py` — `COMPARE_COLUMNS` expanded to 18 columns (adds
  `compare_id`, `locus_id_a`, `locus_id_b`, `contig`, `start_a/b`, `end_a/b`,
  `strictness_a/b`); `COMPARE_SUMMARY_COLUMNS` added
- `src/privy/cli/main.py` — removed XMFA mention from global help text;
  updated compare example usage
- `tests/unit/test_compare_engine.py` — 47 unit tests across 7 classes
- `tests/integration/test_compare.py` — 24 integration tests across 8 classes
  including full CLI command-level tests
- `pyproject.toml` — bumped version to `0.5.0.dev0`

**Phase 7 — privy plot (2026-04-21)**

- `src/privy/plot/themes.py` — expanded colour palette and rcParams:
  - `MATCH_COLOURS` for compare plot (matches all six MatchClass values)
  - `STRICTNESS_ORDER` and `MATCH_ORDER` for consistent plot-axis ordering
  - Publication-ready rcParams (grid, legend, spines, DPI)
- `src/privy/plot/summaries.py` — four diagnostic plot functions:
  - `plot_strictness_bar()` — horizontal bar chart of strictness class distribution
  - `plot_score_distribution()` — stacked histogram of `final_score` by strictness class
  - `plot_support_bar()` — stacked bar of evidence class counts by source type
  - `plot_compare_summary()` — horizontal bar of match class counts from compare.tsv
- `src/privy/plot/loci.py` — locus panel and main dispatcher:
  - `plot_locus_panel()` — ranked lollipop of top-N hits, coloured by strictness class
  - `run_plot()` — dispatches to all applicable plot functions based on `plot_type`;
    `"all"` (default) generates every plot for which inputs are present
- `src/privy/cli/plot.py` — complete rewrite:
  - Simplified to `--hits` (required) + `--evidence`, `--compare` optional inputs
  - `--plot-type all` default; XMFA and raw VCF/BAM/GFA options removed
  - Outputs echoed to stdout unless `--quiet`
- `tests/unit/test_plot.py` — 21 unit tests across 5 classes
- `tests/integration/test_plot_cli.py` — 16 integration tests across 4 classes
- 535 total tests passing
- `pyproject.toml` — bumped version to `0.6.0.dev0`

### Not yet implemented

- `privy annotate` — GFF3/BED feature intersection (v0.7)
- `privy export` — BED/VCF/GFF3 output layer (v0.8)
- Multi-cohort batch mode (v0.9)
