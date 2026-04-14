# Changelog

All notable changes to Panex Privus are documented here.

Format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).
Versioning follows [Semantic Versioning](https://semver.org/).

---

## [Unreleased] — v0.1.0-dev

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

### Not yet implemented

- `privy report` — report generation
- `privy plot` — visualization
- BAM support layer (Phase 4)
- XMFA support layer (Phase 5)
- `privy compare` (Phase 5)
