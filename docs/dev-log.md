# Panex Privus — Development Log

This log records what was built in each phase, the reasoning behind key decisions,
and the state of the codebase at each milestone.  It is written for someone
unfamiliar with the project who needs to understand not just *what* was built
but *why* each design choice was made.

---

## 2026-04-10 — Phase 1: Repository Scaffold

### What this project is

Panex Privus (`privy`) is a comparative genomics toolkit for discovering
*target-private genomic signal*: alleles and regions shared within a focal
cohort and absent from off-target genomes.  The prototypical use case is a
plant pangenome study where a researcher wants to find variants specific to
a target cultivar group that are absent from all other accessions.

The package was designed from the start to be published on GitHub, used by
real researchers, and eventually cited in a paper.  Code quality and design
integrity were treated as first-class requirements from day one.

### Pre-existing material

Before Phase 1, the repository contained:
- `README.md` — a production-quality description of the package
- `docs/architecture.md` — a detailed architectural specification
- `00b_notes_forClaude.txt`, `00c_notes_forClaude.txt` — master design documents
- `example_scripts/` — three working proof-of-concept Python scripts proving
  the core logic in VCF, BAM, and XMFA contexts

These documents served as the authoritative source of truth for all Phase 1
decisions.  No design was invented from scratch.

### What was built

Phase 1 created the full repository scaffold: 63 files across the complete
`src/privy/` package layout, plus tests, configs, and CI.

**Core domain objects (fully implemented, tested):**

| Module | Object | What it does |
|--------|--------|--------------|
| `core/cohort.py` | `CohortDefinition` | Immutable, validated container for target/off-target sample assignments |
| `core/locus.py` | `Locus` | Genomic interval with overlap, distance, and merge logic (0-based half-open) |
| `core/patterns.py` | `AllelePattern`, `StrictnessClass`, `classify_strictness()` | The central classification kernel |
| `core/evidence.py` | `EvidenceRecord`, `ComparisonRecord` | Format-agnostic normalised evidence objects |
| `core/scoring.py` | `ScoredHit`, score component functions | Transparent additive scoring: final = discovery + support − penalty |
| `core/intervals.py` | `merge_loci_to_regions()`, `reciprocal_overlap()` | Interval merging and overlap logic |
| `core/config.py` | `PrivyConfig` + all section models | Pydantic config with YAML loading and priority chain |

**CLI (fully wired, correct option specs):**
- `privy scan` — all ~40 options from the architecture contract
- `privy compare` — full cross-evidence comparison option set
- `privy report` — all reporting options
- `privy plot` — all visualization options

Global options (`--config`, `--outdir`, `--threads`, `--log-level`, `--quiet`, `--version`)
live in `cli/context.py` as a module-level singleton populated by the root
callback and read by all subcommands.

**Stubs (correct interfaces, explicit `NotImplementedError` — not fake):**
- `backends/vcf_scan.py` — says exactly what it will do in Phase 2
- All IO format readers, compare, report, and plot modules

**Tests:** 113 unit tests, all passing.  Tests cover every code path in the
six implemented core modules.

### Key design decisions made in Phase 1

**1. Logic-centered, not format-centered**

The most important architectural decision was to build a logic layer
(`src/privy/core/`) that is independent of file formats.  All evidence
sources (VCF, BAM, GFA, XMFA) map into common domain objects.  This means:
- The scoring and comparison logic has no knowledge of file formats
- Formats can be added without touching the core logic
- Unit tests can test biological correctness without any file I/O

**2. StrictnessClass — missingness is never silent**

The six-value `StrictnessClass` enum is the most important design choice in
the codebase.  It separates three conceptually distinct states:
- *Biological support*: the allele is genuinely present in targets and absent
  from off-targets (`strict_complete`)
- *Technical incompleteness*: the pattern is consistent but some samples lack
  calls (`strict_target_missing`, `strict_offtarget_missing`, `strict_both_missing`)
- *Policy*: the user's thresholds let it pass but strict logic doesn't
  (`relaxed_threshold`)
- *Contradiction*: the private-allele model actually fails (`contradicted`)

A tool that collapses these into pass/fail is not auditable.  This enum
ensures that every output row carries enough information for a researcher to
judge data quality independently.

**3. Pydantic for config with three-tier priority**

Config priority: package defaults → YAML file → CLI flags.  This is the
conventional bioinformatics priority order.  Pydantic ensures that all config
values are type-validated before any analysis runs.  The resolved config is
written to `run.json` so every run is reproducible from its output directory.

**4. Coordinate convention: 0-based half-open [start, end)**

All `Locus` coordinates are 0-based half-open, matching pysam and BED format.
VCF 1-based positions are converted on read.  This is explicit in docstrings
everywhere.  It prevents off-by-one errors at the boundary between formats.

**5. Stubs use `NotImplementedError` with informative messages**

Stubs are honest — they raise `NotImplementedError` with a message explaining
exactly what will be built and in which phase.  They are not silent no-ops
that produce empty outputs, which would mislead users into thinking the tool
ran successfully.

---

## 2026-04-10 — Phase 2: VCF Scan Backend

### Goals

Deliver the first real end-to-end `privy scan` execution:
- Read an indexed multisample VCF with pysam
- Evaluate each alternate allele for target-private status
- Classify strictness via the Phase 1 `classify_strictness()` kernel
- Merge passing loci into candidate regions
- Score and rank hits
- Write `hits.tsv`, `regions.tsv`, `sample_support.tsv`, `qc.tsv`, `run.json`

### What was built

**`src/privy/io/vcf.py` — complete VCF reader**

Implements:
- `validate_vcf_index()` — checks for `.tbi` or `.csi` alongside the VCF
- `get_vcf_samples()` — returns sample names from VCF header
- `get_vcf_contigs()` — returns contig names from VCF header
- `stream_vcf_records()` — generator that opens/closes pysam.VariantFile
  correctly (avoids the context-manager + yield antipattern)
- `is_missing_genotype()` — detects `./. ` and `.|.`
- `has_alt_allele()` — checks if a sample GT contains a specific alternate
- `classify_variant_type()` — SNP/indel/SV from REF/ALT lengths
- `extract_cohort_counts()` — the core per-record counting function that
  feeds directly into `classify_strictness()`

**`src/privy/backends/vcf_scan.py` — complete scan backend**

The `run_vcf_scan()` function implements the full 10-step workflow from
`docs/architecture.md`:
1. Validate inputs
2. Load VCF metadata (samples, contigs)
3. Validate cohort against VCF samples (warn about unrecognised samples)
4. Stream through VCF by contig, chunking is delegated to pysam's indexed
   fetch — no whole-contig loading
5. Per record: apply FILTER/QUAL filters, enumerate alt alleles, count
   cohort support, call `classify_strictness()`
6. Accumulate `HitRecord` objects for passing loci
7. Merge passing loci with `merge_loci_to_regions()`
8. Score using `compute_discovery_score()` and `compute_penalty_score()`
9. Rank with `rank_scored_hits()`
10. Write all output files

`HitRecord` is an internal dataclass (private to the backends module) that
bundles `Locus + AllelePattern + sample genotypes` before scoring.

**`src/privy/utils/metrics.py` — complete ScanStats**

`ScanStats` accumulates per-run counts and produces the `qc.tsv` rows.
Tracked metrics: records evaluated, skipped (filter/qual), alleles evaluated,
alleles passed/contradicted, loci emitted, regions emitted, strictness class
distribution.

### Key design decisions made in Phase 2

**1. Accumulate hits in memory, then score and write**

The scan produces `HitRecord` objects accumulated in a list.  Only after
all contigs are processed are hits scored, ranked, merged, and written.

Rationale: ranking requires all final_scores (to assign ranks 1…N).  Region
merging requires all passing loci.  The `HitRecord` list is much smaller
than the VCF — a genome with 10 million variants might produce tens of
thousands of hits, which is manageable in RAM.  For very large hit sets, a
future pass can stream-score per contig; for Phase 2, this is acceptable.

The VCF *scan loop itself* is streaming — only one VCF record is in memory
at a time.

**2. Locus IDs are sequential `PPX{n:08d}`**

Locus IDs are assigned as a sequential counter (`PPX00000001`, `PPX00000002`,
…) during the scan.  The allele coordinate information is in `allele_key`
(format: `contig:pos:ref:alt`, 1-based VCF POS) and the `Locus` object.

Alternatives considered:
- Hash-based IDs (non-reproducible across different samples in the cohort)
- Coordinate-based IDs (can be very long for large indels/SVs)
Sequential IDs are short, reproducible within a run, and sort naturally.

**3. `source_ids` carries `locus_id` for region reconstruction**

Each `Locus` created during the scan has its own `locus_id` stored in
`source_ids`.  After `merge_loci_to_regions()`, each merged region's
`source_ids` accumulates all constituent locus IDs.  This allows the
backend to look up which `HitRecord` objects contributed to each region
without maintaining a separate index.

**4. Cohort vs. VCF sample validation**

The backend warns (does not fail) when cohort samples are absent from the
VCF header.  Absent samples are treated as entirely missing at every locus,
which will be reflected in the strictness class.  This is the correct
bioinformatics behaviour: a missing sample is a known uncertainty, not an
error.

The backend *does* fail if no target samples at all appear in the VCF — that
is unrecoverable.

**5. `sample_support.tsv` in Phase 2**

`sample_support.tsv` is written in Phase 2 from the genotype information
captured during the scan.  Depth and allele fraction columns are empty
(`NA`) because BAM support is not yet implemented.  The columns exist so the
file schema is stable — Phase 3 (BAM support) will populate those fields
without changing the column structure.

**6. Integration test fixtures use pysam's bundled tabix**

bgzip and tabix are not guaranteed to be on PATH in all environments, but
pysam bundles `pysam.tabix_compress` and `pysam.tabix_index`.  Test VCF
fixtures are created, compressed, and indexed entirely from Python in a
`tmp_path` pytest fixture.  This makes integration tests self-contained.

**7. pysam `VariantRecordSamples.__contains__` raises on unknown names**

A subtle pysam quirk discovered during testing: the `in` operator on
`record.samples` raises `KeyError` rather than returning `False` for sample
names that are absent from the VCF header.  The original implementation used
`if sample not in samples:` which worked for correctly-named samples but
threw on cohort samples not in the VCF.

Fix: replace the `in` guard with a `try/except KeyError` block, fetching the
GT tuple directly and catching `KeyError` or `TypeError`.  This is now the
documented pattern for safe per-sample access in `extract_cohort_counts`.

### Phase 2 test infrastructure

**`tests/conftest.py`** — session-scoped `VCF_TEXT` constant and three
fixtures: `indexed_vcf` (bgzip+tabix via pysam), `small_cohort`
(`CohortDefinition` with T1/T2 targets, O1/O2/O3 off-targets), `default_cfg`
(default `PrivyConfig`).

**`tests/data/small_cohort.vcf`** — plain-text reference copy of the
9-record synthetic VCF.  Each record is annotated in the conftest comment
block explaining the strictness class it should produce.

**`tests/unit/test_vcf_io.py`** — 43 unit tests.  Pure-Python tests
cover `is_missing_genotype`, `has_alt_allele`, `classify_variant_type`, and
`format_allele_key` without pysam.  Pysam-backed tests verify
`get_vcf_samples`, `get_vcf_contigs`, `validate_vcf_index`, and
`extract_cohort_counts` against the `indexed_vcf` fixture.

**`tests/integration/test_vcf_scan.py`** — 27 end-to-end tests across six
test classes (`TestOutputFiles`, `TestHitsTsv`, `TestQcMetrics`,
`TestRegionsTsv`, `TestSampleSupportTsv`, `TestErrorHandling`).  Tests
verify the full scan pipeline: correct hit count (7), strictness class
distribution, filter/qual/multiallelic skipping, coordinate conventions
(0-based half-open in hits.tsv, 1-based POS in allele_key), region merging
(merge_distance=0 → 7 individual regions; merge_distance=500 → 1 merged
region), and error handling (missing VCF, no samples in VCF, unsupported
mode).

**`src/privy/cli/scan.py`** — `run_vcf_scan()` is now wrapped in
`try/except` for `FileNotFoundError`, `ValueError` (exit code 1), and
`NotImplementedError` (exit code 2) so the CLI surfaces actionable error
messages rather than bare tracebacks.

### Phase 2 result

183 tests, all passing.  The complete `privy scan` workflow is implemented
and tested end-to-end from CLI through VCF streaming, strictness
classification, scoring, region merging, and output writing.

---

## 2026-04-14 — Maintenance Pass: Senior Code Review Refinement

### Goals

Refine the Phase 2 codebase without changing its architecture or biological
intent:
- tighten the CLI/backend contract
- improve typing at the pysam boundary
- raise test quality around actual CLI usage
- make the touched code more trustworthy under static analysis

The explicit constraint for this pass was to preserve the existing design
rather than redesign the tool from scratch.

### Weak points identified before changes

**1. The `privy scan` CLI and backend had drifted apart**

`src/privy/cli/scan.py` called `run_vcf_scan(..., write_run_json=...)`,
while the backend parameter was named `write_run_json_flag`.  This did not
surface in the existing test suite because the tests exercised the backend
entry point directly rather than the Typer command.

This was a real correctness issue in the user-facing CLI, not just a style
problem.

**2. Strict typing was enabled, but the VCF boundary was not truly typed**

`src/privy/io/vcf.py` and `src/privy/backends/vcf_scan.py` treated pysam
records as generic `object` values and relied on scattered `type: ignore`
comments.  That allowed the code to run, but it undermined the value of
`mypy --strict` in the most important backend path in the codebase.

**3. Test coverage was stronger for the backend than for the CLI**

The integration suite covered `run_vcf_scan()` thoroughly, but it did not
cover the Typer command surface where users actually interact with the tool.
As a result, interface drift between CLI wiring and backend implementation
could go unnoticed.

**4. Some tests asserted overly broad exceptions**

Several config tests used `pytest.raises(Exception)`.  Those tests would pass
for many unintended failures and therefore provided a weaker signal than they
appeared to.

### What was changed

**`src/privy/cli/scan.py` — repaired CLI/backend contract**

The `run_vcf_scan()` signature was normalised so the backend now accepts
`write_run_json`, matching the CLI call site and the user-facing option name.
The touched annotations in the scan command were also updated to modern
`X | None` and `list[...]` style for consistency with the package's
Python 3.10+ target.

This was intentionally a small repair, not a CLI redesign.

**`src/privy/io/vcf.py` — explicit typed pysam boundary**

The VCF reader now defines a small structural typing layer:
- `Genotype`
- `VariantSampleCall`
- `VariantSamples`
- `VariantRecordLike`

These protocols describe exactly the subset of pysam behaviour that Privus
relies on.  `stream_vcf_records()` now performs the cast once, at the IO
boundary, rather than forcing downstream code to pretend a `VariantRecord`
is an untyped `object`.

This makes the trust boundary explicit: the file reader owns the pysam
adaptation, while the scan backend can operate on a real typed contract.

**`src/privy/backends/vcf_scan.py` — removed scattered typing escapes**

The scan backend was updated to consume the typed VCF interface directly.
This eliminated the main `type: ignore` noise around:
- `record.filter`
- `record.qual`
- `record.alts`
- `record.ref`
- `record.pos`
- `record.chrom`
- `record.samples`

Sample genotype capture for `sample_support.tsv` now uses the shared
`Genotype` alias, which makes the backend's internal data structures easier
to reason about and less dependent on implicit pysam details.

**Tests — CLI coverage and sharper exception assertions**

A new file, `tests/integration/test_scan_cli.py`, was added with two
end-to-end command-level tests:
- `privy scan` runs successfully on the synthetic indexed VCF fixture
- `privy scan` fails cleanly when the cohort definition is incomplete

`tests/unit/test_config.py` was also tightened so validation tests assert
`pydantic.ValidationError` instead of generic `Exception`.

### Key design decisions in this maintenance pass

**1. Fix the boundary, not every downstream use site**

The right place to handle pysam typing uncertainty is the IO layer, where
the external dependency first enters the codebase.  This preserves the
architecture from Phase 1 and Phase 2: format-specific concerns stay at the
edges, while the backend operates on normalised domain data.

**2. Add CLI tests rather than relying only on static analysis**

`mypy` correctly found the CLI/backend mismatch, but static analysis alone
should not be the only protection for a user-facing command.  Adding real
Typer integration tests ensures the command path itself is executable and
correct.

**3. Improve trustworthiness before broad cleanup**

This pass did not attempt a repository-wide Ruff modernization.  There is
still annotation and import cleanup debt in untouched compare/report/plot
modules.  The deliberate choice here was to make the scan path and its tests
more reliable first, rather than mixing targeted correctness work with a
large formatting-only sweep.

### Maintenance result

The full test suite now contains 185 tests, all passing.

`mypy src` passes across all source files, including the VCF scan backend.
The touched files also pass Ruff cleanly.  Broader Ruff cleanup remains as
follow-on work outside this maintenance pass, but the most critical user path
(`privy scan`) is now both runtime-tested and statically coherent.

---

## Phase 3 — GFA standalone scan (2026-04-14)

### Design goal

GFA (Graphical Fragment Assembly) is a first-class primary input, not an annotation
layer.  When a user runs `privy scan --gfa`, the scanner discovers target-private graph
segments and writes the same six output files as the VCF backend.  The outputs are
directly comparable via `privy compare`.

### Discovery model

A graph segment is "target-private" when:
- Target-sample paths/walks traverse it, and
- Off-target paths/walks do not traverse it (or traverse a different bubble arm).

The StrictnessClass framework applies unchanged:
- `strict_complete` — all targets traverse; all off-targets are present at the locus
  but traverse an alternative segment
- `strict_target_missing` — off-target exclusion holds; one or more targets have no
  walk/path coverage at the locus at all
- `strict_offtarget_missing` — target support holds; one or more off-targets have no
  coverage
- `strict_both_missing` — pattern consistent with private status but missingness in both
- `contradicted` — off-targets also traverse this segment above the threshold

The key insight: **missing vs. absent** requires separate detection.
- *Absent* (informative): sample has coverage at the locus but traverses a different
  segment.  Counted as "called, does not support."
- *Missing* (uninformative): sample has no walk or path covering the locus at all.
  Requires coordinate information to detect.

### Coordinate requirement

Panex Privus can only place a segment on the output coordinate grid if it carries
`SN:Z:`, `SO:i:`, and `LN:i:` optional tags (minigraph/PGGB standard output).
Segments without these tags are silently skipped.  W-line coordinates are used
for missingness detection: if a sample has a W-line whose `seq_start`/`seq_end`
interval overlaps the locus, they are "present"; otherwise "missing."

### Missingness detection logic (implemented)

`get_samples_present_at_locus(graph, contig, start, end)`:
1. Check all W-lines: if `walk.seq_id == contig and walk.seq_start < end and walk.seq_end > start` → present
2. Check P-lines: if any traversed segment overlaps `[start, end)` via the position index → present

`extract_cohort_segment_counts()` then classifies each cohort sample as:
- **support** — traverses the segment
- **absent** (present at locus, not counted as missing) — present but on a different arm
- **missing** — not present at the locus at all

### Files created

**`src/privy/io/gfa.py`** — full GFA1/1.1 parser replacing the Phase 4 stub:
- Typed dataclasses: `GfaSegment`, `GfaLink`, `GfaPath`, `GfaWalk`, `GfaWalkStep`, `GfaGraph`
- `parse_gfa(path)` — reads entire file, builds all indices
- Four inverted indices built on parse: `segment_to_paths`, `segment_to_walks`,
  `sample_to_paths`, `sample_to_walks`
- Position index: `_contig_segments[contig]` = sorted `(start, end, seg_name)` tuples
- Public query functions: `query_segments_at_locus()`, `get_samples_traversing_segment()`,
  `get_samples_present_at_locus()`, `extract_cohort_segment_counts()`
- Path-name conventions: `SAMPLE#HAP#CONTIG` (pangenome tools) and plain names

**`src/privy/backends/gfa_scan.py`** — new standalone scan backend:
- `GfaHitRecord` — mirrors `HitRecord` but carries `segment_name` and `sample_traversal`
  (traverses/absent/missing map) instead of VCF fields
- `run_gfa_scan()` — same 7-step pipeline as `run_vcf_scan()`:
  parse → validate → filter → scan segments → score → merge → write
- `_scan_segments()` — iterates coordinate-indexed segments; applies `min_segment_length`
  filter; calls `extract_cohort_segment_counts()` → `build_allele_pattern()` → `classify_strictness()`
- All six output writers; `variant_type = "graph_region"`, locus IDs `GPX########`

**`src/privy/core/config.py`** — `GfaConfig` extended:
- `min_segment_length: int = 1` — filter very short bubble arms
- `path_name_format: str = "pangenome"` — documents convention in use

**`src/privy/cli/scan.py`** — routing updated:
- `if vcf is not None` → `run_vcf_scan()` (unchanged)
- `elif gfa is not None` → `run_gfa_scan()` (new branch)
- `else` → `NotImplementedError` (XMFA not yet implemented)
- Error message updated: "At least one primary input is required: --vcf, --gfa, or --xmfa"

**`tests/data/small_cohort.gfa`** — GFA1.1 W-line fixture:
- 5 samples: T1, T2 (targets), O1, O2, O3 (off-targets)
- 7 segments with SN/SO/LN tags on chr1
- 2 bubbles: chr1:8-18 (`strict_complete`) and chr1:60-67 (`strict_target_missing`;
  T2 has no walk covering this region)
- 3 backbone segments (all samples traverse → `contradicted`, never hit)

**`tests/unit/test_gfa_io.py`** — 53 tests across 9 classes:
`TestParsing`, `TestParsePLine`, `TestInvertedIndices`, `TestGetGfaSamples`,
`TestQuerySegmentsAtLocus`, `TestGetSamplesTraversingSegment`,
`TestGetSamplesPresentAtLocus`, `TestExtractCohortSegmentCounts`, `TestParseErrors`

**`tests/integration/test_gfa_scan.py`** — 37 tests across 8 classes:
`TestOutputFiles`, `TestHitsTsv`, `TestRegionsTsv`, `TestQcMetrics`,
`TestSampleSupportTsv`, `TestMinSegmentLength`, `TestContigFilter`, `TestErrorHandling`

### Phase 3 result

Total test suite: **280 tests, all passing.**

The GFA and VCF backends are architecturally symmetric: same StrictnessClass logic,
same scoring pipeline, same output schema.  A biologist can run `privy scan --vcf` and
`privy scan --gfa` independently and compare the two result sets with `privy compare`
once that backend is implemented.

---

## 2026-04-14 — Maintenance Pass: Scan CLI Truthfulness and GFA Ergonomics

### Goals

Refine the user-facing `privy scan` command after the standalone GFA backend
landed, without changing the architecture:
- make CLI flags truthfully affect the resolved config
- expose primary-GFA controls more clearly
- turn placeholder CLI affordances into real functionality
- add command-level tests for the new GFA workflow

This was a refinement pass, not a redesign.

### Weak points identified before changes

**1. Several `privy scan` options were accepted but not applied**

The command surface had grown enough that many flags were parsed by Typer and
shown in `--help`, but only a subset were actually copied into the resolved
`PrivyConfig`.  This is especially risky in a scientific CLI because a user
can reasonably assume that a flag appearing in help output also changes the
run recorded in `run.json`.

**2. GFA was now a primary backend, but the CLI still read partly like a support layer**

The GFA backend itself was functional, but the `--gfa` help text and key scan
options did not yet fully reflect that users can run `privy scan --gfa ...`
as a first-class primary analysis path.

**3. `--cohort-file` existed but had no real behaviour**

The flag appeared in the command interface and architecture docs, but it was
not used when building the effective cohort.  A dormant flag is worse than no
flag, because it creates false confidence in reproducibility.

### What was changed

**`src/privy/cli/scan.py` — explicit CLI override application**

The scan command now centralises config resolution through helper functions:
- `_apply_cli_overrides()`
- `_provided_updates()`
- `_was_provided()`

These use Click parameter-source inspection so only values explicitly supplied
on the command line override the YAML/default config.  This preserves the
intended priority model:

package defaults → config file → explicit CLI flags

This is more correct than blindly copying every Typer option value into the
config, because boolean defaults such as `--pass-only` or
`--report-graph-complexity` should not silently overwrite YAML settings unless
the user actually passed them.

**`src/privy/cli/scan.py` — real cohort-file support**

`--cohort-file` is now implemented for both:
- YAML (`targets`, `off_targets`, optional `ignored_samples`)
- TSV (`sample_id`/`sample`, `cohort_role`/`role`)

CLI cohort flags still take precedence over the cohort file, and the cohort
file takes precedence over the config file.  This keeps behaviour predictable
and matches the documented configuration philosophy of the project.

**`src/privy/cli/scan.py` — primary-GFA ergonomics**

The scan command now exposes:
- clearer `--gfa` help text ("primary backend when used without `--vcf`")
- `--min-segment-length` wired through to `cfg.gfa.min_segment_length`

This makes the standalone GFA workflow more self-sufficient from the command
line and reduces the need to edit YAML for a common graph-specific filter.

**`src/privy/backends/gfa_scan.py` — small structural cleanup**

Removed a dead `skipped_no_coords` counter and an unused `GfaSegment` import.
These were minor, but they made the backend look less deliberate than it is.

### Tests added in this pass

`tests/integration/test_scan_cli.py` now covers:
- successful `privy scan --gfa ...`
- boolean scan-option override reflected in `run.json`
- `--min-segment-length` affecting the GFA scan output
- cohort loading from YAML via `--cohort-file`
- cohort loading from TSV via `--cohort-file`

This matters because the scan CLI is now the place where multiple primary
backends, config layers, and cohort-definition pathways meet.  Backend tests
alone are not sufficient protection for that boundary.

### Maintenance result

The full suite now contains **280 tests, all passing**.

`mypy src` remains clean across all source files.  The touched scan/GFA files
also pass Ruff cleanly.  Broader Ruff modernization work still remains in
untouched compare/report/plot modules, but the main discovery command now
better matches what it claims to do in its help text and recorded config.

---

## 2026-04-20 — Phase 4: privy report (v0.3)

### Design goal

`privy report` converts raw scan outputs into a ranked, interpretable summary
package that a researcher can hand directly to a collaborator.  It operates
entirely on the TSV and JSON files already produced by `privy scan` — it does
not re-open the VCF or GFA, and it makes no biological decisions.  Its job is
interpretation and presentation, not discovery.

### Inputs and outputs

`privy report` reads:

| Input | Required | Purpose |
|-------|:--------:|---------|
| `hits.tsv` | Yes | Primary data source — all loci and scores |
| `regions.tsv` | No | Candidate region table |
| `evidence.tsv` | No | Per-locus evidence for source support summary |
| `qc.tsv` | No | Run-level metrics for the QC section |
| `run.json` | No | Provenance and parameters for the summary |

It writes:

| Output | Contents |
|--------|----------|
| `summary.tsv` | Run-level key/value table |
| `ranked_hits.tsv` | Top-N hits with explicit `rank` column |
| `strictness_summary.tsv` | Count and % per strictness class |
| `support_summary.tsv` | Evidence grouped by source × class (when evidence provided) |
| `contradiction_summary.tsv` | Contradiction metrics from QC/compare |
| `report.md` | Human-readable Markdown report |
| `report.html` | HTML version (`--format html` or `--format both`) |

### Architecture

**`src/privy/report/summary.py`** — `run_report()` orchestrator.

The function:
1. Loads all available input TSVs with `read_tsv()`.
2. Computes five structured summaries from the loaded data.
3. Writes the summary TSVs.
4. Assembles a `sections: dict[str, Any]` bundle.
5. Calls `render_markdown_report(sections, title, outdir)`.
6. Conditionally calls `render_html_report(md_path, outdir)`.

The five summary computations are isolated pure functions:
- `_rank_hits(rows, top_n)` — sort by `final_score` descending, take top N.
- `_compute_strictness_summary(rows)` — Counter over `strictness_class` column,
  with a canonical ordering (strict_complete first) and safe zero-filling.
- `_compute_support_summary(evidence_rows)` — Counter over
  `(source_type, evidence_class)` pairs from `evidence.tsv`.
- `_compute_contradiction_summary(qc_rows, compare_rows)` — extracts
  `alleles_contradicted` from `qc.tsv`; if `compare.tsv` is provided, also
  counts `match_class == "contradicted"` rows.
- `_compute_run_summary(hit_rows, region_rows, qc_rows, run_meta, cfg)` —
  assembles the summary.tsv rows from all available sources.

**`src/privy/report/markdown.py`** — `render_markdown_report()`.

Builds the report as a list of string lines, then joins and writes.  Each
section is a small private function (`_section_run_summary`, `_section_qc`,
etc.) that appends Markdown to the list.  Sections that require optional data
(QC, regions, support) are gated by checking `sections.get(key) is not None`.

The `_md_table()` helper produces standard pipe tables with `|`-escaped cell
values so allele keys (which contain `:` but not `|`) are safe.

**`src/privy/report/html.py`** — `render_html_report()`.

Reads `report.md`, converts using `markdown.markdown(text, extensions=["tables",
"fenced_code"])`, wraps in a minimal self-contained HTML template with inline
CSS, and writes `report.html`.  The `markdown` library is imported lazily
inside the function body (not at module top) so that test files that don't
trigger HTML rendering don't incur the import cost.

### New column schemas in `io/tsv.py`

Three column schemas were added:
- `RANKED_HITS_COLUMNS = ["rank", *HITS_COLUMNS]` — adds an explicit rank
  column to the existing hits schema.  Rank is 1-based and is assigned by the
  report module, not by the scan.
- `STRICTNESS_SUMMARY_COLUMNS = ["strictness_class", "n_loci", "pct_hits"]`
- `SUPPORT_SUMMARY_COLUMNS = ["source_type", "evidence_class", "n_records",
  "pct_of_source"]`

`summary.tsv` and `contradiction_summary.tsv` reuse the existing `QC_COLUMNS`
schema (`["metric", "value", "description"]`).

### Dependency

`Markdown>=3.4` was added to runtime deps.  It is pure Python, has no C
extensions, and brings the `tables` and `fenced_code` built-in extensions
needed for the HTML renderer.  `types-Markdown` was added to dev deps for
mypy stub completeness.

### Key design decisions

**1. Report operates on TSV files, not in-memory scan objects.**

The report command is a separate pass from the scan.  It reads `hits.tsv`
from disk rather than receiving in-memory `HitRecord` objects.  This means:
- The user can run `privy report` at any time after a scan, on any machine.
- No coupling between the report and the scan implementation details.
- The report can be re-run with different `--top-n` or `--format` without
  repeating the full scan.

**2. ranked_hits.tsv adds an explicit `rank` column.**

`hits.tsv` from the scan is already sorted by `final_score`, but the rank is
implicit.  `ranked_hits.tsv` makes it explicit (1, 2, 3...) so users can sort
the file by other columns (e.g. coordinate) and still recover the original rank.

**3. Sections are gated, not omitted by emptiness.**

Sections like QC and regions only appear in the report if the corresponding
input file was provided.  If `--qc` is not given, the Filtering and QC section
is simply absent — not replaced with an "empty" table.  This gives the user
explicit control and avoids misleading empty sections.

**4. HTML conversion uses the `Markdown` library, not string templating.**

The `Markdown` library converts pipe tables to `<table>` elements correctly.
A hand-rolled HTML template from the same section data would be a second code
path to maintain and diverge.  By deriving HTML from the Markdown source, we
guarantee that both outputs contain the same content.

### Phase 4 result

**358 tests, all passing.**  The full `privy report` workflow is implemented
and tested end-to-end from CLI through TSV loading, summary computation,
TSV output writing, Markdown rendering, and HTML conversion.

`privy report` is now a usable command.  A researcher can produce a shareable
summary from a completed scan in one line:

```bash
privy report \
  --hits results/hits.tsv \
  --regions results/regions.tsv \
  --qc results/qc.tsv \
  --format both \
  --outdir report/
```
