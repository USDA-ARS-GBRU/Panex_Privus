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
