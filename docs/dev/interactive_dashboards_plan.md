---
title: Interactive Dashboard Development Plan
description: Working plan for adding self-contained interactive HTML dashboards to Panex Privus.
---

# Privy Interactive Dashboard Development Plan

This is a working design document for adding a `privy interactive` command to
Panex Privus. The first deliverable should be the focal interval browser
prototyped in the Gm15 Harosoy-sharp/Kingawa analysis; landscape, scan, and
pangenome dashboards can follow the same command namespace once the shared HTML
machinery is stable.

## Implementation Status

Started in development:

- `privy interactive --focus REGION` CLI registered.
- One self-contained HTML file is written per focus region.
- Multi-region runs also write `index.html`.
- Direct `--vcf` extraction writes one focal `*.sites.tsv` per focus region.
- A precomputed focal genotype `--sites-tsv` is still supported for reproducible
  rebuilds and debugging.
- Gene GFF3, exon/CDS/intron/promoter rendering, optional generic GFF3 tracks,
  feature-ranking TSVs, and run metadata are implemented for the first focus
  browser slice.
- Functional TSV joins and user-defined phenotype keyword groups are
  implemented for variant-supported feature lists.

Next implementation step: add more of the original Gm15 dashboard polish to the
package path, especially optional landscape or scan overlays and browser-side
refinements. Curated preset keyword packs may become a later convenience layer.

## Proposed CLI Namespace

Recommended first-pass command:

```bash
privy interactive \
  --focus Gm15:1-4000000 \
  --focus Gm12:2340000-2440000 \
  --vcf introg01.full.vcf.gz \
  --gff3 Wm82.gene_exons.gff3.gz \
  --samples Harosoy Harosoy-sharp Kingawa \
  --track-gff RepeatMasker=Wm82.repeats.gff3.gz \
  --track-gff SSR=Wm82.ssr.gff3 \
  --outdir interactive/
```

`--focus` is preferred over `--range` for this feature because it signals a
biological interval that will be interpreted and annotated, not merely a table
filter. `--range` is technically accurate, but it sounds generic and does not
carry the report/dashboard intent. `--region` is the common bioinformatics term
and could be supported as an alias later, but `--focus` is a good Panex Privus
term for user-facing interactive review.

Possible future shortcuts:

```bash
privy interactive --landscape landscape_outdir/ --outdir interactive/
privy interactive --scan scan_outdir/ --outdir interactive/
privy interactive --pangenome pangenome_outdir/ --outdir interactive/
```

Longer-term, if the command becomes crowded, the same model can become
subcommands without breaking the concept:

```bash
privy interactive focus --region Gm15:1-4000000
privy interactive landscape --input-dir landscape_outdir/
privy interactive scan --input-dir scan_outdir/
privy interactive pangenome --input-dir pangenome_outdir/
```

For now, a top-level `privy interactive --focus ...` command feels closest to
the workflow researchers will actually type.

## First Deliverable: `--focus`

The first implementation should reproduce the self-contained focal interval
browser as a Panex Privus command. It should handle one or many intervals and
write outputs that are easy to email, upload as supplementary material, or send
to co-authors.

### MVP User Stories

1. A researcher can run `privy interactive --focus CHR:START-END ...` and get
   an HTML browser for each focal interval.
2. The browser can be opened directly as a static file, with no local server,
   CDN, internet access, or Python runtime.
3. The browser displays gene models, exons, CDS, introns, strand-aware promoter
   windows, focal variants, and optional generic GFF3 tracks.
4. The browser distinguishes target-private SNPs, INDEL/complex records, and
   size/symbol-based SV-like records.
5. Candidate feature lists are variant-supported, not annotation-only.
6. Repeated `--focus` intervals produce either one indexed dashboard with tabs
   or one index page plus one HTML per focus region.
7. Every output has reproducibility metadata: command, inputs, sample roles,
   filters, interval definitions, and display assumptions.

### Initial Command Shape

```bash
privy interactive \
  --focus CONTIG:START-END \
  --focus CONTIG:START-END \
  --vcf cohort.vcf.gz \
  --gff3 genes.gff3.gz \
  --samples OFFTARGET DERIVED DONOR \
  --track-gff LABEL=features.gff3.gz \
  --functional-tsv annotations.tsv \
  --sample-abbrev HS=Harosoy-sharp \
  --promoter-bp 2000 \
  --sv-size-threshold 50 \
  --outdir interactive/
```

The first implementation can keep the three-sample `--samples OFFTARGET DERIVED
DONOR` model because it directly supports the donor/recurrent/derived
introgression use case that motivated this feature. A later generalization can
also accept normal Privy grouped cohort arguments:

```bash
--targets SAMPLE [SAMPLE ...] --off-targets SAMPLE [SAMPLE ...]
```

That generalization should map naturally onto target-private group signal but
may need separate UI language from the donor-like introgression browser.

### Useful Optional Inputs

- `--sites-tsv`: use a precomputed focal genotype TSV instead of extracting from
  VCF. Useful for debugging, reproducibility, and publication rebuilds.
- `--landscape-dir`: overlay windows, pairwise similarity, and candidate blocks
  from a `privy landscape` output directory.
- `--scan-dir`: overlay VCF/GFA scan hits, regions, strictness classes, and
  evidence summaries.
- `--report-html`: embed an existing report section below the browser.
- `--keyword-group NAME=term1,term2,...`: define phenotype-oriented candidate
  lists.
- `--title` and `--subtitle`: control collaborator-facing document text.

### Output Layout

For one interval:

```text
interactive/
  focus_Gm15_1_4000000.html
  focus_Gm15_1_4000000.features.tsv
  focus_Gm15_1_4000000.sites.tsv
  focus_Gm15_1_4000000.json
```

For multiple intervals:

```text
interactive/
  index.html
  focus_Gm15_1_4000000.html
  focus_Gm15_1_4000000.features.tsv
  focus_Gm12_2340000_2440000.html
  focus_Gm12_2340000_2440000.features.tsv
  interactive.json
```

Decision: default behavior is one static HTML file per `--focus` region. If
the run contains more than one focus region, also write an `index.html` that
links to each region dashboard and its audit files. Many expected use cases
will focus on a single region, so the single-focus path should feel like the
primary workflow rather than a special case. For novice-facing documentation,
recommend starting with a focus region around 4 Mbp or smaller; larger regions
can be tractable when variant and annotation density are low, but static HTML
size and browser responsiveness become data-dependent.

## Dashboard Design: Focus Browser

### First Screen

The first viewport should be the actual browser, not a landing page. The header
can be compact:

- project title
- region selector when multiple focus regions are present
- region coordinates
- key counts: genes, target-private variants, added GFF3 features

### Browser Tracks

Recommended vertical order:

1. coordinate ruler
2. optional donor-like block and breakpoint shading
3. target-private variant density
4. PASS SNP point layer
5. INDEL/complex point layer
6. SV-like point layer
7. genes with exon/CDS/intron structure
8. strand-aware promoters
9. added GFF3 tracks, one row per `--track-gff` label
10. optional landscape/scan overlays when provided

The density and variant tracks should stay compact. The browser is for scanning
relationships among evidence layers, so annotation tracks should not be pushed
far below the fold.

### Interactions

- zoom presets: whole focus, target-private block, breakpoint, selected feature
- manual coordinate entry
- gene search
- track toggles
- hover/click detail panel
- click-and-drag horizontal panning when zoomed in
- feature group dropdown with a fixed-height internal scroll list
- click a feature list row to zoom to that interval

### Feature Lists

Candidate feature ranking should be evidence-first:

- gene body
- exon
- CDS
- intron
- promoter
- repeat/SSR/other generic GFF3 feature

Rows should require target-private variant support. Functional annotation can
modify rank or group membership, but should not create a candidate in the
absence of relevant variation.

Feature TSV columns should include:

- rank
- feature type
- feature ID/name
- contig/start/end
- gene, if applicable
- target-private total
- target-private SNP count
- target-private INDEL/complex count
- target-private SV-like count
- all focal variant count
- rank score
- functional category
- representative function
- screening note
- extra track/class metadata

## Internal Architecture

Suggested modules:

```text
src/privy/cli/interactive.py
src/privy/interactive/__init__.py
src/privy/interactive/focus.py
src/privy/interactive/genotypes.py
src/privy/interactive/gff_tracks.py
src/privy/interactive/models.py
src/privy/interactive/render.py
src/privy/interactive/templates.py
```

`src/privy/cli/main.py` would register:

```python
from privy.cli import interactive

app.add_typer(interactive.app, name="interactive")
```

Keep the HTML/JavaScript template self-contained and generated by package code.
Avoid an external web framework for the first implementation. This should be
closer to `privy report --format html` than to a hosted app.

### Reuse and Refactor Targets

The current reusable prototype lives in the Codex skill:

```text
/Users/rothconrad/.codex/skills/privy-focused-interval-report/scripts/build_interval_browser_report.py
```

When moving it into Panex Privus, refactor it into testable functions instead
of copying it as one large script:

- interval parsing
- VCF/focal genotype extraction
- gene model parsing
- generic GFF3 parsing
- variant classification
- target-private pattern classification
- feature scoring
- HTML rendering
- metadata writing

Use existing Panex Privus conventions where possible:

- 0-based half-open intervals internally
- 1-based closed coordinates in user-facing CLI and browser labels
- `privy.io.tsv.read_tsv` / `TsvWriter`
- `privy.io.gff` parsing patterns, expanded as needed for gene models and
  generic tracks
- existing Typer CLI structure and global state

## Future Dashboard Modes

These should live in the same namespace but can be separate implementation
phases after `--focus`.

### `privy interactive --landscape`

Purpose: explore local ancestry/background and genotype-similarity windows.

Potential views:

- sample-by-window private burden heatmap
- missingness heatmap
- nearest background/local background map
- pairwise similarity matrix
- candidate introgression blocks
- linked brushing: click a window to show sample genotypes, nearest background,
  missingness, and relevant block calls
- contig selector and genome-wide mini-map

Inputs:

- `windows.tsv`
- `sample_windows.tsv`
- `similarity.tsv`
- `candidate_introgression_blocks.tsv`
- `background_blocks.tsv`
- `landscape.json`

### `privy interactive --scan`

Purpose: inspect target-private hits and candidate regions from VCF/GFA/BAM/XMFA
scan outputs.

Potential views:

- ranked hit table with strictness filters
- locus/region browser
- strictness class distribution
- score component breakdown
- evidence support bars by sample
- contradiction and missingness summaries
- links from hits to focus-browser intervals

Inputs:

- `hits.tsv`
- `regions.tsv`
- `evidence.tsv`
- `sample_support.tsv`
- `qc.tsv`
- `run.json`

### `privy interactive --pangenome`

Purpose: explore pangenome composition and sample/feature coverage summaries.

Potential views:

- growth curves
- feature composition
- coverage histograms
- sample presence/absence or support heatmaps
- feature table with search/filter
- links from pangenome features to scan/focus views when coordinate-backed

Inputs:

- `pangenome.json`
- `feature_summary.tsv`
- `composition.tsv`
- `coverage_histogram.tsv`
- `growth_curves.tsv`

## Development Milestones

### Milestone 1: Design Spike

- Add this plan to `docs/dev`.
- Decide final user-facing name: `--focus` versus `--region` alias.
- Decide whether repeated focus regions produce one multi-tab HTML or one HTML
  per focus region plus an index. Decision: one HTML per focus region plus an
  index only for multi-region runs.

### Milestone 2: Focus MVP

- Add `privy interactive` CLI.
- Parse repeated `--focus` intervals.
- Accept `--sites-tsv` for precomputed focal genotype tables.
- Render one self-contained focus HTML from a TSV and GFF3.
- Write feature TSV and JSON metadata.
- Add unit tests for parsing, GFF features, variant classes, feature ranking,
  and HTML structure.

### Milestone 3: VCF Extraction Path

- Add direct `--vcf` extraction for focus regions. Done.
- Support `--samples OFFTARGET DERIVED DONOR`. Done.
- Preserve extraction filters in metadata. Done.
- Write `*.sites.tsv` for auditability. Done.

### Milestone 4: Track Expansion

- Add repeated `--track-gff LABEL=PATH`.
- Add functional TSV joins and keyword groups.
- Add optional landscape/scan overlays from existing Privy output directories.

### Milestone 5: Multi-Focus Index

- Write `index.html` summarizing all intervals.
- Add links to each focus HTML and feature TSV.
- Add cross-focus summary counts and input provenance.

### Milestone 6: Landscape/Scan/Pangenome Dashboards

- Build shared interactive HTML components.
- Add mode-specific data loaders and renderers.
- Keep each mode self-contained and reproducible.

## Open Questions

1. Should `--focus` accept only coordinate strings, or also BED/GFF3 interval
   files through `--focus-file`?
2. Should `--samples OFFTARGET DERIVED DONOR` remain the first focus model, or
   should MVP use standard `--targets/--off-targets` with an optional
   `--donor-sample`?
3. Should added GFF3 tracks be embedded in every focus-region HTML or
   centralized into one multi-region data bundle?
4. Should the first implementation support only VCF-derived focal variants, or
   should GFA scan hits be displayable in the focus browser from day one?
5. How much of the original Markdown/HTML scientific report should the command
   generate automatically versus only the interactive browser?
