---
title: Team tutorial
description: A hands-on, concepts-first tutorial for the Panex Privus comparative-pangenome layer.
---

# Team tutorial — comparative pangenome layer

This teaches the new comparative-pangenome tools by **running them on built-in
synthetic data** (no downloads, no cluster). By the end you'll understand what
each command does, what its outputs mean, and how they all answer one question:

> **What is genetically unique to my target group, where does it sit, and is it a
> usable marker?**

Work through it top to bottom; every command is runnable as written. When you're
comfortable, the [testing guide](testing-guide.md) and
[command guide](comparative-pangenome.md) will make full sense.

---

## The mental model (read this once)

- A **pangenome graph** (GFA) packs many genomes into one structure. Each genome
  is a **path** that walks through shared **segments** (nodes). Two genomes share
  a segment exactly when they share that sequence.
- Genomes are named **PanSN** style: `sample#haplotype#contig` (e.g.
  `Benning#0#Gm01`). A polyploid sample has several haplotype paths.
- You split samples into a **target** cohort (your group of interest) and an
  **off-target** cohort (everything else). Privy's core job is finding signal
  **present in the target and absent from the off-target** — "target-private".
- Two coordinate worlds: a base is either `(segment, offset)` (**unstable** — changes
  if the graph is rebuilt) or `(reference, position)` (**stable** — what you cite).
  Privy converts between them, which is how it **projects** anything to any reference.

The new layer adds four ideas on top of that core:

| Idea | Command | One-liner |
|---|---|---|
| **Synteny** | `privy synteny` | where genomes are collinear, and where they're rearranged (inversion/translocation/duplication) — and which regions are target-private |
| **Projection** | `privy project` | take a region in one genome → find it in every other genome |
| **Microhaplotypes** | `privy microhap` | local multi-allelic markers; flag alleles private to the target |
| **Breeder pop-gen** | `privy popgen` | diversity, target-vs-off-target differentiation, diagnostic markers, genomic-prediction matrices |

…plus `privy plot` (static figures) and `privy dashboard` (interactive HTML).

---

## Setup

```bash
conda activate privy            # an env with privy installed (see Installation)
mkdir privy-tutorial && cd privy-tutorial
```

We'll generate small, realistic-looking graphs with Privy's synthetic generator —
the same fixtures the tests use. Make one with a **target-private deletion**
(targets keep an `s2,s3` block that off-targets deleted):

```bash
python - <<'PY'
from privy.synthetic import presence_absence_pangenome, inversion_pangenome
presence_absence_pangenome(seg_len=400).write("pa.gfa")     # target-private signal
inversion_pangenome(seg_len=400).write("inv.gfa")           # a nice rearrangement for plots
print("wrote pa.gfa and inv.gfa")
PY
```

`pa.gfa` has four genomes: `sample0`/`sample1` (targets, carry the block) and
`sample2`/`sample3` (off-targets, deleted it). `sample0#0#chr1` is our reference.

---

## Part 1 — Synteny: where, and what's rearranged or private

```bash
privy synteny --gfa pa.gfa --reference sample0#0#chr1 \
  --targets sample1 --off-targets sample2,sample3 --outdir out/synteny/
```

Open `out/synteny/synteny_blocks.tsv` and `synteny_regions.tsv`:

- **`synteny_blocks.tsv`** — one row per collinear/rearranged block between a query
  genome and the reference. The `block_type` column is `collinear`, `inversion`,
  `translocation`, or `duplication`.
- **`synteny_regions.tsv`** — blocks merged on the reference, with a
  **`target_private`** column. The region covering the deleted `s2,s3` block is
  `True`: present in the targets, absent from the off-targets. **That's a
  target-private structural region.**
- `synteny.json` summarizes block-type counts and how many regions are private.

Try the inversion graph to see typed rearrangements:

```bash
privy synteny --gfa inv.gfa --reference sample0#0#chr1 --outdir out/synteny_inv/
grep inversion out/synteny_inv/synteny_blocks.tsv     # the planted inversion, typed
```

PAF mode (anchors from `odgi untangle` / minimap2 / wfmash) instead of a graph:
`privy synteny --paf anchors.paf --outdir ...`.

---

## Part 2 — Projection: find a region in every genome

Define a region by its graph segments and project it onto all genomes:

```bash
privy project --gfa pa.gfa --node-set s2,s3 --outdir out/project/
column -t out/project/projection.tsv
```

`projection.tsv` lists, per genome, where `s2,s3` lands (`present=True` with
coordinates) or that it's absent. You'll see it **present in `sample0`/`sample1`
and absent in `sample2`/`sample3`** — the same private signal, now as coordinates
you can hand to a browser. You can also project a coordinate range:
`--region sample0#0#chr1:800-1600`.

---

## Part 3 — Microhaplotypes: multi-allelic markers + private alleles

```bash
privy microhap --gfa pa.gfa --reference sample0#0#chr1 \
  --targets sample0,sample1 --off-targets sample2,sample3 --outdir out/microhap/
column -t out/microhap/microhaplotypes.tsv
```

- A **microhaplotype** is a short locus with several linked variants → a
  multi-allelic marker (richer than a yes/no SNP). Each genome's local sequence is
  an **allele**, content-hashed (MD5) so identical sequence = identical allele id.
- `microhaplotypes.tsv` shows `n_alleles`, `aaf` (combined alternative allele
  frequency), and **`target_private`** (an allele present only in the target cohort).
- `allele_matrix.tsv` is loci × genomes with integer allele indices — the input to
  allele-aware analyses.

---

## Part 4 — Breeder population genetics: is it a usable marker?

```bash
privy popgen --gfa pa.gfa --reference sample0#0#chr1 \
  --targets sample0,sample1 --off-targets sample2,sample3 --outdir out/popgen/
column -t out/popgen/popgen_loci.tsv
cat out/popgen/popgen.json
```

- `popgen_loci.tsv` — per locus: diversity (`gene_diversity`, `effective_alleles`),
  `fis`, and **differentiation** between cohorts (`gst`, `jost_d`, `is_diagnostic`).
  A locus with **`gst` ≈ 1 and `is_diagnostic = True`** perfectly separates target
  from off-target — a **candidate selection marker**.
- `popgen.json` — `genome_wide_fst`, number of diagnostic markers, and per-cohort
  private-allele counts/richness.
- `grm.tsv` (VanRaden genomic relationship matrix) and `dosage_matrix.tsv` are
  **genomic-prediction inputs** — hand them to rrBLUP / BGLR / sommer. Privy
  produces the inputs; it does not fit the prediction model.
- `pca.tsv` (when ≥2 loci) gives population-structure coordinates.

Polyploids are handled automatically: a sample's several haplotype paths are
grouped into one genotype, so allele **dosage** runs 0..ploidy.

---

## Part 5 — Visualize

Static, publication figures:

```bash
privy plot --plot-set synteny --input-dir out/synteny_inv/
# -> riparian.png, dotplot.png, block_density.png
```

Interactive, self-contained dashboard (open the `.html` in any browser — no
server, works offline, e-mailable):

```bash
privy dashboard --synteny out/synteny_inv/                 # synteny only
privy dashboard --synteny out/synteny/ --microhap out/microhap/   # + allele panel
open out/synteny_inv/synteny_dashboard.html
```

Hover a braid to highlight it in the dotplot, toggle "private regions only", read
the per-block detail panel, browse the microhaplotype/allele panel, and click
**Export SVG** for a figure.

---

## The through-line (why these are one tool, not five)

For the `pa.gfa` example, the *same* target-private signal shows up at every layer:

1. **synteny** → a `target_private = True` region (structural: the targets have a
   block the off-targets lack),
2. **project** → that block is present in the targets, absent in the off-targets,
   with coordinates,
3. **microhap** → a `target_private = True` multi-allelic locus there,
4. **popgen** → that locus is a fully diagnostic marker (`gst = 1`), exported in
   the dosage matrix / GRM for prediction,
5. **dashboard** → all of it, highlighted and shareable.

That consistency is exactly what the
[end-to-end test](testing-guide.md) checks automatically.

---

## Mini-glossary

- **GFA** — pangenome graph file format (segments, links, paths/walks).
- **PanSN** — `sample#haplotype#contig` naming convention.
- **target / off-target** — your group of interest vs everything else.
- **target-private** — present in the target cohort, absent from the off-target.
- **synteny block** — a run of shared, ordered segments between two genomes.
- **collinear / inversion / translocation / duplication** — block types (how a
  region is arranged relative to the reference).
- **projection** — mapping a region/coordinate from one genome onto another via the graph.
- **microhaplotype** — a short multi-allelic locus (several linked variants).
- **allele (MD5)** — a local haplotype sequence identified by its content hash.
- **AAF** — combined alternative allele frequency (multi-allelic analogue of MAF).
- **F_ST / Jost's D** — how strongly a locus differentiates two cohorts (1 = fully diagnostic).
- **F_IS** — heterozygote deficit/excess within genotypes.
- **dosage** — copies of an allele in a (polyploid) genotype, 0..ploidy.
- **GRM** — genomic relationship matrix; a genomic-prediction input.

---

## Where to go next

- [Comparative pangenome command guide](comparative-pangenome.md) — every flag and output.
- [Testing guide](testing-guide.md) — run the suite and read the tests.
- [Sapelo2 live runs](sapelo2-runs.md) — when you're ready to run on real crop data.
- [Validation & benchmarking plan](benchmarking.md) — how we'll prove accuracy.

> Reminder: this layer is **synthetic-validated and under active development** —
> great for learning the tools and reviewing behaviour, but treat results on real
> data as candidate signal until the benchmarking phase is done.
