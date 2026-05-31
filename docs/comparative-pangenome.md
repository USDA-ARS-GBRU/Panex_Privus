# Comparative pangenome workflows

> Status: under active development (the comparative-pangenome layer is new in the
> `feature/comparative-pangenome-p0` line). Commands and outputs may change.

Panex Privus consumes an existing pangenome graph (GFA from minigraph-cactus or
PGGB) and adds a comparative layer on top of its target-private core:

- **`privy project`** — project a region or graph node-set to *any* reference in the graph
- **`privy synteny`** — graph- or PAF-derived synteny blocks with typed rearrangements + target-private regions
- **`privy microhap`** — multi-allelic microhaplotype loci + target-private alleles
- **`privy popgen`** — allelic diversity, cohort differentiation, diagnostic markers, GP-ready matrices
- **`privy plot --plot-set synteny`** — static riparian + dotplot figures
- **`privy dashboard --synteny`** — a self-contained interactive HTML dashboard

All commands are 0-based half-open internally and identify genomes by PanSN path
names (`sample#haplotype#contig`). Cohorts (`--targets` / `--off-targets`) accept
either PanSN sample names or full path ids, comma-separated.

## Try it locally with synthetic data

No downloads needed — `privy.synthetic` generates tiny crop-like graphs:

```python
from privy.synthetic import inversion_pangenome
inversion_pangenome(seg_len=400).write("demo.gfa")
```

```bash
# graph-native synteny (typed blocks + private regions)
privy synteny --gfa demo.gfa --reference sample0#0#chr1 \
  --targets sample3 --off-targets sample1,sample2 --outdir out/

# static figures + interactive dashboard from the synteny output
privy plot --plot-set synteny --input-dir out/
privy dashboard --synteny out/      # → out/synteny_dashboard.html (open in any browser)
```

## A typical workflow

1. **Build/obtain a graph.** Privy consumes GFA; build with minigraph-cactus or
   PGGB (external), or convert an existing graph with `vg convert -f graph.gbz > graph.gfa`.
2. **Synteny.** `privy synteny --gfa graph.gfa --reference <PanSN path> --targets … --off-targets …`
   → `synteny_blocks.tsv`, `synteny_regions.tsv` (with `target_private`), `synteny.json`.
3. **Project a region of interest** to every genome:
   `privy project --gfa graph.gfa --region <path>:<start>-<end>` → `projection.tsv`.
   Or define it in node space: `--node-set s12,s13,s14`.
4. **Microhaplotypes + private alleles.**
   `privy microhap --gfa graph.gfa --reference <path> --targets … --off-targets …`
   → `microhaplotypes.tsv`, `allele_matrix.tsv`.
5. **Breeder population genetics.**
   `privy popgen --gfa graph.gfa --reference <path> --targets … --off-targets …`
   → `popgen_loci.tsv` (per-locus diversity + G_ST/Jost's D + `is_diagnostic`),
   `dosage_matrix.tsv`, `grm.tsv` (VanRaden, GP-ready), `popgen.json`
   (genome-wide F_ST, # diagnostic markers).
6. **Visualize.** `privy plot --plot-set synteny` for publication figures;
   `privy dashboard --synteny` for an interactive, shareable single-file HTML.

## Interpreting the breeder output

- A `popgen_loci.tsv` locus with **`gst` ≈ 1.0 / `is_diagnostic = True`** is a fully
  diagnostic marker: an allele present in the target cohort and absent from the
  off-target cohort — a candidate selection marker.
- `grm.tsv` (VanRaden genomic relationship matrix) and `dosage_matrix.tsv`
  (samples × loci, 0..ploidy) are inputs for genomic prediction — hand them to
  rrBLUP / BGLR / sommer. Privy produces the inputs; it does not fit GP models.
- Polyploid genotypes are assembled from a sample's PanSN haplotype paths, so
  dosage runs 0..ploidy automatically.

## Interoperability (optional)

Privy ingests standard outputs rather than running heavy tools itself:

- **PAF** from `odgi untangle` / minimap2 / wfmash → `privy synteny --paf aln.paf`
- **GBZ → GFA** via `vg convert` before any command
- (planned) `vg deconstruct` VCF, PHG hVCF, DArTag MADC

See also: [Sapelo2 live-run instructions](sapelo2-runs.md) for running on real crop
pangenomes at scale.
