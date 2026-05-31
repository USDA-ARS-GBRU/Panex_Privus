# Validation & benchmarking plan

> The comparative-pangenome layer is validated on synthetic fixtures with known
> answers (in CI). This page is the turnkey plan for **real-data validation and
> benchmarking against established tools** — run on the UGA Sapelo2 cluster
> (see [Sapelo2 live runs](sapelo2-runs.md)). It backs the manuscript's accuracy claims.

## Principle

For each Privy component, compare its output to a community-standard tool on the
*same* real input, quantify agreement, and document where (and why) they differ.
Privy's methods are MVP/heuristic in places (translocation reclassification,
core-segment microhaplotype definition, no context-Jaccard repeat handling), so
the goal is to characterize behavior honestly, not to claim parity everywhere.

## Component → reference tool → metric

| Privy component | Compare against | Input | Agreement metric |
|---|---|---|---|
| `privy project` (coordinate projection) | `odgi position` | same GFA + query positions | exact-match rate of projected `(contig, pos)`; report mismatches |
| `privy synteny` typed blocks | **SyRI** (+ MUMmer/minimap2 aln) and **GENESPACE** | same assemblies/graph | block boundary overlap (reciprocal); INV/TRANS/DUP confusion matrix |
| `privy synteny --paf` chainer | **MCScanX** / jcvi MCscan | same PAF/anchors | collinear-block recall & precision vs MCScanX blocks |
| `privy microhap` loci + copy-number | `vg deconstruct` snarls; alfalfa/blueberry copy-number from the microhap paper | same graph | locus concordance; copy-number agreement on known polyploids |
| `privy popgen` diversity (He, Ne) | **hierfstat** / **adegenet** | same allele matrix (VCF) | per-locus He difference (expect ~0) |
| `privy popgen` F_ST / G_ST / Jost D | **hierfstat** (`wc`), **mmod**, scikit-allel | same cohorts | per-locus + genome-wide difference |
| polyploid F_ST | **StAMPP** | same dosage | difference; confirm frequency-level (no gametic phase) match |
| VanRaden GRM | **rrBLUP** `A.mat`, **AGHmatrix** | same dosage matrix | max abs element difference (expect ~1e-6) |
| private allelic richness | **ADZE** | same cohorts | rarefied-richness difference |
| PCA | `smartpca` (EIGENSOFT) | same dosage | PC correlation (\|r\| ≈ 1, up to sign) |
| DAPC | adegenet `dapc` / DAPCy | same labels | assignment-accuracy agreement |

## Procedure (per dataset)

1. Pick a real crop pangenome (soybean primary; wheat/cotton allopolyploid;
   blueberry autotetraploid) — see [Sapelo2 live runs](sapelo2-runs.md).
2. Run the relevant `privy` command and the reference tool on identical input.
3. Compute the agreement metric; record in a results table (one row per
   dataset × component).
4. Investigate every material discrepancy; classify as (a) expected
   (documented heuristic), (b) bug (fix + add a regression test), or (c) tool
   difference (document).
5. Capture timing + peak memory for the scale section.

## Scale & performance

- Run on whole-chromosome and whole-genome graphs; record wall-time + peak RSS.
- Identify hot paths (likely: GFA parsing, `PathCoordinateModel` build, chaining DP)
  and confirm streaming/index behavior holds; document practical size limits.

## Manuscript outputs

- A benchmark table (the matrix above, filled with real numbers).
- A figure: a real riparian + dotplot + block-density (publication SVG via the
  dashboard's Export SVG or `privy plot`), plus a private-region/diagnostic-marker
  example mirroring the GENESPACE Mo18W/B73 private-inversion use case.
- A methods section citing the reimplemented algorithms (MCScanX/DAGchainer, SyRI
  typing, Nei/Weir-Cockerham/Jost, VanRaden, Kalinowski/ADZE, PHG hVCF, PanSN/rGFA).

## Status

- ✅ Synthetic-fixture correctness (CI, 940+ tests; closed-form checks for pop-gen).
- ✅ Full-pipeline end-to-end integration test (synteny→project→microhap→popgen→dashboard).
- ☐ Real-data runs on Sapelo2 (requires cluster).
- ☐ Tool-vs-tool benchmarks (requires odgi/SyRI/GENESPACE/hierfstat/StAMPP/ADZE installed).
- ☐ Scale/performance on GB-scale graphs.
- ☐ Manuscript figures + benchmark tables.
