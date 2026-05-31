# Real-data runs on the UGA Sapelo2 cluster

> Development uses tiny synthetic fixtures (`privy.synthetic`) — no downloads.
> Real-world, GB-scale validation happens on the **UGA Sapelo2** cluster. This
> page gives per-crop recipes: where to get data, how to make a Privy-ready GFA,
> the exact `privy` commands, and SLURM resource hints.
>
> **Verify licenses, accessions, and file layouts at run time** — public crop
> resources move and re-version. Sizes are approximate.

## General setup

```bash
# one-time: a dedicated environment (see installation docs)
mamba create -n privy -c conda-forge -c bioconda python=3.11 pysam samtools bcftools htslib
conda activate privy
pip install -U .            # from a Panex Privus checkout
# optional graph/aligner tools (Tier-2, only for preparing inputs):
mamba install -c bioconda vg odgi minimap2 gfatools
```

Privy consumes a **GFA** graph. If a resource ships GBZ or HAL, convert first:

```bash
vg convert -f graph.gbz > graph.gfa            # GBZ -> GFA
# or extract one chromosome to keep memory modest:
vg chunk -x graph.gbz -p "REF#0#chr1" -O gfa > chr1.gfa
```

### Generic SLURM wrapper (Sapelo2)

```bash
#!/bin/bash
#SBATCH --job-name=privy
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G            # raise for whole-genome graphs (see per-crop notes)
#SBATCH --time=24:00:00
#SBATCH --output=privy_%j.out
conda activate privy
privy synteny --gfa "$GFA" --reference "$REF" \
  --targets "$TARGETS" --off-targets "$OFFTARGETS" --outdir results/synteny/
privy microhap --gfa "$GFA" --reference "$REF" \
  --targets "$TARGETS" --off-targets "$OFFTARGETS" --outdir results/microhap/
privy popgen --gfa "$GFA" --reference "$REF" \
  --targets "$TARGETS" --off-targets "$OFFTARGETS" --outdir results/popgen/
privy dashboard --synteny results/synteny/
```

Tip: start on a single chromosome subgraph before a whole-genome run.

## Soybean — *Glycine max* (primary system; paleopolyploid)

- **Data:** SoyBase pangenome collections (https://www.soybase.org/collections/) +
  resequenced-accession VCFs; Wm82 reference (a2/a4) on SoyBase / Phytozome.
- **Make a GFA:** no ready public GFA — build a small minigraph-cactus or PGGB graph
  from a handful of accessions vs Wm82 (see the GigaScience 2025 soybean graph recipe).
- **Run:** `REF=Wm82#0#Gm01` (PanSN-name your assemblies accordingly), targets/off-targets
  = your high-/low-trait accessions. Mem ~32–64 G per chromosome subgraph.

## Bread wheat — *Triticum aestivum* (allohexaploid; top allopolyploid test)

- **Data:** Wheat 10+/Panache graph, Zenodo 10.5281/zenodo.6085239 — download only
  `15-wheat10+.gfa.gz` (~4.4 GB; **CC-BY 4.0**); skip the ~127 GB Giraffe indexes.
- **Prep:** `gunzip 15-wheat10+.gfa.gz`; `vg chunk` one chromosome for a tractable first run.
- **Run:** set `REF` to a cultivar's PanSN path; mem 64–128 G per chromosome; expect
  homeolog structure (A/B/D) — interpret with subgenome awareness.

## Cotton — *Gossypium hirsutum* (allotetraploid AD; lighter allopolyploid)

- **Data:** Jin et al. 2023 SV-pangenome (CottonGen, https://www.cottongen.org/) —
  11 assemblies + 182k-SV VCF. Build a small graph from a few assemblies.
- **Run:** two subgenomes (A/D); good case for the polyploid dosage + GP-matrix outputs.

## Blueberry — *Vaccinium* (autotetraploid; published graph)

- **Data:** GDV pangenome (https://www.vaccinium.org/pangenome_graphs); code
  github.com/Aeyocca/VaccPan; SRA PRJNA687008. Graph is HAL/browser-hosted — export
  GFA with `hal2vg` then `vg convert -f`.
- **Run:** true autotetraploid — exercises 0..4 dosage and observed heterozygosity.

## Maize — NAM founders, B73 (paleopolyploid; comparative classic)

- **Data:** Gramene pan-maize / Hufford 2021 (26 NAM genomes; Zenodo 4781590);
  the GENESPACE Mo18W/B73 example (github.com/jtlovell/GENESPACE_data). FASTA/GFF;
  build a graph from NAM assemblies. Strong case for the riparian/dotplot figures.

## Amplicon / microhaplotype panels (DArTag)

- Breeding Insight open DArTag panels + HapApp (github.com/Breeding-Insight/HapApp_utils,
  MIT). A lighter, allele-level on-ramp — a future `io/madc` importer will read MADC
  reports directly; for now interoperate via VCF.

## After the run

Pull `results/**/synteny_dashboard.html` back to your laptop (it is fully
self-contained — open in any browser, no server). The `popgen_loci.tsv`
`is_diagnostic` markers and `grm.tsv` are the breeder-actionable outputs.
