---
title: Testing guide
description: How to run and read the Panex Privus test suite (comparative-pangenome layer).
---

# Testing guide

This walks you through running and reading the test suite for the
comparative-pangenome layer, so you can review every component as we validate it.
Everything here runs **locally on synthetic data** — no cluster, no downloads.
Real-data validation is separate (see [Sapelo2 live runs](sapelo2-runs.md) and
the [validation plan](benchmarking.md)).

## 1. Set up a dev environment

```bash
mamba create -n privy-dev -c conda-forge -c bioconda python=3.11 pysam samtools bcftools htslib
conda activate privy-dev
cd Panex_Privus
git checkout feature/comparative-pangenome-p0
pip install -e ".[dev]"        # core + test/lint/type tooling
pip install -e ".[dev,full]"   # add optional extras (scipy, scikit-learn, gfapy, …) to run every test
```

> Without the `full` extra, one DAPC test is **skipped** (scikit-learn absent) and
> everything else still passes — that is expected and confirms graceful degradation.

## 2. Run the whole suite

```bash
pytest                         # ~945 tests, a few seconds
pytest -q                      # quiet
pytest --cov=privy --cov-report=term-missing   # with coverage
```

Expected: **944 passed, 1 skipped** (or 945 passed with the `full` extra).

## 3. Run by area (review one piece at a time)

| Area | Command |
|---|---|
| Anchors / I/O | `pytest tests/unit/test_paf_io.py tests/unit/test_hvcf_io.py tests/unit/test_madc_deconstruct_io.py` |
| Graph coordinates & projection | `pytest tests/unit/test_synteny_coordinates.py tests/unit/test_synteny_projection.py` |
| Synteny model & typed blocks | `pytest tests/unit/test_synteny_model.py tests/unit/test_synteny_graph_blocks.py tests/unit/test_synteny_build.py` |
| Synteny chainer & gene mode | `pytest tests/unit/test_synteny_chain.py tests/unit/test_synteny_genes.py` |
| Microhaplotypes | `pytest tests/unit/test_microhap.py` |
| Polyploid dosage | `pytest tests/unit/test_polyploid.py` |
| Population genetics | `pytest tests/unit/test_popgen.py tests/unit/test_popgen_private.py tests/unit/test_popgen_relationship.py tests/unit/test_popgen_structure.py` |
| Chromosome structure | `pytest tests/unit/test_structure.py` |
| Interval / density utilities | `pytest tests/unit/test_intervals_classes.py` |
| Static figures | `pytest tests/unit/test_plot_synteny.py` |
| Synthetic data generator | `pytest tests/unit/test_synthetic.py` |
| Dashboard builder | `pytest tests/unit/test_synteny_dashboard.py` |
| CLI (integration) | `pytest tests/integration/test_{project,synteny,microhap,popgen,plot_synteny,dashboard}_cli.py` |
| **Full pipeline** | `pytest tests/integration/test_end_to_end.py -v` |

The **end-to-end test** is the one to read first — it runs every command on one
synthetic graph and asserts the target-private signal threads consistently from
synteny → projection → microhaplotypes → pop-gen → dashboard.

## 4. Lint and type checks (same gates used in development)

```bash
ruff check src tests           # style/lint (line length 100, import order, etc.)
mypy src/privy                 # strict static typing
```

Both should report clean.

## 5. Eyeball the interactive dashboard

The dashboard is tested for data injection and self-containment, but not for
browser rendering — so render one yourself:

```bash
python - <<'PY'
from pathlib import Path
from privy.synthetic import inversion_pangenome
from privy.backends.synteny import run_synteny
from privy.interactive.synteny_dashboard import build_synteny_dashboard

d = Path("demo"); d.mkdir(exist_ok=True)
gfa = inversion_pangenome(seg_len=400).write(d / "demo.gfa")
run_synteny(gfa, reference="sample0#0#chr1",
            targets=["sample3"], off_targets=["sample1", "sample2"], outdir=d / "syn")
print(build_synteny_dashboard(d / "syn"))
PY
open demo/syn/synteny_dashboard.html    # or double-click; works offline, no server
```

You should see stacked genome tracks with a coloured inversion braid, a linked
dotplot, a block-type legend, the target-private toggle, and an "Export SVG" button.

## 6. What "passing" does and does not mean

- ✅ It means the algorithms behave correctly on **synthetic graphs with known
  answers** (e.g. a planted inversion is typed as an inversion; a deleted block is
  flagged target-private; a fully diagnostic locus gives F_ST = 1).
- ❌ It does **not** yet mean validated on real crop pangenomes, nor benchmarked
  against established tools, nor profiled at genome scale. That is the next phase —
  see the [validation & benchmarking plan](benchmarking.md).

## 7. New to the tools first?

If you are reviewing the *behaviour* rather than the code, start with the
[team tutorial](team-guide.md): it teaches the concepts and walks each command
hands-on with the synthetic data, which makes the tests far easier to read.
