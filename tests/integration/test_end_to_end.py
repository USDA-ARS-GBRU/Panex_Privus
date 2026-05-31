"""End-to-end pipeline integration test (P8 reproducible example).

Runs the full comparative-pangenome pipeline on a single synthetic graph that
carries a known target-private signal (targets keep an s2,s3 block that
off-targets deleted), and asserts the signal threads consistently across every
command: synteny -> project -> microhap -> popgen -> dashboard.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

from typer.testing import CliRunner

from privy.cli.main import app
from privy.synthetic import presence_absence_pangenome

runner = CliRunner()

REF = "sample0#0#chr1"
TARGETS = "sample1"
OFFTARGETS = "sample2,sample3"


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_full_pipeline_threads_private_signal(tmp_path):
    gfa = presence_absence_pangenome(seg_len=10).write(tmp_path / "pan.gfa")

    # 1. synteny -> target-private structural region
    syn = tmp_path / "synteny"
    r = runner.invoke(app, [
        "synteny", "--gfa", str(gfa), "--reference", REF,
        "--targets", TARGETS, "--off-targets", OFFTARGETS, "--outdir", str(syn),
    ])
    assert r.exit_code == 0, r.output
    regions = _read_tsv(syn / "synteny_regions.tsv")
    assert any(row["target_private"] == "True" for row in regions)

    # 2. project the private region's segments to all genomes
    proj = tmp_path / "proj"
    r = runner.invoke(app, [
        "project", "--gfa", str(gfa), "--node-set", "s2,s3", "--outdir", str(proj),
    ])
    assert r.exit_code == 0, r.output
    pj = _read_tsv(proj / "projection.tsv")
    present = {row["target_path"] for row in pj if row["present"] == "True"}
    # the s2,s3 block is present in the targets and absent from the off-targets
    assert "sample0#0#chr1" in present and "sample1#0#chr1" in present
    assert "sample2#0#chr1" not in present and "sample3#0#chr1" not in present

    # 3. microhap -> a target-private multi-allelic locus
    mh = tmp_path / "microhap"
    r = runner.invoke(app, [
        "microhap", "--gfa", str(gfa), "--reference", REF,
        "--targets", "sample0,sample1", "--off-targets", OFFTARGETS, "--outdir", str(mh),
    ])
    assert r.exit_code == 0, r.output
    loci = _read_tsv(mh / "microhaplotypes.tsv")
    assert any(row["target_private"] == "True" for row in loci)

    # 4. popgen -> the locus is a fully diagnostic marker; GP matrices written
    pg = tmp_path / "popgen"
    r = runner.invoke(app, [
        "popgen", "--gfa", str(gfa), "--reference", REF,
        "--targets", "sample0,sample1", "--off-targets", OFFTARGETS, "--outdir", str(pg),
    ])
    assert r.exit_code == 0, r.output
    meta = json.loads((pg / "popgen.json").read_text())
    assert meta["genome_wide_fst"] == 1.0
    assert meta["n_diagnostic_loci"] >= 1
    assert (pg / "grm.tsv").exists() and (pg / "dosage_matrix.tsv").exists()

    # 5. dashboard (synteny + microhap) -> self-contained HTML carrying the signal
    r = runner.invoke(app, [
        "dashboard", "--synteny", str(syn), "--microhap", str(mh),
    ])
    assert r.exit_code == 0, r.output
    html = (syn / "synteny_dashboard.html").read_text(encoding="utf-8")
    assert '{"__privy_placeholder__": true}' not in html   # data injected
    assert '"microhaplotypes"' in html                     # allele layer present
    assert "src=" not in html.split("<body")[0]            # self-contained


def test_plot_set_synteny_renders_all_three_figures(tmp_path):
    gfa = presence_absence_pangenome(seg_len=10).write(tmp_path / "pan.gfa")
    syn = tmp_path / "synteny"
    runner.invoke(app, [
        "synteny", "--gfa", str(gfa), "--reference", REF,
        "--targets", TARGETS, "--off-targets", OFFTARGETS, "--outdir", str(syn),
    ])
    r = runner.invoke(app, ["plot", "--plot-set", "synteny", "--input-dir", str(syn)])
    assert r.exit_code == 0, r.output
    for fig in ("riparian.png", "dotplot.png", "block_density.png"):
        assert (syn / fig).exists()
