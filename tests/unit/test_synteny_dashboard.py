"""Unit tests for src/privy/interactive/synteny_dashboard.py."""

from __future__ import annotations

from pathlib import Path

import pytest

from privy.backends.synteny import run_synteny
from privy.interactive.synteny_dashboard import build_synteny_dashboard
from privy.synthetic import presence_absence_pangenome


def _synteny_dir(tmp_path: Path) -> Path:
    gfa = presence_absence_pangenome(seg_len=10).write(tmp_path / "g.gfa")
    out = tmp_path / "syn"
    run_synteny(
        gfa, reference="sample0#0#chr1",
        targets=["sample1"], off_targets=["sample2", "sample3"],
        outdir=out,
    )
    return out


class TestBuildDashboard:
    def test_writes_self_contained_html(self, tmp_path):
        out = build_synteny_dashboard(_synteny_dir(tmp_path))
        assert out.name == "synteny_dashboard.html"
        text = out.read_text(encoding="utf-8")
        # the placeholder JSON literal is replaced (the JS still references the key name)
        assert '{"__privy_placeholder__": true}' not in text
        # self-contained: no external script/style references
        assert "src=" not in text.split("<body")[0]  # no external scripts in head
        assert "<script" in text                       # inlined bundle present
        # injected data present
        assert "sample0" in text
        assert '"blocks"' in text
        assert "target_private" in text

    def test_injected_data_has_private_region(self, tmp_path):
        out = build_synteny_dashboard(_synteny_dir(tmp_path))
        text = out.read_text(encoding="utf-8")
        # presence/absence fixture -> a target-private region -> n_target_private >= 1
        assert '"n_target_private":1' in text or '"n_target_private": 1' in text

    def test_outdir_override(self, tmp_path):
        figs = tmp_path / "dash"
        out = build_synteny_dashboard(_synteny_dir(tmp_path), outdir=figs)
        assert out.parent == figs
        assert out.exists()

    def test_missing_blocks_tsv_raises(self, tmp_path):
        empty = tmp_path / "empty"
        empty.mkdir()
        with pytest.raises(FileNotFoundError, match="synteny_blocks.tsv"):
            build_synteny_dashboard(empty)
