"""Integration tests for privy plot — run_plot() and CLI."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

from privy.cli.main import app
from privy.core.config import default_config
from privy.plot.loci import run_plot

GFA_PATH = Path(__file__).parent.parent / "data" / "small_cohort.gfa"

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

_HITS_HEADER = (
    "locus_id\tcontig\tstart\tend\tvariant_type\tallele_key\t"
    "target_support_n\ttarget_total_n\tofftarget_support_n\t"
    "offtarget_total_n\ttarget_missing_n\tofftarget_missing_n\t"
    "strictness_class\tdiscovery_score\tsupport_score\tpenalty_score\tfinal_score\n"
)

_EVIDENCE_HEADER = (
    "locus_id\tsource_type\tsample_id\tevidence_class\t"
    "metric_name\tmetric_value\tdetails\n"
)

_COMPARE_HEADER = (
    "compare_id\tlocus_id_a\tlocus_id_b\tsource_a\tsource_b\tcontig\t"
    "start_a\tend_a\tstart_b\tend_b\tmatch_class\tcoordinate_overlap\t"
    "state_compatibility\tstrictness_a\tstrictness_b\t"
    "support_summary\tcontradiction_summary\tcomparison_score\n"
)


def _write_hits(tmp_path: Path, n: int = 5) -> Path:
    path = tmp_path / "hits.tsv"
    strictness = ["strict_complete", "strict_target_missing",
                  "relaxed_threshold", "contradicted", "strict_both_missing"]
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(_HITS_HEADER)
        for i in range(n):
            fh.write(
                f"PPX{i:06d}\tchr1\t{i*1000}\t{i*1000+100}\tsnp\t"
                f"chr1:{i*1000}:A:T\t3\t3\t0\t2\t0\t0\t"
                f"{strictness[i % len(strictness)]}\t1.0\t0.0\t0.0\t"
                f"{round(1.0 - i * 0.1, 2)}\n"
            )
    return path


def _write_evidence(tmp_path: Path) -> Path:
    path = tmp_path / "evidence.tsv"
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(_EVIDENCE_HEADER)
        for src in ("vcf", "bam"):
            for ec in ("support", "absence"):
                fh.write(f"PPX000001\t{src}\tS1\t{ec}\tdepth\t12.0\t\n")
    return path


def _write_compare(tmp_path: Path) -> Path:
    path = tmp_path / "compare.tsv"
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(_COMPARE_HEADER)
        for i, mc in enumerate(("supported", "source_specific", "contradicted")):
            fh.write(
                f"CMP{i:06d}\tPPX{i:06d}\t"
                f"{'GPX' + str(i).zfill(6) if mc != 'source_specific' else 'NA'}\t"
                f"vcf\tgfa\tchr1\t{i*1000}\t{i*1000+100}\t"
                f"{'NA' if mc == 'source_specific' else str(i*1000)}\t"
                f"{'NA' if mc == 'source_specific' else str(i*1000+100)}\t"
                f"{mc}\t{'0.9' if mc != 'source_specific' else '0.0'}\t"
                f"True\tstrict_complete\t"
                f"{'strict_complete' if mc != 'source_specific' else 'NA'}\t"
                f"NA\tNA\t0.9\n"
            )
    return path


@pytest.fixture()
def hits(tmp_path: Path) -> Path:
    return _write_hits(tmp_path)


@pytest.fixture()
def evidence(tmp_path: Path) -> Path:
    return _write_evidence(tmp_path)


@pytest.fixture()
def compare(tmp_path: Path) -> Path:
    return _write_compare(tmp_path)


@pytest.fixture()
def runner() -> CliRunner:
    return CliRunner()


# ---------------------------------------------------------------------------
# TestRunPlotAllMode
# ---------------------------------------------------------------------------

class TestRunPlotAllMode:
    def test_all_generates_three_base_plots(
        self, tmp_path: Path, hits: Path
    ) -> None:
        outdir = tmp_path / "out"
        generated = run_plot(
            hits=hits, regions=None, evidence=None, vcf=None, bam=None,
            bam_manifest=None, gfa=None, xmfa=None, compare=None,
            cfg=default_config(), locus_id=None, region_id=None,
            top_n=None, contig=None, region=None, plot_type="all",
            width=8.0, height=4.0, dpi=72, output_format="png",
            show_labels=True, outdir=outdir,
        )
        names = {p.name for p in generated}
        assert "locus_panel.png" in names
        assert "strictness_bar.png" in names
        assert "score_distribution.png" in names

    def test_all_with_evidence_adds_support_bar(
        self, tmp_path: Path, hits: Path, evidence: Path
    ) -> None:
        outdir = tmp_path / "out"
        generated = run_plot(
            hits=hits, regions=None, evidence=evidence, vcf=None, bam=None,
            bam_manifest=None, gfa=None, xmfa=None, compare=None,
            cfg=default_config(), locus_id=None, region_id=None,
            top_n=None, contig=None, region=None, plot_type="all",
            width=8.0, height=4.0, dpi=72, output_format="png",
            show_labels=True, outdir=outdir,
        )
        assert any(p.name == "support_bar.png" for p in generated)

    def test_all_with_compare_adds_compare_summary(
        self, tmp_path: Path, hits: Path, compare: Path
    ) -> None:
        outdir = tmp_path / "out"
        generated = run_plot(
            hits=hits, regions=None, evidence=None, vcf=None, bam=None,
            bam_manifest=None, gfa=None, xmfa=None, compare=compare,
            cfg=default_config(), locus_id=None, region_id=None,
            top_n=None, contig=None, region=None, plot_type="all",
            width=8.0, height=4.0, dpi=72, output_format="png",
            show_labels=True, outdir=outdir,
        )
        assert any(p.name == "compare_summary.png" for p in generated)


# ---------------------------------------------------------------------------
# TestRunPlotIndividualTypes
# ---------------------------------------------------------------------------

class TestRunPlotIndividualTypes:
    def _run(self, hits: Path, tmp_path: Path, plot_type: str, **extra: object) -> list[Path]:
        outdir = tmp_path / "out"
        return run_plot(
            hits=hits, regions=None, evidence=extra.get("evidence"),  # type: ignore[arg-type]
            vcf=None, bam=None, bam_manifest=None, gfa=None, xmfa=None,
            compare=extra.get("compare"),  # type: ignore[arg-type]
            cfg=default_config(), locus_id=None, region_id=None,
            top_n=None, contig=None, region=None, plot_type=plot_type,
            width=8.0, height=4.0, dpi=72, output_format="png",
            show_labels=True, outdir=outdir,
        )

    def test_locus_panel_only(self, tmp_path: Path, hits: Path) -> None:
        generated = self._run(hits, tmp_path, "locus_panel")
        assert len(generated) == 1
        assert generated[0].name == "locus_panel.png"

    def test_strictness_bar_only(self, tmp_path: Path, hits: Path) -> None:
        generated = self._run(hits, tmp_path, "strictness_bar")
        assert len(generated) == 1
        assert generated[0].name == "strictness_bar.png"

    def test_score_distribution_only(self, tmp_path: Path, hits: Path) -> None:
        generated = self._run(hits, tmp_path, "score_distribution")
        assert len(generated) == 1
        assert generated[0].name == "score_distribution.png"

    def test_support_bar_without_evidence_raises(
        self, tmp_path: Path, hits: Path
    ) -> None:
        with pytest.raises(ValueError, match="--evidence"):
            self._run(hits, tmp_path, "support_bar")

    def test_compare_summary_without_compare_raises(
        self, tmp_path: Path, hits: Path
    ) -> None:
        with pytest.raises(ValueError, match="--compare"):
            self._run(hits, tmp_path, "compare_summary")

    def test_support_bar_with_evidence(
        self, tmp_path: Path, hits: Path, evidence: Path
    ) -> None:
        generated = self._run(hits, tmp_path, "support_bar", evidence=evidence)
        assert len(generated) == 1
        assert generated[0].name == "support_bar.png"

    def test_compare_summary_with_compare(
        self, tmp_path: Path, hits: Path, compare: Path
    ) -> None:
        generated = self._run(hits, tmp_path, "compare_summary", compare=compare)
        assert len(generated) == 1
        assert generated[0].name == "compare_summary.png"


# ---------------------------------------------------------------------------
# TestRunPlotOutputFormats
# ---------------------------------------------------------------------------

class TestRunPlotOutputFormats:
    def test_svg_output(self, tmp_path: Path, hits: Path) -> None:
        outdir = tmp_path / "out"
        generated = run_plot(
            hits=hits, regions=None, evidence=None, vcf=None, bam=None,
            bam_manifest=None, gfa=None, xmfa=None, compare=None,
            cfg=default_config(), locus_id=None, region_id=None,
            top_n=None, contig=None, region=None, plot_type="locus_panel",
            width=8.0, height=4.0, dpi=72, output_format="svg",
            show_labels=True, outdir=outdir,
        )
        assert generated[0].suffix == ".svg"
        assert generated[0].exists()

    def test_files_are_non_empty(self, tmp_path: Path, hits: Path) -> None:
        outdir = tmp_path / "out"
        generated = run_plot(
            hits=hits, regions=None, evidence=None, vcf=None, bam=None,
            bam_manifest=None, gfa=None, xmfa=None, compare=None,
            cfg=default_config(), locus_id=None, region_id=None,
            top_n=None, contig=None, region=None, plot_type="all",
            width=8.0, height=4.0, dpi=72, output_format="png",
            show_labels=True, outdir=outdir,
        )
        for path in generated:
            assert path.stat().st_size > 200, f"{path.name} is suspiciously small"


# ---------------------------------------------------------------------------
# TestPlotCli
# ---------------------------------------------------------------------------

class TestPlotCli:
    def test_successful_run(
        self, tmp_path: Path, hits: Path, runner: CliRunner
    ) -> None:
        outdir = tmp_path / "out"
        result = runner.invoke(app, [
            "plot",
            "--plot-set", "scan",
            "--hits", str(hits),
            "--plot-type", "strictness_bar",
            "--dpi", "72",
            "--outdir", str(outdir),
        ])
        assert result.exit_code == 0, result.output
        assert (outdir / "strictness_bar.png").exists()

    def test_missing_hits_exits_nonzero(
        self, tmp_path: Path, runner: CliRunner
    ) -> None:
        result = runner.invoke(app, [
            "plot",
            "--hits", str(tmp_path / "nonexistent.tsv"),
            "--outdir", str(tmp_path / "out"),
        ])
        assert result.exit_code != 0

    def test_all_mode_generates_multiple_files(
        self, tmp_path: Path, hits: Path, runner: CliRunner
    ) -> None:
        outdir = tmp_path / "out"
        result = runner.invoke(app, [
            "plot",
            "--hits", str(hits),
            "--plot-type", "all",
            "--dpi", "72",
            "--outdir", str(outdir),
        ])
        assert result.exit_code == 0, result.output
        plots = list(outdir.glob("*.png"))
        assert len(plots) >= 3

    def test_with_compare_file(
        self, tmp_path: Path, hits: Path, compare: Path, runner: CliRunner
    ) -> None:
        outdir = tmp_path / "out"
        result = runner.invoke(app, [
            "plot",
            "--hits", str(hits),
            "--compare", str(compare),
            "--plot-type", "compare_summary",
            "--dpi", "72",
            "--outdir", str(outdir),
        ])
        assert result.exit_code == 0, result.output
        assert (outdir / "compare_summary.png").exists()

    def test_landscape_plot_set_renders_existing_tables(
        self, indexed_vcf: Path, tmp_path: Path, runner: CliRunner
    ) -> None:
        outdir = tmp_path / "landscape-data"
        result = runner.invoke(app, [
            "landscape",
            "--vcf", str(indexed_vcf),
            "--targets", "T1", "T2",
            "--window-records", "3",
            "--step-records", "3",
            "--outdir", str(outdir),
        ])
        assert result.exit_code == 0, result.output
        assert not (outdir / "missingness_heatmap.svg").exists()

        plot_result = runner.invoke(app, [
            "plot",
            "--plot-set", "landscape",
            "--input-dir", str(outdir),
            "--output-format", "svg",
        ])

        assert plot_result.exit_code == 0, plot_result.output
        assert (outdir / "missingness_heatmap.svg").exists()
        assert (outdir / "private_burden_heatmap.svg").exists()
        assert (outdir / "local_background_map.svg").exists()
        assert (outdir / "similarity_cluster_map.svg").exists()
        data = json.loads((outdir / "landscape.json").read_text())
        assert data["parameters"]["write_plots"] is True
        assert data["parameters"]["plot_format"] == "svg"
        assert "missingness_heatmap.svg" in data["outputs"]

    def test_pangenome_plot_set_renders_existing_tables(
        self, tmp_path: Path, runner: CliRunner
    ) -> None:
        outdir = tmp_path / "pangenome-data"
        result = runner.invoke(app, [
            "pangenome",
            "--gfa", str(GFA_PATH),
            "--targets", "T1", "T2",
            "--permutations", "1",
            "--outdir", str(outdir),
        ])
        assert result.exit_code == 0, result.output
        assert not (outdir / "pangenome_growth.svg").exists()

        plot_result = runner.invoke(app, [
            "plot",
            "--plot-set", "pangenome",
            "--input-dir", str(outdir),
            "--output-format", "svg",
        ])

        assert plot_result.exit_code == 0, plot_result.output
        assert (outdir / "pangenome_growth.svg").exists()
        assert (outdir / "pangenome_coverage.svg").exists()
        assert (outdir / "pangenome_composition.svg").exists()
        data = json.loads((outdir / "pangenome.json").read_text())
        assert data["parameters"]["write_plots"] is True
        assert data["parameters"]["plot_format"] == "svg"
        assert "pangenome_growth.svg" in data["outputs"]
