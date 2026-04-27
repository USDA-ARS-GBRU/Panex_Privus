"""Unit tests for privy plot functions."""

from __future__ import annotations

from pathlib import Path

# ---------------------------------------------------------------------------
# Minimal test data
# ---------------------------------------------------------------------------

def _hits_rows(n: int = 5) -> list[dict[str, str]]:
    strictness_cycle = [
        "strict_complete", "strict_target_missing",
        "strict_offtarget_missing", "relaxed_threshold", "contradicted",
    ]
    return [
        {
            "locus_id": f"PPX{i:06d}",
            "contig": "chr1",
            "start": str(i * 1000),
            "end": str(i * 1000 + 100),
            "variant_type": "snp",
            "allele_key": f"chr1:{i*1000}:A:T",
            "strictness_class": strictness_cycle[i % len(strictness_cycle)],
            "final_score": str(round(1.0 - i * 0.1, 2)),
        }
        for i in range(n)
    ]


def _evidence_rows() -> list[dict[str, str]]:
    rows = []
    for src in ("vcf", "bam"):
        for ec in ("support", "absence", "uninformative"):
            rows.append({
                "locus_id": "PPX000001",
                "source_type": src,
                "evidence_class": ec,
                "metric_name": "depth",
                "metric_value": "12.0",
                "sample_id": "S1",
                "details": "",
            })
    return rows


def _compare_rows() -> list[dict[str, str]]:
    match_classes = ["supported", "partially_supported", "source_specific"]
    return [
        {
            "compare_id": f"CMP{i:06d}",
            "locus_id_a": f"PPX{i:06d}",
            "locus_id_b": f"GPX{i:06d}" if mc != "source_specific" else "NA",
            "source_a": "vcf",
            "source_b": "gfa",
            "match_class": mc,
            "coordinate_overlap": "0.9" if mc != "source_specific" else "0.0",
            "comparison_score": "0.9",
        }
        for i, mc in enumerate(match_classes)
    ]


def _pangenome_growth_rows() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for group, base in (("full", 5), ("target", 3), ("off_target", 4)):
        for trial in (1, 2):
            for n in (1, 2, 3):
                rows.append({
                    "group": group,
                    "trial": trial,
                    "n": n,
                    "sample_added": f"S{n}",
                    "features": base + n + trial,
                    "bp": (base + n + trial) * 10,
                    "new_features": 1,
                    "new_bp": 10,
                    "singleton_features": 2,
                    "singleton_bp": 20,
                })
    return rows


def _pangenome_coverage_rows() -> list[dict[str, object]]:
    return [
        {"group": group, "coverage": coverage, "n_features": 4 - coverage, "n_bp": 10}
        for group in ("full", "target", "off_target")
        for coverage in range(4)
    ]


def _pangenome_composition_rows() -> list[dict[str, object]]:
    return [
        {"group": group, "category": category, "n_features": i + 1, "n_bp": (i + 1) * 10}
        for group in ("full", "target", "off_target")
        for i, category in enumerate(("core", "accessory", "private", "absent"))
    ]


# ---------------------------------------------------------------------------
# TestPlotStrictnessBar
# ---------------------------------------------------------------------------

class TestPlotStrictnessBar:
    def test_creates_file(self, tmp_path: Path) -> None:
        from privy.plot.summaries import plot_strictness_bar
        out = plot_strictness_bar(_hits_rows(), tmp_path, output_format="png")
        assert out.exists()
        assert out.stat().st_size > 500

    def test_output_path_has_correct_name(self, tmp_path: Path) -> None:
        from privy.plot.summaries import plot_strictness_bar
        out = plot_strictness_bar(_hits_rows(), tmp_path)
        assert out.name == "strictness_bar.png"

    def test_svg_format(self, tmp_path: Path) -> None:
        from privy.plot.summaries import plot_strictness_bar
        out = plot_strictness_bar(_hits_rows(), tmp_path, output_format="svg")
        assert out.suffix == ".svg"
        assert out.exists()

    def test_empty_rows(self, tmp_path: Path) -> None:
        from privy.plot.summaries import plot_strictness_bar
        out = plot_strictness_bar([], tmp_path)
        assert out.exists()

    def test_single_class(self, tmp_path: Path) -> None:
        from privy.plot.summaries import plot_strictness_bar
        rows = [{"strictness_class": "strict_complete", "final_score": "1.0"}]
        out = plot_strictness_bar(rows, tmp_path)
        assert out.exists()


# ---------------------------------------------------------------------------
# TestPlotScoreDistribution
# ---------------------------------------------------------------------------

class TestPlotScoreDistribution:
    def test_creates_file(self, tmp_path: Path) -> None:
        from privy.plot.summaries import plot_score_distribution
        out = plot_score_distribution(_hits_rows(), tmp_path)
        assert out.exists()
        assert out.stat().st_size > 500

    def test_output_name(self, tmp_path: Path) -> None:
        from privy.plot.summaries import plot_score_distribution
        out = plot_score_distribution(_hits_rows(), tmp_path)
        assert out.name == "score_distribution.png"

    def test_empty_rows_creates_placeholder(self, tmp_path: Path) -> None:
        from privy.plot.summaries import plot_score_distribution
        out = plot_score_distribution([], tmp_path)
        assert out.exists()

    def test_malformed_score_skipped(self, tmp_path: Path) -> None:
        from privy.plot.summaries import plot_score_distribution
        rows = [
            {"strictness_class": "strict_complete", "final_score": "not_a_float"},
            {"strictness_class": "strict_complete", "final_score": "0.9"},
        ]
        out = plot_score_distribution(rows, tmp_path)
        assert out.exists()


# ---------------------------------------------------------------------------
# TestPlotSupportBar
# ---------------------------------------------------------------------------

class TestPlotSupportBar:
    def test_creates_file(self, tmp_path: Path) -> None:
        from privy.plot.summaries import plot_support_bar
        out = plot_support_bar(_evidence_rows(), tmp_path)
        assert out.exists()
        assert out.stat().st_size > 500

    def test_output_name(self, tmp_path: Path) -> None:
        from privy.plot.summaries import plot_support_bar
        out = plot_support_bar(_evidence_rows(), tmp_path)
        assert out.name == "support_bar.png"

    def test_empty_rows_creates_placeholder(self, tmp_path: Path) -> None:
        from privy.plot.summaries import plot_support_bar
        out = plot_support_bar([], tmp_path)
        assert out.exists()


# ---------------------------------------------------------------------------
# TestPlotCompareSummary
# ---------------------------------------------------------------------------

class TestPlotCompareSummary:
    def test_creates_file(self, tmp_path: Path) -> None:
        from privy.plot.summaries import plot_compare_summary
        out = plot_compare_summary(_compare_rows(), tmp_path)
        assert out.exists()
        assert out.stat().st_size > 500

    def test_output_name(self, tmp_path: Path) -> None:
        from privy.plot.summaries import plot_compare_summary
        out = plot_compare_summary(_compare_rows(), tmp_path)
        assert out.name == "compare_summary.png"

    def test_empty_rows_creates_placeholder(self, tmp_path: Path) -> None:
        from privy.plot.summaries import plot_compare_summary
        out = plot_compare_summary([], tmp_path)
        assert out.exists()


# ---------------------------------------------------------------------------
# TestPlotLocusPanel
# ---------------------------------------------------------------------------

class TestPlotLocusPanel:
    def test_creates_file(self, tmp_path: Path) -> None:
        from privy.plot.loci import plot_locus_panel
        out = plot_locus_panel(_hits_rows(10), tmp_path)
        assert out.exists()
        assert out.stat().st_size > 500

    def test_output_name(self, tmp_path: Path) -> None:
        from privy.plot.loci import plot_locus_panel
        out = plot_locus_panel(_hits_rows(), tmp_path)
        assert out.name == "locus_panel.png"

    def test_top_n_respected(self, tmp_path: Path) -> None:
        from privy.plot.loci import plot_locus_panel
        out = plot_locus_panel(_hits_rows(20), tmp_path, top_n=5)
        assert out.exists()

    def test_empty_rows_creates_placeholder(self, tmp_path: Path) -> None:
        from privy.plot.loci import plot_locus_panel
        out = plot_locus_panel([], tmp_path)
        assert out.exists()

    def test_svg_format(self, tmp_path: Path) -> None:
        from privy.plot.loci import plot_locus_panel
        out = plot_locus_panel(_hits_rows(), tmp_path, output_format="svg")
        assert out.suffix == ".svg"
        assert out.exists()

    def test_show_labels_false(self, tmp_path: Path) -> None:
        from privy.plot.loci import plot_locus_panel
        out = plot_locus_panel(_hits_rows(5), tmp_path, show_labels=False)
        assert out.exists()


# ---------------------------------------------------------------------------
# TestPangenomePlots
# ---------------------------------------------------------------------------

class TestPangenomePlots:
    def test_growth_plot_creates_file(self, tmp_path: Path) -> None:
        from privy.plot.pangenome import plot_pangenome_growth
        out = plot_pangenome_growth(_pangenome_growth_rows(), tmp_path)
        assert out.exists()
        assert out.stat().st_size > 500

    def test_coverage_plot_creates_file(self, tmp_path: Path) -> None:
        from privy.plot.pangenome import plot_pangenome_coverage
        out = plot_pangenome_coverage(_pangenome_coverage_rows(), tmp_path)
        assert out.exists()
        assert out.stat().st_size > 500

    def test_composition_plot_creates_file(self, tmp_path: Path) -> None:
        from privy.plot.pangenome import plot_pangenome_composition
        out = plot_pangenome_composition(_pangenome_composition_rows(), tmp_path)
        assert out.exists()
        assert out.stat().st_size > 500

    def test_all_pangenome_plots_creates_three_files(self, tmp_path: Path) -> None:
        from privy.plot.pangenome import plot_all_pangenome
        generated = plot_all_pangenome(
            coverage_rows=_pangenome_coverage_rows(),
            composition_rows=_pangenome_composition_rows(),
            growth_rows=_pangenome_growth_rows(),
            outdir=tmp_path,
        )
        assert {p.name for p in generated} == {
            "pangenome_growth.png",
            "pangenome_coverage.png",
            "pangenome_composition.png",
        }
