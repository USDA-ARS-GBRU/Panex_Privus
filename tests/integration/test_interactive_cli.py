"""Integration tests for ``privy interactive``."""

from __future__ import annotations

from pathlib import Path

from typer.testing import CliRunner

from privy.cli.main import app


def _write_sites(path: Path) -> Path:
    path.write_text(
        "\t".join(
            [
                "contig",
                "pos",
                "ref",
                "alt",
                "variant_type",
                "Harosoy_gt",
                "Harosoy-sharp_gt",
                "Kingawa_gt",
                "target_private_alt_pattern",
                "overlapping_gene_names",
                "nearest_gene_name",
            ]
        )
        + "\n"
        + "chr1\t120\tA\tG\tsnp\t0/0\t1/1\t1/1\ttrue\tGeneA\tGeneA\n"
        + "chr1\t180\tT\tTA\tnon_snp\t0/0\t1/1\t1/1\ttrue\tGeneA\tGeneA\n"
        + "chr1\t600\tC\tT\tsnp\t0/1\t0/1\t1/1\tfalse\tNA\tGeneB\n",
        encoding="utf-8",
    )
    return path


def _write_gff3(path: Path) -> Path:
    path.write_text(
        "##gff-version 3\n"
        "chr1\tprivy\tgene\t100\t300\t.\t+\t.\tID=geneA;Name=GeneA\n"
        "chr1\tprivy\tmRNA\t100\t300\t.\t+\t.\tID=txA;Parent=geneA;longest=1\n"
        "chr1\tprivy\texon\t100\t150\t.\t+\t.\tID=exon1;Parent=txA\n"
        "chr1\tprivy\texon\t220\t300\t.\t+\t.\tID=exon2;Parent=txA\n"
        "chr1\tprivy\tCDS\t120\t145\t.\t+\t0\tID=cds1;Parent=txA\n"
        "chr1\tprivy\tgene\t500\t800\t.\t-\t.\tID=geneB;Name=GeneB\n",
        encoding="utf-8",
    )
    return path


def _write_track(path: Path) -> Path:
    path.write_text(
        "##gff-version 3\n"
        "chr1\tprivy\trepeat_region\t110\t190\t.\t+\t.\tID=repeat1;Name=RepeatOne;class=Simple_repeat\n",
        encoding="utf-8",
    )
    return path


def _write_functional_tsv(path: Path) -> Path:
    path.write_text(
        "gene\tgene_id\trepresentative_predicted_function\t"
        "functional_category_keywords\tscreening_priority\tscreening_note\n"
        "GeneA\tgeneA\tAuxin-responsive trichome regulator\t"
        "trichome_or_epidermal_development\thigh_screening_interest\t"
        "Variant-supported trichome candidate\n",
        encoding="utf-8",
    )
    return path


def _write_vcf(path: Path) -> Path:
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=1000>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        "Harosoy\tHarosoy-sharp\tKingawa\n"
        "chr1\t120\t.\tA\tG\t60\tPASS\t.\tGT\t0/0\t1/1\t1/1\n"
        "chr1\t180\t.\tT\tTA\t60\tPASS\t.\tGT\t0/0\t1/1\t1/1\n"
        "chr1\t600\t.\tC\tT\t60\tPASS\t.\tGT\t0/1\t0/1\t1/1\n",
        encoding="utf-8",
    )
    return path


def _write_scan_source(path: Path, prefix: str) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    path.joinpath("hits.tsv").write_text(
        "locus_id\tcontig\tstart\tend\tvariant_type\tallele_key\t"
        "target_support_n\ttarget_total_n\tofftarget_support_n\tofftarget_total_n\t"
        "target_missing_n\tofftarget_missing_n\tstrictness_class\t"
        "discovery_score\tsupport_score\tpenalty_score\tfinal_score\n"
        f"{prefix}001\tchr1\t99\t100\tsnp\tchr1:100:A:G\t2\t2\t0\t3\t0\t0\t"
        "strict_complete\t1.0\t0.2\t0.0\t1.2\n"
        f"{prefix}002\tchr2\t199\t205\tindel\tchr2:200:AT:A\t1\t2\t0\t3\t1\t0\t"
        "strict_target_missing\t0.8\t0.0\t0.1\t0.7\n",
        encoding="utf-8",
    )
    path.joinpath("regions.tsv").write_text(
        "region_id\tcontig\tstart\tend\tn_loci\tvariant_types\t"
        "dominant_strictness_class\ttarget_consistency\tofftarget_exclusion\tfinal_score\n"
        f"{prefix}R1\tchr1\t90\t120\t2\tsnp\tstrict_complete\t1.0\t1.0\t1.2\n",
        encoding="utf-8",
    )
    path.joinpath("qc.tsv").write_text(
        "metric\tvalue\tdescription\n"
        "records_evaluated\t20\tTotal records processed\n"
        "loci_emitted\t2\tLoci written to hits.tsv\n",
        encoding="utf-8",
    )
    path.joinpath("evidence.tsv").write_text(
        "locus_id\tsource_type\tsample_id\tevidence_class\tmetric_name\tmetric_value\tdetails\n"
        f"{prefix}001\tvcf\t\tsupport\tallele_pattern\t1.0\tpasses\n",
        encoding="utf-8",
    )
    return path


def _write_landscape(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    path.joinpath("windows.tsv").write_text(
        "window_id\tcontig\twindow_index\tstart\tend\tmidpoint\twindow_mode\t"
        "n_variants\tspan_bp\tdensity_variants_per_kb\ttarget_mean_missing_rate\t"
        "offtarget_mean_missing_rate\ttarget_mean_nonref_rate\tofftarget_mean_nonref_rate\t"
        "target_private_alt_n\tofftarget_private_alt_n\ttarget_private_alt_rate\t"
        "offtarget_private_alt_rate\ttop_nearest_background\ttop_nearest_background_n\n"
        "LW1\tchr1\t1\t0\t100\t50\tbp\t10\t100\t0.100\t0.000\t0.000\t"
        "0.500\t0.100\t2\t0\t0.200\t0.000\tOffA\t1\n"
        "LW2\tchr1\t2\t100\t200\t150\tbp\t20\t100\t0.200\t0.100\t0.000\t"
        "0.600\t0.200\t1\t0\t0.050\t0.000\tOffA\t1\n"
        "LW3\tchr2\t1\t0\t100\t50\tbp\t5\t100\t0.050\t0.000\t0.000\t"
        "0.200\t0.300\t0\t1\t0.000\t0.200\tTargetA\t1\n",
        encoding="utf-8",
    )
    path.joinpath("sample_windows.tsv").write_text(
        "window_id\tcontig\twindow_index\tstart\tend\tmidpoint\twindow_mode\tn_variants\t"
        "sample\tcohort_role\tcalled_n\tmissing_n\tmissing_rate\thet_n\thet_rate\t"
        "nonref_n\tnonref_rate\tminor_genotype_n\tminor_genotype_rate\trare_alt_n\t"
        "rare_alt_rate\tprivate_alt_n\tprivate_alt_rate\tmedian_call_freq\t"
        "nearest_background\tnearest_background_role\tnearest_similarity\t"
        "similarity_compared_variants\n"
        "LW1\tchr1\t1\t0\t100\t50\tbp\t10\tTargetA\ttarget\t10\t0\t0.000\t"
        "0\t0.000\t5\t0.500\t0\t0.000\t0\t0.000\t2\t0.200\tNA\tOffA\toff_target\t0.900\t10\n"
        "LW1\tchr1\t1\t0\t100\t50\tbp\t10\tOffA\toff_target\t10\t0\t0.000\t"
        "0\t0.000\t1\t0.100\t0\t0.000\t0\t0.000\t0\t0.000\tNA\tTargetA\ttarget\t0.900\t10\n"
        "LW2\tchr1\t2\t100\t200\t150\tbp\t20\tTargetA\ttarget\t18\t2\t0.100\t"
        "0\t0.000\t12\t0.600\t0\t0.000\t0\t0.000\t1\t0.050\tNA\tOffA\toff_target\t0.800\t18\n"
        "LW2\tchr1\t2\t100\t200\t150\tbp\t20\tOffA\toff_target\t20\t0\t0.000\t"
        "0\t0.000\t4\t0.200\t0\t0.000\t0\t0.000\t0\t0.000\tNA\tTargetA\ttarget\t0.800\t20\n",
        encoding="utf-8",
    )
    path.joinpath("candidate_introgression_blocks.tsv").write_text(
        "block_id\tsample\tcontig\tstart\tend\tn_windows\tcandidate_donor\t"
        "candidate_donor_role\tmean_donor_similarity\tmean_nearest_target_similarity\t"
        "mean_similarity_delta\tmax_missing_rate\tmean_private_alt_rate\t"
        "mean_nonref_rate\tevidence_class\tinterpretation\n"
        "IB1\tTargetA\tchr1\t0\t200\t2\tOffA\toff_target\t0.850\t0.300\t"
        "0.550\t0.100\t0.125\t0.550\tofftarget_closer_than_target\t"
        "Target sample is locally closest to an off-target sample.\n",
        encoding="utf-8",
    )
    path.joinpath("background_blocks.tsv").write_text(
        "block_id\tsample\tcohort_role\tcontig\tstart\tend\tn_windows\t"
        "nearest_background\tnearest_background_role\tmean_similarity\n"
        "LB1\tTargetA\ttarget\tchr1\t0\t200\t2\tOffA\toff_target\t0.850\n",
        encoding="utf-8",
    )
    path.joinpath("filter_summary.tsv").write_text(
        "metric\tvalue\tdescription\n"
        "records_seen\t3\tRecords considered\n",
        encoding="utf-8",
    )
    path.joinpath("similarity.tsv").write_text(
        "window_id\tcontig\twindow_index\tstart\tend\tsample_a\tsample_b\t"
        "similarity\tcompared_variants\n"
        "LW1\tchr1\t1\t0\t100\tTargetA\tOffA\t0.900\t10\n"
        "LW2\tchr1\t2\t100\t200\tTargetA\tOffA\t0.800\t18\n",
        encoding="utf-8",
    )
    path.joinpath("landscape.json").write_text(
        '{"analysis":"landscape","parameters":{"window_bp":100,"step_bp":100}}\n',
        encoding="utf-8",
    )
    return path


def _write_pangenome(path: Path, source_type: str) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    prefix = source_type.upper()
    path.joinpath("feature_summary.tsv").write_text(
        "feature_id\tsource_type\tfeature_type\tcontig\tstart\tend\tlength\t"
        "total_present_n\ttarget_present_n\ttarget_total_n\tofftarget_present_n\t"
        "offtarget_total_n\tfull_category\ttarget_category\tofftarget_category\t"
        "target_private\tofftarget_private\n"
        f"{prefix}F1\t{source_type}\tsegment\tchr1\t0\t100\t100\t2\t2\t2\t0\t3\t"
        "accessory\tcore\tabsent\tTrue\tFalse\n"
        f"{prefix}F2\t{source_type}\tsegment\tchr1\t100\t160\t60\t5\t2\t2\t3\t3\t"
        "core\tcore\tcore\tFalse\tFalse\n"
        f"{prefix}F3\t{source_type}\tsegment\tchr2\t0\t40\t40\t3\t0\t2\t3\t3\t"
        "accessory\tabsent\tcore\tFalse\tTrue\n",
        encoding="utf-8",
    )
    path.joinpath("composition.tsv").write_text(
        "group\tcategory\tn_features\tn_bp\n"
        "full\tabsent\t0\t0\n"
        "full\tprivate\t0\t0\n"
        "full\taccessory\t2\t140\n"
        "full\tcore\t1\t60\n"
        "target\tabsent\t1\t40\n"
        "target\tprivate\t0\t0\n"
        "target\taccessory\t0\t0\n"
        "target\tcore\t2\t160\n"
        "off_target\tabsent\t1\t100\n"
        "off_target\tprivate\t0\t0\n"
        "off_target\taccessory\t0\t0\n"
        "off_target\tcore\t2\t100\n",
        encoding="utf-8",
    )
    path.joinpath("coverage_histogram.tsv").write_text(
        "group\tcoverage\tn_features\tn_bp\n"
        "full\t0\t0\t0\n"
        "full\t1\t0\t0\n"
        "full\t2\t1\t100\n"
        "full\t3\t1\t40\n"
        "target\t0\t1\t40\n"
        "target\t1\t0\t0\n"
        "target\t2\t2\t160\n",
        encoding="utf-8",
    )
    path.joinpath("growth_curves.tsv").write_text(
        "group\ttrial\tn\tsample_added\tfeatures\tbp\tnew_features\tnew_bp\t"
        "singleton_features\tsingleton_bp\n"
        "full\t1\t1\tT1\t2\t160\t2\t160\t2\t160\n"
        "full\t1\t2\tO1\t3\t200\t1\t40\t3\t200\n"
        "target\t1\t1\tT1\t2\t160\t2\t160\t2\t160\n"
        "target\t1\t2\tT2\t2\t160\t0\t0\t0\t0\n",
        encoding="utf-8",
    )
    path.joinpath("pangenome.json").write_text(
        '{"analysis":"pangenome","source_type":"'
        + source_type
        + '","summary":{"n_features":3,"n_samples":5,'
        '"n_target_samples":2,"n_offtarget_samples":3},'
        '"samples":{"target":["T1","T2"],"off_target":["O1","O2","O3"],'
        '"full":["T1","T2","O1","O2","O3"]},'
        '"parameters":{"permutations":1}}\n',
        encoding="utf-8",
    )
    return path


def test_interactive_focus_writes_one_html_per_region(tmp_path: Path) -> None:
    sites = _write_sites(tmp_path / "sites.tsv")
    gff3 = _write_gff3(tmp_path / "genes.gff3")
    track = _write_track(tmp_path / "repeats.gff3")
    outdir = tmp_path / "interactive"

    result = CliRunner().invoke(
        app,
        [
            "interactive",
            "--focus",
            "chr1:1-1000",
            "--sites-tsv",
            str(sites),
            "--gff3",
            str(gff3),
            "--samples",
            "Harosoy",
            "Harosoy-sharp",
            "Kingawa",
            "--track-gff",
            f"Repeats={track}",
            "--sample-abbrev",
            "HS=Harosoy-sharp",
            "--outdir",
            str(outdir),
        ],
    )

    assert result.exit_code == 0, result.output
    html = outdir / "focus_chr1_1_1000.html"
    features = outdir / "focus_chr1_1_1000.features.tsv"
    metadata = outdir / "focus_chr1_1_1000.json"
    assert html.exists()
    assert features.exists()
    assert metadata.exists()
    assert not (outdir / "index.html").exists()
    text = html.read_text(encoding="utf-8")
    assert "Interactive Genome Browser" in text
    assert "RepeatOne" in text
    assert "GeneA" in features.read_text(encoding="utf-8")


def test_interactive_focus_writes_index_for_multiple_regions(tmp_path: Path) -> None:
    sites = _write_sites(tmp_path / "sites.tsv")
    gff3 = _write_gff3(tmp_path / "genes.gff3")
    outdir = tmp_path / "interactive"

    result = CliRunner().invoke(
        app,
        [
            "interactive",
            "--focus",
            "chr1:1-400",
            "--focus",
            "chr1:401-1000",
            "--sites-tsv",
            str(sites),
            "--gff3",
            str(gff3),
            "--samples",
            "Harosoy",
            "Harosoy-sharp",
            "Kingawa",
            "--outdir",
            str(outdir),
        ],
    )

    assert result.exit_code == 0, result.output
    assert (outdir / "index.html").exists()
    assert (outdir / "focus_chr1_1_400.html").exists()
    assert (outdir / "focus_chr1_401_1000.html").exists()


def test_interactive_focus_extracts_sites_from_vcf(tmp_path: Path) -> None:
    vcf = _write_vcf(tmp_path / "cohort.vcf")
    gff3 = _write_gff3(tmp_path / "genes.gff3")
    outdir = tmp_path / "interactive"

    result = CliRunner().invoke(
        app,
        [
            "interactive",
            "--focus",
            "chr1:1-1000",
            "--vcf",
            str(vcf),
            "--gff3",
            str(gff3),
            "--samples",
            "Harosoy",
            "Harosoy-sharp",
            "Kingawa",
            "--outdir",
            str(outdir),
        ],
    )

    assert result.exit_code == 0, result.output
    sites = outdir / "focus_chr1_1_1000.sites.tsv"
    assert sites.exists()
    assert "target_private_alt_pattern" in sites.read_text(encoding="utf-8")
    assert (outdir / "focus_chr1_1_1000.html").exists()


def test_interactive_keyword_group_adds_candidate_dropdown_group(tmp_path: Path) -> None:
    sites = _write_sites(tmp_path / "sites.tsv")
    gff3 = _write_gff3(tmp_path / "genes.gff3")
    functional = _write_functional_tsv(tmp_path / "functional.tsv")
    outdir = tmp_path / "interactive"

    result = CliRunner().invoke(
        app,
        [
            "interactive",
            "--focus",
            "chr1:1-1000",
            "--sites-tsv",
            str(sites),
            "--gff3",
            str(gff3),
            "--functional-tsv",
            str(functional),
            "--keyword-group",
            "Trichome=trichome,auxin",
            "--samples",
            "Harosoy",
            "Harosoy-sharp",
            "Kingawa",
            "--outdir",
            str(outdir),
        ],
    )

    assert result.exit_code == 0, result.output
    html = (outdir / "focus_chr1_1_1000.html").read_text(encoding="utf-8")
    assert '"title":"Trichome"' in html
    assert "Auxin-responsive trichome regulator" in html


def test_interactive_scan_writes_dashboard_from_combined_scan_dir(tmp_path: Path) -> None:
    scan = tmp_path / "scan"
    _write_scan_source(scan / "vcf", "V")
    _write_scan_source(scan / "gfa", "G")
    compare = scan / "compare"
    compare.mkdir()
    compare.joinpath("compare.tsv").write_text(
        "a_locus_id\tb_locus_id\tmatch_class\n"
        "V001\tG001\tsupported\n"
        "V002\t\tmissing_data\n",
        encoding="utf-8",
    )
    outdir = tmp_path / "interactive"

    result = CliRunner().invoke(
        app,
        [
            "interactive",
            "--scan",
            str(scan),
            "--max-hits",
            "1",
            "--max-regions",
            "1",
            "--outdir",
            str(outdir),
        ],
    )

    assert result.exit_code == 0, result.output
    html = outdir / "scan_dashboard.html"
    metadata = outdir / "scan_dashboard.json"
    assert html.exists()
    assert metadata.exists()
    text = html.read_text(encoding="utf-8")
    assert "Privy Interactive Scan Dashboard" in text
    assert "Strictness Classes" in text
    assert "V001" in text
    assert "G001" in text
    assert "supported" in text


def test_interactive_landscape_writes_dashboard(tmp_path: Path) -> None:
    landscape = _write_landscape(tmp_path / "landscape")
    outdir = tmp_path / "interactive"

    result = CliRunner().invoke(
        app,
        [
            "interactive",
            "--landscape",
            str(landscape),
            "--max-windows",
            "10",
            "--max-sample-windows",
            "10",
            "--outdir",
            str(outdir),
        ],
    )

    assert result.exit_code == 0, result.output
    html = outdir / "landscape_dashboard.html"
    metadata = outdir / "landscape_dashboard.json"
    assert html.exists()
    assert metadata.exists()
    text = html.read_text(encoding="utf-8")
    assert "Privy Interactive Landscape Dashboard" in text
    assert "Sample-By-Window Heatmap" in text
    assert "IB1" in text
    assert "TargetA" in text


def test_interactive_pangenome_writes_dashboard_from_combined_dir(tmp_path: Path) -> None:
    pangenome = tmp_path / "pangenome"
    _write_pangenome(pangenome / "vcf", "vcf")
    _write_pangenome(pangenome / "gfa", "gfa")
    outdir = tmp_path / "interactive"

    result = CliRunner().invoke(
        app,
        [
            "interactive",
            "--pangenome",
            str(pangenome),
            "--max-features",
            "2",
            "--max-private-features",
            "1",
            "--outdir",
            str(outdir),
        ],
    )

    assert result.exit_code == 0, result.output
    html = outdir / "pangenome_dashboard.html"
    metadata = outdir / "pangenome_dashboard.json"
    assert html.exists()
    assert metadata.exists()
    text = html.read_text(encoding="utf-8")
    assert "Privy Interactive Pangenome Dashboard" in text
    assert "Composition" in text
    assert "VCFF1" in text
    assert "GFAF1" in text


def test_interactive_requires_sites_tsv_or_vcf() -> None:
    result = CliRunner().invoke(app, ["interactive", "--focus", "chr1:1-1000"])

    assert result.exit_code == 1
    assert "Provide either --sites-tsv or --vcf" in result.output


def test_interactive_requires_one_dashboard_mode(tmp_path: Path) -> None:
    scan = _write_scan_source(tmp_path / "scan", "S")

    result = CliRunner().invoke(
        app,
        [
            "interactive",
            "--focus",
            "chr1:1-1000",
            "--scan",
            str(scan),
        ],
    )

    assert result.exit_code == 2
    assert "Provide exactly one dashboard mode" in result.output
