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
