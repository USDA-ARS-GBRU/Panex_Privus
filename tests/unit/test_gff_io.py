"""Unit tests for src/privy/io/gff.py."""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from privy.io.gff import (
    AnnotationIndex,
    GffRecord,
    build_annotation_index,
    load_contig_alias,
    parse_gff3,
    query_genes,
    query_sub_feature,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

GFF3_FIXTURE = Path(__file__).parent.parent / "data" / "small_cohort.gff3"


def _write_gff(tmp_path: Path, lines: list[str], gz: bool = False) -> Path:
    content = "\n".join(lines) + "\n"
    if gz:
        p = tmp_path / "test.gff3.gz"
        with gzip.open(p, "wt", encoding="utf-8") as fh:
            fh.write(content)
    else:
        p = tmp_path / "test.gff3"
        p.write_text(content, encoding="utf-8")
    return p


_MINIMAL_GFF3 = [
    "##gff-version 3",
    "chr1\t.\tgene\t101\t200\t.\t+\t.\tID=G1;Name=Gene1",
    "chr1\t.\tmRNA\t101\t200\t.\t+\t.\tID=G1.t1;Parent=G1",
    "chr1\t.\texon\t101\t150\t.\t+\t.\tParent=G1.t1",
    "chr1\t.\tCDS\t101\t150\t.\t+\t0\tParent=G1.t1",
    "chr1\t.\tfive_prime_UTR\t151\t160\t.\t+\t.\tParent=G1.t1",
    "chr1\t.\tthree_prime_UTR\t190\t200\t.\t-\t.\tParent=G1.t1",
]


# ---------------------------------------------------------------------------
# Class 1: parse_gff3
# ---------------------------------------------------------------------------

class TestParseGff3:
    def test_parses_gene_record(self, tmp_path):
        p = _write_gff(tmp_path, _MINIMAL_GFF3)
        records = list(parse_gff3(p))
        gene = next(r for r in records if r.feature_type == "gene")
        assert gene.seqid == "chr1"
        assert gene.start == 100   # 1-based 101 → 0-based 100
        assert gene.end == 200     # 1-based 200 inclusive → 0-based exclusive 200
        assert gene.strand == "+"
        assert gene.gene_id == "Gene1"

    def test_coordinate_conversion(self, tmp_path):
        p = _write_gff(tmp_path, _MINIMAL_GFF3)
        records = list(parse_gff3(p))
        cds = next(r for r in records if r.feature_type == "CDS")
        assert cds.start == 100   # 1-based 101 → 0-based 100
        assert cds.end == 150     # 1-based 150 → 0-based exclusive 150

    def test_feature_type_filter(self, tmp_path):
        p = _write_gff(tmp_path, _MINIMAL_GFF3)
        records = list(parse_gff3(p, feature_types=frozenset({"gene"})))
        assert all(r.feature_type == "gene" for r in records)
        assert len(records) == 1

    def test_skips_comment_lines(self, tmp_path):
        p = _write_gff(tmp_path, _MINIMAL_GFF3)
        records = list(parse_gff3(p))
        # No records with empty seqid
        assert all(r.seqid != "" for r in records)

    def test_parses_gz_file(self, tmp_path):
        p = _write_gff(tmp_path, _MINIMAL_GFF3, gz=True)
        records = list(parse_gff3(p))
        assert len(records) >= 1
        assert records[0].feature_type in {"gene", "mRNA", "exon", "CDS",
                                            "five_prime_UTR", "three_prime_UTR"}

    def test_gene_id_from_name_attr(self, tmp_path):
        p = _write_gff(tmp_path, _MINIMAL_GFF3)
        genes = [r for r in parse_gff3(p) if r.feature_type == "gene"]
        assert genes[0].gene_id == "Gene1"

    def test_gene_id_falls_back_to_id_attr(self, tmp_path):
        lines = [
            "##gff-version 3",
            "chr1\t.\tgene\t1\t100\t.\t+\t.\tID=MyGene",
        ]
        p = _write_gff(tmp_path, lines)
        genes = [r for r in parse_gff3(p) if r.feature_type == "gene"]
        assert genes[0].gene_id == "MyGene"

    def test_returns_iterator(self, tmp_path):
        p = _write_gff(tmp_path, _MINIMAL_GFF3)
        result = parse_gff3(p)
        import types
        assert isinstance(result, types.GeneratorType)

    def test_real_fixture_parses(self):
        records = list(parse_gff3(GFF3_FIXTURE))
        feature_types = {r.feature_type for r in records}
        assert "gene" in feature_types
        assert "CDS" in feature_types


# ---------------------------------------------------------------------------
# Class 2: build_annotation_index
# ---------------------------------------------------------------------------

class TestBuildAnnotationIndex:
    def test_genes_populated(self, tmp_path):
        p = _write_gff(tmp_path, _MINIMAL_GFF3)
        idx = build_annotation_index(p)
        assert "chr1" in idx.genes
        assert len(idx.genes["chr1"]) == 1

    def test_gene_coords_correct(self, tmp_path):
        p = _write_gff(tmp_path, _MINIMAL_GFF3)
        idx = build_annotation_index(p)
        g = idx.genes["chr1"][0]
        assert g[0] == 100  # start
        assert g[1] == 200  # end
        assert g[2] == "Gene1"
        assert g[3] == "+"

    def test_cds_indexed(self, tmp_path):
        p = _write_gff(tmp_path, _MINIMAL_GFF3)
        idx = build_annotation_index(p)
        assert "CDS" in idx.sub_features.get("chr1", {})

    def test_utr_indexed(self, tmp_path):
        p = _write_gff(tmp_path, _MINIMAL_GFF3)
        idx = build_annotation_index(p)
        sf = idx.sub_features.get("chr1", {})
        assert "five_prime_UTR" in sf or "three_prime_UTR" in sf

    def test_mrna_not_indexed(self, tmp_path):
        p = _write_gff(tmp_path, _MINIMAL_GFF3)
        idx = build_annotation_index(p)
        sf = idx.sub_features.get("chr1", {})
        assert "mRNA" not in sf

    def test_sorted_by_start(self, tmp_path):
        lines = [
            "##gff-version 3",
            "chr1\t.\tgene\t200\t300\t.\t+\t.\tID=G2;Name=Gene2",
            "chr1\t.\tgene\t50\t100\t.\t+\t.\tID=G1;Name=Gene1",
        ]
        p = _write_gff(tmp_path, lines)
        idx = build_annotation_index(p)
        starts = [g[0] for g in idx.genes["chr1"]]
        assert starts == sorted(starts)

    def test_multiple_contigs(self, tmp_path):
        lines = [
            "##gff-version 3",
            "chr1\t.\tgene\t1\t100\t.\t+\t.\tID=G1;Name=G1",
            "chr2\t.\tgene\t1\t100\t.\t+\t.\tID=G2;Name=G2",
        ]
        p = _write_gff(tmp_path, lines)
        idx = build_annotation_index(p)
        assert "chr1" in idx.genes
        assert "chr2" in idx.genes

    def test_real_fixture(self):
        idx = build_annotation_index(GFF3_FIXTURE)
        assert "chr1" in idx.genes
        assert len(idx.genes["chr1"]) >= 2


# ---------------------------------------------------------------------------
# Class 3: query_genes and query_sub_feature
# ---------------------------------------------------------------------------

class TestQueryFunctions:
    def setup_method(self):
        self.idx = build_annotation_index(GFF3_FIXTURE)

    def test_query_genes_hit(self):
        # VCF pos 100 → 0-based [99, 100); gene covers [50, 350)
        result = query_genes(self.idx, "chr1", 99, 100)
        assert len(result) >= 1
        assert result[0][2] == "GeneA"

    def test_query_genes_miss(self):
        # Position 700 → 0-based [699, 700); no gene
        result = query_genes(self.idx, "chr1", 699, 700)
        assert result == []

    def test_query_genes_unknown_contig(self):
        result = query_genes(self.idx, "chrX", 0, 100)
        assert result == []

    def test_query_cds_hit(self):
        # VCF pos 100 → 0-based [99, 100); CDS is [90, 110)
        assert query_sub_feature(self.idx, "chr1", "CDS", 99, 100) is True

    def test_query_cds_miss(self):
        # Position 250 → 0-based [249, 250); no CDS there
        assert query_sub_feature(self.idx, "chr1", "CDS", 249, 250) is False

    def test_query_utr_hit(self):
        # VCF pos 200 → 0-based [199, 200); five_prime_UTR covers [190, 210)
        assert query_sub_feature(self.idx, "chr1", "five_prime_UTR", 199, 200) is True

    def test_query_utr_miss(self):
        assert query_sub_feature(self.idx, "chr1", "five_prime_UTR", 699, 700) is False

    def test_query_unknown_feature_type(self):
        assert query_sub_feature(self.idx, "chr1", "nonexistent_feat", 99, 100) is False


# ---------------------------------------------------------------------------
# Class 4: load_contig_alias
# ---------------------------------------------------------------------------

class TestLoadContigAlias:
    def test_basic_load(self, tmp_path):
        p = tmp_path / "alias.tsv"
        p.write_text("Gm01\tchr1\nGm02\tchr2\n", encoding="utf-8")
        alias = load_contig_alias(p)
        assert alias == {"Gm01": "chr1", "Gm02": "chr2"}

    def test_skips_comment_lines(self, tmp_path):
        p = tmp_path / "alias.tsv"
        p.write_text("# header\nGm01\tchr1\n", encoding="utf-8")
        alias = load_contig_alias(p)
        assert "# header" not in alias
        assert alias == {"Gm01": "chr1"}

    def test_skips_blank_lines(self, tmp_path):
        p = tmp_path / "alias.tsv"
        p.write_text("Gm01\tchr1\n\nGm02\tchr2\n", encoding="utf-8")
        alias = load_contig_alias(p)
        assert len(alias) == 2

    def test_empty_file(self, tmp_path):
        p = tmp_path / "alias.tsv"
        p.write_text("", encoding="utf-8")
        alias = load_contig_alias(p)
        assert alias == {}
