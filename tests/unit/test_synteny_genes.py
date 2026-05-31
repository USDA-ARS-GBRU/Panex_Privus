"""Unit tests for src/privy/synteny/genes.py (gene-anchored synteny)."""

from __future__ import annotations

from privy.synteny.chain import ChainParams
from privy.synteny.genes import GeneRow, build_gene_synteny, gene_anchors
from privy.synteny.model import AnchorSource, BlockType

# 5 orthogroups; query genes ascending on chrQ.
_QUERY = [
    GeneRow(f"q{i}", "Q", "chrQ", i * 1000, i * 1000 + 500, "+")
    for i in range(1, 6)
]
_ORTHO = {f"q{i}": f"og{i}" for i in range(1, 6)}
_ORTHO.update({f"t{i}": f"og{i}" for i in range(1, 6)})

PARAMS = ChainParams(min_anchors=3, max_gap=10_000)


class TestGeneAnchors:
    def test_builds_orthogroup_pairs(self):
        target = [GeneRow(f"t{i}", "T", "chrT", i * 1000, i * 1000 + 500, "+") for i in range(1, 6)]
        anchors = gene_anchors(_QUERY, target, _ORTHO)
        assert len(anchors) == 5
        assert all(a.source is AnchorSource.GENE for a in anchors)
        assert {a.name for a in anchors} == {f"og{i}" for i in range(1, 6)}

    def test_unmapped_genes_skipped(self):
        target = [GeneRow("tX", "T", "chrT", 0, 500, "+")]   # no orthogroup
        assert gene_anchors(_QUERY, target, _ORTHO) == []


class TestBuildGeneSynteny:
    def test_collinear_block(self):
        target = [GeneRow(f"t{i}", "T", "chrT", i * 1000, i * 1000 + 500, "+") for i in range(1, 6)]
        blocks = build_gene_synteny(_QUERY, target, _ORTHO, PARAMS)
        assert len(blocks) == 1
        assert blocks[0].block_type is BlockType.COLLINEAR
        assert blocks[0].n_anchors == 5

    def test_inversion_block(self):
        # target genes in reverse genomic order AND reverse strand -> inversion
        target = [
            GeneRow(f"t{i}", "T", "chrT", (6 - i) * 1000, (6 - i) * 1000 + 500, "-")
            for i in range(1, 6)
        ]
        blocks = build_gene_synteny(_QUERY, target, _ORTHO, PARAMS)
        assert len(blocks) == 1
        assert blocks[0].block_type is BlockType.INVERSION
        assert blocks[0].strand == "-"

    def test_min_anchors_filter(self):
        target = [GeneRow(f"t{i}", "T", "chrT", i * 1000, i * 1000 + 500, "+") for i in range(1, 6)]
        blocks = build_gene_synteny(_QUERY, target, _ORTHO, ChainParams(min_anchors=10))
        assert blocks == []
