"""Gene-anchored synteny: orthogroup-constrained gene pairs → collinear blocks.

The classic GENESPACE/MCScanX route: anchor two genomes by genes that share an
orthogroup, then chain those anchors into synteny blocks.  This module builds the
gene-pair :class:`~privy.synteny.model.Anchor` objects; the existing
:func:`~privy.synteny.chain.chain_anchors` does the chaining (so graph-, PAF-, and
gene-anchored modes all share one chainer).

Orthogroups can come from the pangenome graph (shared segments → free), an ingested
OrthoFinder run, or any gene→group mapping.  Gene coordinates can be read with
``privy.io.gff``; this engine takes plain :class:`GeneRow` records to stay decoupled.

All coordinates are 0-based half-open.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

from privy.synteny.chain import ChainParams, chain_anchors
from privy.synteny.model import Anchor, AnchorSource, GenomeInterval, SyntenyBlock


@dataclass(frozen=True)
class GeneRow:
    """A gene model: id, genome, contig, span, strand."""

    gene_id: str
    genome: str
    contig: str
    start: int
    end: int
    strand: str = "+"


def gene_anchors(
    query_genes: Sequence[GeneRow],
    target_genes: Sequence[GeneRow],
    orthogroups: Mapping[str, str],
) -> list[Anchor]:
    """Build gene-pair anchors between two genomes from a shared orthogroup map.

    For every orthogroup present in both genomes, each query×target gene pair
    becomes an :class:`Anchor` (relative strand ``+`` when the genes share
    orientation, else ``-``).  Genes absent from *orthogroups* are skipped.
    """
    q_by_og: dict[str, list[GeneRow]] = {}
    for gene in query_genes:
        og = orthogroups.get(gene.gene_id)
        if og is not None:
            q_by_og.setdefault(og, []).append(gene)
    t_by_og: dict[str, list[GeneRow]] = {}
    for gene in target_genes:
        og = orthogroups.get(gene.gene_id)
        if og is not None:
            t_by_og.setdefault(og, []).append(gene)

    anchors: list[Anchor] = []
    for og in q_by_og.keys() & t_by_og.keys():
        for q in q_by_og[og]:
            for t in t_by_og[og]:
                strand = "+" if q.strand == t.strand else "-"
                anchors.append(
                    Anchor(
                        query=GenomeInterval(q.genome, q.contig, q.start, q.end),
                        target=GenomeInterval(t.genome, t.contig, t.start, t.end),
                        strand=strand,
                        source=AnchorSource.GENE,
                        name=og,
                    )
                )
    return anchors


def build_gene_synteny(
    query_genes: Sequence[GeneRow],
    target_genes: Sequence[GeneRow],
    orthogroups: Mapping[str, str],
    params: ChainParams | None = None,
) -> list[SyntenyBlock]:
    """Orthogroup-anchored synteny blocks between two genomes (anchors → chain)."""
    return chain_anchors(gene_anchors(query_genes, target_genes, orthogroups), params)
