"""Chromosome-structure context for plant/crop genomes.

Pure-Python k-mer detection of telomeric and centromeric/satellite repeat arrays,
and a deterministic chromosome binning (telomere / arm / pericentromere /
centromere) that can be carried as a context track alongside synteny and
private-allele results.

Operates on a DNA sequence string (e.g. a reference path reconstructed from the
graph, or an assembly contig).  A FASTA-input CLI is intentionally out of scope
for now — Privy's primary inputs are GFA/VCF.
"""

from __future__ import annotations

from privy.structure.karyotype import (
    StructureBin,
    bin_chromosome,
    find_centromere,
    find_telomeres,
    kmer_blocks,
)

__all__ = [
    "StructureBin",
    "bin_chromosome",
    "find_centromere",
    "find_telomeres",
    "kmer_blocks",
]
