"""Microhaplotype (allele-space) layer for Panex Privus.

A *microhaplotype* is a short genomic locus carrying several tightly linked
variants — a multi-allelic marker.  In a pangenome graph it is the allele-space
analog of a syntenic region: synteny says *where* (collinear regions, projectable
to any reference), microhaplotypes say the multi-allelic *what* inside them
(local phased haplotype alleles).  Target-private microhaplotypes are Privy's
existing private-allele core at higher resolution.
"""

from __future__ import annotations

from privy.microhap.detect import detect_microhaplotypes
from privy.microhap.model import Microhaplotype

__all__ = ["Microhaplotype", "detect_microhaplotypes"]
