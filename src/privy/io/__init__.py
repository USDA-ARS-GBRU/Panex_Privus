"""IO subpackage for Panex Privus.

Provides format-specific readers and writers.  All readers produce
format-agnostic domain objects (:class:`~privy.core.locus.Locus`,
:class:`~privy.core.evidence.EvidenceRecord`).  All writers consume
those objects.

Modules:
    vcf     — VCF/BCF streaming reader (pysam/tabix)
    bam     — BAM depth and allele queries (pysam)
    gfa     — GFA graph context reader
    xmfa    — XMFA alignment block reader
    bed     — BED region file reader/writer
    tsv     — TSV output writers (hits, regions, evidence, qc, …)
    jsonio  — JSON run metadata writer/reader
"""
