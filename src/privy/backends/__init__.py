"""Backend subpackage for Panex Privus.

Backends implement the data-source-specific analysis passes.  Each backend
consumes domain objects from :mod:`privy.core` and produces
:class:`~privy.core.evidence.EvidenceRecord` instances.

Backends are deliberately separate from the IO layer:
    - IO modules read and write raw file formats.
    - Backend modules apply biological logic to produce typed evidence.

Modules:
    vcf_scan      — VCF-first private-allele discovery (Phase 2)
    bam_support   — BAM read-level support queries (Phase 3)
    gfa_support   — GFA graph-context annotation (Phase 4)
    xmfa_support  — XMFA alignment corroboration (Phase 5)
"""
