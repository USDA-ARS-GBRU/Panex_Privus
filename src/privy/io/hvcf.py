"""PHG hVCF read/write — interop bridge to the Practical Haplotype Graph ecosystem.

hVCF is VCF v4.2 with **symbolic ALT alleles whose IDs are the MD5 checksum of the
haplotype sequence** (the PHG convention Privy already adopts for microhaplotype
allele ids).  Writing microhaplotypes as hVCF lets PHG / rPHG / BrAPI tools consume
Privy's loci; reading hVCF lets Privy ingest PHG haplotype calls.

This is a pragmatic subset: one record per locus (a reference range), one column per
genome/haplotype path (haploid GT = the 1-based ALT index).  Coordinates are
0-based half-open internally; POS is written 1-based per VCF.

Spec: https://phg.maizegenetics.net/hvcf_specifications/
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path

from privy.microhap.model import Microhaplotype


@dataclass(frozen=True)
class HvcfRecord:
    """One hVCF locus: a reference range with per-genome haplotype (MD5) alleles."""

    locus_id: str
    contig: str
    start: int                       # 0-based
    end: int                         # 0-based exclusive
    alleles: dict[str, str] = field(default_factory=dict)   # genome/path -> MD5 allele id


def write_hvcf(
    loci: Sequence[Microhaplotype],
    path: Path,
    *,
    samples: Sequence[str] | None = None,
) -> Path:
    """Write microhaplotype loci to an hVCF file.

    Args:
        samples: ordered genome/path ids → VCF sample columns (default: union over
            *loci*, first-seen order).
    """
    if samples is None:
        ordered: list[str] = []
        for mh in loci:
            for genome in mh.alleles:
                if genome not in ordered:
                    ordered.append(genome)
        samples = ordered

    all_alleles = sorted({a for mh in loci for a in mh.alleles.values()})
    lines = [
        "##fileformat=VCFv4.2",
        '##INFO=<ID=END,Number=1,Type=Integer,Description="Reference range end (1-based)">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Haplotype (ALT) index">',
    ]
    lines += [f'##ALT=<ID={a},Description="haplotype">' for a in all_alleles]
    lines.append(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples)
    )

    for mh in loci:
        order: list[str] = []
        for genome in samples:
            allele = mh.alleles.get(genome)
            if allele is not None and allele not in order:
                order.append(allele)
        if not order:
            continue
        idx = {a: i + 1 for i, a in enumerate(order)}
        alt = ",".join(f"<{a}>" for a in order)
        gts = [
            "." if (a := mh.alleles.get(genome)) is None else str(idx[a])
            for genome in samples
        ]
        row = [
            mh.contig, str(mh.start + 1), mh.locus_id, "N", alt, ".", "PASS",
            f"END={mh.end}", "GT", *gts,
        ]
        lines.append("\t".join(row))

    out = Path(path)
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out


def read_hvcf(path: Path) -> list[HvcfRecord]:
    """Read an hVCF file into :class:`HvcfRecord` objects (one per locus)."""
    samples: list[str] = []
    records: list[HvcfRecord] = []
    with open(path, encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n\r")
            if not line or line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.split("\t")[9:]
                continue
            fields = line.split("\t")
            if len(fields) < 9:
                continue
            contig = fields[0]
            start = int(fields[1]) - 1
            locus_id = fields[2]
            alt_alleles = [a.strip("<>") for a in fields[4].split(",")]
            info = dict(
                kv.split("=", 1) for kv in fields[7].split(";") if "=" in kv
            )
            end = int(info.get("END", str(start)))
            alleles: dict[str, str] = {}
            for sample, gt in zip(samples, fields[9:], strict=False):
                if gt in (".", ""):
                    continue
                # Haploid GT = ALT index; take the first sub-allele if phased/polyploid.
                first = gt.replace("|", "/").split("/")[0]
                if first == "." or not first.isdigit():
                    continue
                alleles[sample] = alt_alleles[int(first) - 1]
            records.append(
                HvcfRecord(locus_id=locus_id, contig=contig, start=start, end=end, alleles=alleles)
            )
    return records
